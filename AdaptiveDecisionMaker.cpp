#include<algorithm>
#include<AdaptiveDecisionMaker.h>

// Dakota header for read/write
// Note: This header also has functions for 
// read/write from/to file.
#include<dakota_data_io.hpp>

// dakota util include
#include<util/util_math_tools.hpp>

using namespace std;
using namespace Dakota;
using namespace Eigen;
using namespace dakota;
using namespace dakota::util;
using namespace dakota::surrogates;

//------------------------------------------------------------------------------
//! Helper function
double EuclideanDistance(const RealVector &v1, const RealVector &v2)
{
  assert(v1.length() == v2.length());
  double dist=0;
  int dim=v1.length();

  for(int i=0; i<dim; ++i)
    dist += pow(v1[i]-v2[i],2);

  return sqrt(dist);

}

//------------------------------------------------------------------------------

AdaptiveDecisionMaker::AdaptiveDecisionMaker()
                     : true_evals(), approx_evals(), error_vals(),
                       readyToPredict(false)
{
  // empty ctor
}

//------------------------------------------------------------------------------

void AdaptiveDecisionMaker::GetNearestNeighbors(const RealVector &cont_vars, 
                                                IntVector &into, 
                                                size_t num_neighbors,
                                                bool force)
{

  if(into.length() != num_neighbors) 
    into.size(num_neighbors);

  // AN: Throw an error if force find is true but true-db is empty.
  // This distinction is made so that when decision maker is called
  // for error calculations one can check if database with true 
  // evaluations is empty.
  if(true_evals.empty() and force) {
    Cerr << "ADAPTIVE SOGA: Trying to find nearest neighbors "
         << "in an empty database.\n";
    abort_handler(METHOD_ERROR);
  }

  // Send garbage values when true-db is empty.
  if(true_evals.empty()) {
    for(int i=0; i<num_neighbors; ++i) into[i] = -1;
    return;
  }

  if(true_evals.size() < num_neighbors) {
    Cerr << "ADAPTIVE SOGA: Number of true evaluations stored "
         << "(" << true_evals.size() << ") "
         << "is less than number of neighbors requested "
         << "(" << num_neighbors << ").\n";
    abort_handler(METHOD_ERROR);
  }

  // AN: Currently doing naive search.
  // TODO: Could use KDTree.
  
  vector<pair<double, int>> dist2targ;
  IntRealVectorMap::const_iterator tr_it = true_evals.begin();
  const IntRealVectorMap::const_iterator tr_e = true_evals.end();
  for(; tr_it!=tr_e; ++tr_it) {
    double d = EuclideanDistance(tr_it->second, cont_vars);
    dist2targ.push_back({d, tr_it->first});
  }

  sort(dist2targ.begin(), dist2targ.end());

  // fill up first "num_neighbors" into IntVector
  for(int i=0; i<num_neighbors; ++i) {
    into[i] = dist2targ[i].second;
  }

}

//------------------------------------------------------------------------------

bool
AdaptiveDecisionMaker::GetEvaluationDecision(const RealVector &cont_vars)
{
  
  // Send "TRUE" decision when true-db is empty.
  if(true_evals.empty()) return true;

  // predict only when model is switched on.
  if(readyToPredict) {
    
    MatrixXd query(1, cont_vars.length());
    VectorXd var_query = gp_model.variance(query);

    // predict the error when the model is certain enough
    if(var_query(0) < 1e-2) {
      VectorXd pred_query = gp_model.value(query);
      return (pred_query(0) < 1e-2);
    }
    else 
      return false;

  }

  // if we reached here then the model was 
  // not used to make predictions.
  return false;

}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::RecordEvaluationError(const int eval_id,
                                             const RealVector &cont_vars,
                                             const double &error)
{

  const RealVector &stored_vars = true_evals[eval_id];

  // Check if variables match.  This is done to ensure error values stored
  // in error_vals are mapped to expected continuous variables stored in 
  // true_evals.
  if(cont_vars != stored_vars) {
    Cerr << "ADAPTIVE SOGA: Trying to map error values for unknown design "
         << "variables.\n";
    abort_handler(METHOD_ERROR);
  }

  error_vals[eval_id] = error;

}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::RecordEvaluationDecision(int eval_id, 
                                                const RealVector &cont_vars, 
                                                const bool eval_type)
{

  if(eval_type == true) {
    // First check if evaluation id has already been mapped or not.
    IntRealVectorMap::iterator it = true_evals.find(eval_id);
    if(it != true_evals.end()) {
      Cerr << "ADAPTIVE SOGA: Duplicate evaluation ID " << eval_id
           << ", previously mapped to variables: " << it->second;
      abort_handler(METHOD_ERROR);
    }

    true_evals[eval_id] = cont_vars;
  }
  else {
    // Firt check if evaluation id has already been mapped or not.
    IntRealVectorMap::iterator it = approx_evals.find(eval_id);
    if(it != approx_evals.end()) {
      Cerr << "ADAPTIVE SOGA: Duplicate evaluation ID " << eval_id
           << ", previously mapped to variables: " << it->second;
      abort_handler(METHOD_ERROR);
    }

    approx_evals[eval_id] = cont_vars;
  }

}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::LoadGaussianProcesssOptions()
{

  ParameterList options;
  gp_model.get_options(options);

  options.set("verbosity", 1);
  // rng options.
  options.set("gp seed", 42);
  // scaling options
  options.set("scaler name", "none");
  options.set("standardize response", false);

  // kernel options
  options.set("kernel type", "squared exponential");

  // polynomial trend options.
  options.sublist("Trend").sublist("estimate trend", true);

  // nugget options.
  options.sublist("Nugget").set("estimate nugget", true);
  options.sublist("Nugget").sublist("Bounds").set("lower bound", 0);
  options.sublist("Nugget").sublist("Bounds").set("upper bound", 1);

  // hyperparameter optimization options.
  options.set("num restarts", 20);

  gp_model.set_options(options);

}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::LoadParameters(const IntRealVectorMap &from,
                                      MatrixXd &into) const
{

  int var_dim = from.begin()->second.length();

  // prepare for GP model.
  // Dakota GP model wants parameters in an Eigen::MatrixXd
  // with rows as the number of points and columns assert
  // number of dimensions (or features).
  
  // allocate size to MatrixXd
  into.resize(from.size(), var_dim);

  IntRealVectorMap::const_iterator it = from.begin();
  const IntRealVectorMap::const_iterator it_e = from.end();

  int counter = 0;
  for(; it!=it_e; it++) {

    const RealVector &from_var = it->second;
    for(int d=0; d<var_dim; ++d) 
      into(counter, d) = from_var[d];
    counter++;

  }

}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::LoadResponses(const IntRealMap &from,
                                     VectorXd &into) const
{

  // prepare for GP model.
  // Dakota GP model wants responses (i.e., true or observed values)
  // in an Eigen::VectorXd with rows as the number of points.

  // allocate size to VectorXd
  into.resize(from.size());
  
  IntRealMap::const_iterator it = from.begin();
  const IntRealMap::const_iterator it_e = from.end();

  int counter = 0;
  for(; it!=it_e; it++) {
    into(counter) = it->second;
    counter++;
  }

}

//------------------------------------------------------------------------------

double
AdaptiveDecisionMaker::BuildGaussianProcessModel(const MatrixXd &samples,
                                                 const VectorXd &values,
                                                 const double split_ratio,
                                                 const size_t seed)
{

  assert(samples.rows() == values.rows());

  if(split_ratio == 0)  {
    Cerr << "ADAPTIVE SOGA: Train/Test split ratio not provided.\n";
    abort_handler(METHOD_ERROR);
  }

  // split data set into training and testing sets.
  int var_dim     = samples.cols();
  int num_samples = samples.rows();
  int num_holdout = (int)(num_samples*split_ratio);
  int num_train   = num_samples - num_holdout;

  MatrixXd training_samples(num_train, var_dim);
  VectorXd training_values(num_train);
  MatrixXd testing_samples(num_holdout, var_dim);
  VectorXd testing_values(num_holdout);

  VectorXi index_permutations = 
    VectorXi::LinSpaced(num_samples, 0, num_samples-1);

  // randomize indices (using dakota's random_permutation)
  if(seed == 0) 
    random_permutation(num_samples, (unsigned int)std::time(0), 
                       index_permutations);
  else
    random_permutation(num_samples, (unsigned int)seed, 
                       index_permutations);

  int train_index = 0;
  int test_index = 0;
  for(int i=0; i<num_samples; ++i) {

    // first "num_train" samples are used for training.
    if(i < num_train) {
      int samples_index = index_permutations(i);
      training_samples.row(train_index) = samples.row(samples_index);
      training_values(train_index, 0) = values(samples_index, 0);
      train_index++;
    }
    // next populate test set.
    else {
      int samples_index = index_permutations(i);
      testing_samples.row(test_index) = samples.row(samples_index);
      testing_values(test_index, 0) = values(samples_index, 0);
      test_index++;
    }
  }
  
  // build the model
  gp_model.build(training_samples, training_values);
  
  // test the model
  VectorXd prediction_values = 
    gp_model.value(testing_samples);

  // calculate loss 
  double loss = (prediction_values - testing_values).norm();
  loss /= num_holdout;

  return std::sqrt(loss);

}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::CrossValidateGausssianModel() const
{

  // cross-validation
  // Eigen::VectorXd cross_val_metrics = gp_model.cross_validate(...);

}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::Train()
{

  // NOTE: since we use the default GP constructor, the default
  // GP options have already been setup. Here we override a few
  // GP options based on our requirements.
  LoadGaussianProcesssOptions();

  // prepare to train 
  MatrixXd parameters; // continuous design variables
  VectorXd responses; // error values for the designs

  LoadParameters(true_evals, parameters);
  LoadResponses(error_vals, responses);

  double loss = BuildGaussianProcessModel(parameters, responses);

  // switch the model on once loss is below a threshold.
  if(loss < 1e-2) 
    readyToPredict = true;

}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

