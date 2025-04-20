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

AdaptiveDecisionMaker::AdaptiveDecisionMaker(const ProblemDescDB &problem_db_)
                     : problem_db(problem_db_), true_evals(), approx_evals(), 
                       error_vals(), ready_to_predict(false)
{

  short dakota_level = problem_db.get_short("method.output");
  //! AN: Note these output levels correspond to dakota::surrogates.
  //!     We follow the same convention.
  switch(dakota_level) {
    case SILENT_OUTPUT:
    case QUIET_OUTPUT:
      verbose = 0; /* no output */
      break;
    case NORMAL_OUTPUT:
      verbose = 1; /* normal output */
      break;
    case VERBOSE_OUTPUT:
    case DEBUG_OUTPUT:
      verbose = 2; /* maximum output */
      break;
    default:
      verbose = 0;
      break;
  }

  String options_file = problem_db.get_string("method.advanced_options_file");

  //! Read Gaussian Process options from user specified file.
  ReadOptionsFile(options_file);

  //! NOTE: since we use the default GP constructor, the default
  //! GP options have already been setup. Here we override a few
  //! GP options based on our requirements.
  LoadGaussianProcesssOptions();

}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::ReadOptionsFile(const String filename)
{

  if(filename.empty() and verbose > 0) {
    Cout << "Adaptive SOGA Warning: No options file for " 
         << "the decision maker was provided. Default values "
         << "will be used.\n";
  }

  //! Start parsing options for the decision maker

}

//------------------------------------------------------------------------------

void AdaptiveDecisionMaker::GetNearestNeighbors(const RealVector &cont_vars, 
                                                IntVector &into, 
                                                size_t num_int_points,
                                                bool force)
{

  if(into.length() != num_int_points) 
    into.size(num_int_points);

  //! One point is reserved for target evaluation id.
  int num_interp_points = num_int_points-1;

  // AN: Throw an error if force find is true but true-db is empty.
  // This distinction is made so that when decision maker is called
  // for error calculations one can check if database with true 
  // evaluations is empty.
  if(true_evals.empty() and force) {
    Cerr << "Adaptive SOGA Error: Trying to find nearest neighbors "
         << "in an empty database.\n";
    abort_handler(METHOD_ERROR);
  }

  // Send garbage values when true-db is empty.
  if(true_evals.empty() or true_evals.size() < num_interp_points) {
    for(int i=0; i<num_int_points; ++i) into[i] = -1;
    return;
  }

  //if(true_evals.size() < num_interp_points) {
  //  Cerr << "Adaptive SOGA Error: Number of true evaluations stored "
  //       << "(" << true_evals.size() << ") "
  //       << "is less than number of neighbors requested "
  //       << "(" << num_interp_points << ").\n";
  //  abort_handler(METHOD_ERROR);
  //}

  // AN: Currently doing naive search.
  // TODO: Could use KDTree.
  
  vector<pair<double, int>> dist2targ;
  IntRealVectorMap::const_iterator tr_it = true_evals.begin();
  const IntRealVectorMap::const_iterator tr_e = true_evals.end();
  for(; tr_it!=tr_e; ++tr_it) {
    double d = EuclideanDistance(tr_it->second, cont_vars);
    //if(d == 0) continue; // skip the point itself.
    dist2targ.push_back({d, tr_it->first});
  }

  sort(dist2targ.begin(), dist2targ.end());

  // fill up first "num_int_points" into IntVector
  // Note: first point will be the current point
  // itself (i.e., target)
  for(int i=0; i<num_int_points; ++i) {
    into[i] = dist2targ[i].second;
  }

}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::GetEvaluationType(const RealVector &cont_vars,
                                         String &into,
                                         const bool error_sim)
{

  //! These are simulations for computing errors
  //! for true design evaluations.
  if(error_sim) {
    into = "ERROR";
    return;
  }

  //! These are actual simulations (high- or low-fidelity)
  GetEvaluationDecision(cont_vars, into); 

}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::GetEvaluationDecision(const RealVector &cont_vars,
                                             String &into)
{
  
  //! Set the base case
  into = "TRUE";

  // Return when databse is empty.
  if(true_evals.empty()) {
    return;
  }

  if(ready_to_predict) {
    
    MatrixXd query(1, cont_vars.length());
    VectorXd var_query = gp_model.variance(query);

    // predict the error when the model is certain enough
    if(var_query(0) < 1e-2) {
      VectorXd pred_query = gp_model.value(query);
      into = (pred_query(0) < 1e-2) ? "APPROX" : "TRUE";
    }

  }

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
    Cerr << "Adaptive SOGA Error: Trying to map error values for unknown design "
         << "variables.\n";
    abort_handler(METHOD_ERROR);
  }

  error_vals[eval_id] = error;

}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::RecordEvaluationDecision(int eval_id, 
                                                const RealVector &cont_vars, 
                                                const String &eval_type)
{

  if(eval_type == "TRUE") {
    // First check if evaluation id has already been mapped or not.
    IntRealVectorMap::iterator it = true_evals.find(eval_id);
    if(it != true_evals.end() and verbose > 1) {
      Cout << "Adaptive SOGA Error: Duplicate evaluation ID " << eval_id
           << ", previously mapped to variables: " << it->second
           << ". Skipping ...\n";
    }

    true_evals[eval_id] = cont_vars;
  }
  else if (eval_type == "APPROX") {
    // Firt check if evaluation id has already been mapped or not.
    IntRealVectorMap::iterator it = approx_evals.find(eval_id);
    if(it != approx_evals.end() and verbose > 1) {
      Cout << "Adaptive SOGA Error: Duplicate evaluation ID " << eval_id
           << ", previously mapped to variables: " << it->second
           << ". Skipping ...\n";
    }

    approx_evals[eval_id] = cont_vars;
  }
  else {
    Cerr << "Adaptive SOGA Error: I do not understand evaluation type "
         << eval_type << ".\n";
    abort_handler(METHOD_ERROR);
  }

}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::LoadGaussianProcesssOptions()
{

  ParameterList options;
  gp_model.get_options(options);

  options.set("verbosity", verbose);

  // rng options.
  options.set("gp seed", 42);
  // scaling options
  options.set("scaler name", "none");
  options.set("standardize response", false);

  // kernel options
  options.set("kernel type", "squared exponential");

  // polynomial trend options.
  options.sublist("Trend").set("estimate trend", true);

  // nugget options.
  options.sublist("Nugget").set("estimate nugget", true);
  options.sublist("Nugget").sublist("Bounds").set("lower bound", 0.0);
  options.sublist("Nugget").sublist("Bounds").set("upper bound", 1.0);

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
                                                 double split_ratio,
                                                 const size_t seed)
{

  assert(samples.rows() == values.rows());

  if(split_ratio == 0 and verbose > 0)  {
    Cout << "Adaptive SOGA Warning: Train/Test split ratio not provided. "
         << "Using default (20%).\n";
    split_ratio = 0.2;
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

  // prepare to train 
  MatrixXd parameters; // continuous design variables
  VectorXd responses; // error values for the designs

  LoadParameters(true_evals, parameters);
  LoadResponses(error_vals, responses);

  double loss = BuildGaussianProcessModel(parameters, responses);

  if(verbose>1) {
    Cout << "Adaptive SOGA: Computed loss for Gaussian Process Regression "
         << "model is " << loss*100 << "%.\n";
  }

  // switch the model on once loss is below a threshold.
  if(loss < 5e-2) { 
    if(verbose>0)
      Cout << "Adaptive SOGA: Gaussian Process Regression model is ready "
           << "to predict.\n";
    ready_to_predict = true;
  }

}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

