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
                     : problem_db(problem_db_), id2var(), id2type(), 
                       id2error(), ready_to_predict(false), num_train_calls(0)
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

void AdaptiveDecisionMaker::GetEvaluationAndNeighbors(const RealVector &cont_vars, 
                                                      String &into_type,
                                                      IntVector &into_neighbors, 
                                                      size_t num_points,
                                                      bool flag)
{

  //! First get evaluation type "ERROR", "TRUE", "APPROX"
  GetEvaluationType(cont_vars, into_type, flag);

  //! Avoid neighbor calculation if evaluation is "TRUE"
  if(into_type == "TRUE") {
    for(int i=0; i<num_points; ++i) 
      into_neighbors[i] = -1;
    return;
  }

  //! Fill up neighbors for "ERROR" and "APPROX" types
  GetNearestNeighbors(cont_vars, into_neighbors, num_points);

}

//------------------------------------------------------------------------------

void AdaptiveDecisionMaker::GetNearestNeighbors(const RealVector &cont_vars, 
                                                IntVector &into, 
                                                size_t num_points)
{

  if(into.length() != num_points) 
    into.size(num_points);

  //! One point is reserved for target evaluation id.
  int num_interp_points = num_points-1;

  //! Throw an error if force find is true but true-db is empty.
  if(id2var.empty()) {
    Cerr << "Adaptive SOGA Error: Trying to find nearest neighbors "
         << "in an empty database.\n";
    abort_handler(METHOD_ERROR);
  }

  //! Throw an error true-db does not have enough points
  if(id2var.size() < num_interp_points) {
    Cerr << "Adaptive SOGA Error: Number of true evaluations stored "
         << "(" << id2var.size() << ") "
         << "is less than number of neighbors requested "
         << "(" << num_interp_points << ").\n";
    abort_handler(METHOD_ERROR);
  }

  // AN: Currently doing naive search.
  // TODO: Could use KDTree.
  
  vector<pair<double, int>> dist2targ;
  IntRealVectorMap::const_iterator tr_it = id2var.begin();
  const IntRealVectorMap::const_iterator tr_e = id2var.end();
  for(; tr_it!=tr_e; ++tr_it) {
    double d = EuclideanDistance(tr_it->second, cont_vars);
    //if(d == 0) continue; // skip the point itself.
    dist2targ.push_back({d, tr_it->first});
  }

  sort(dist2targ.begin(), dist2targ.end());

  // fill up first "num_points" into IntVector
  // Note: first point will be the current point
  // itself (i.e., target)
  for(int i=0; i<num_points; ++i) {
    into[i] = dist2targ[i].second;
  }

}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::GetEvaluationType(const RealVector &cont_vars,
                                         String &into,
                                         bool flag)
{

  //! These are simulations for computing errors
  //! for true design evaluations.
  if(flag) {
    into = "ERROR";
    return;
  }

  //! Set the base case
  into = "TRUE";

  //! Return when databse is empty.
  if(id2var.empty()) 
    return;
  
  //! Return when GP model is not ready.
  if(!ready_to_predict)
    return;

  //! Check if there is an existing decision for this variable
  IntRealVectorMap::const_iterator it = id2var.begin();
  const IntRealVectorMap::const_iterator e = id2var.end();
  for(; it!=e; ++it) {
    if(it->second == cont_vars) {
      into = id2type[it->first];
      return;
    }
  }

  //! Use GP model to predict the evaluation decision.
  MatrixXd query(1, cont_vars.length());
  for(int i=0; i<cont_vars.length(); ++i)
    query(0,i) = cont_vars[i];

  VectorXd query_variance   = gp_model.variance(query);
  VectorXd query_prediction = gp_model.value(query);

  if(std::abs(3*query_variance(0)) < 1e-2) 
    into = (query_prediction(0) < 1e-2) ? "APPROX" : "TRUE";

  if(verbose>1) { 
    //cont_vars.print(Cout) << 
    Cout << "\n";
    Cout << "variables : " << query << "\n";
    Cout << "variance  : " << query_variance(0) << "\n";
    Cout << "prediction: " << query_prediction(0) << "\n"; 
    Cout << "decision  : " << into << "\n";
    Cout << "\n";
  }

}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::RecordEvaluationError(const int eval_id,
                                             const RealVector &cont_vars,
                                             const double &error)
{

  const RealVector &stored_vars = id2var[eval_id];

  //! Check if variables match.  This is done to ensure error values stored
  //! in id2error are mapped to expected continuous variables stored in 
  //! id2var.
  if(cont_vars != stored_vars) {
    Cerr << "Adaptive SOGA Error: Trying to map error values "
         << "for unknown design variables.\n";
    abort_handler(METHOD_ERROR);
  }

  id2error[eval_id] = error;

}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::RecordEvaluationDecision(int eval_id, 
                                                const RealVector &cont_vars, 
                                                const String &eval_type)
{

  if(eval_type == "TRUE") {
    //! First check if evaluation id has already been mapped or not.
    IntRealVectorMap::iterator it = id2var.find(eval_id);
    if(it != id2var.end() and verbose>1) {
      Cout << "Adaptive SOGA Warning: Duplicate evaluation ID " << eval_id
           << ", previously mapped to variables: " << it->second
           << ". Skipping ...\n";
    }

    id2var[eval_id] = cont_vars;
  }
  else if(eval_type == "ERROR") {
    Cout << "Adaptive SOGA Warning: Error evaluations are not stored in "
         << "the decision maker. Ignoring...\n";
    return;
  }
  else if(eval_type != "APPROX") {
    Cerr << "Adaptive SOGA Error: I do not understand evaluation type "
         << eval_type << ".\n";
    abort_handler(METHOD_ERROR);
  }

  //! Store the evaluation ids for only for "TRUE" and "APPROX" types
  id2type[eval_id] = eval_type;

}

//------------------------------------------------------------------------------

bool
AdaptiveDecisionMaker::NeedToComputeErrors()
{
  bool ret = id2var.size() != id2error.size();
  return ret;
}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::LoadGaussianProcesssOptions()
{

  ParameterList options;
  gp_model.get_options(options);

  options.set("verbosity", verbose);

  //! rng options.
  options.set("gp seed", 42);
  //! scaling options
  options.set("scaler name", "none");
  options.set("standardize response", false);

  //! kernel options
  options.set("kernel type", "squared exponential");

  //! polynomial trend options.
  options.sublist("Trend").set("estimate trend", true);

  //! nugget options.
  options.sublist("Nugget").set("estimate nugget", true);
  options.sublist("Nugget").sublist("Bounds").set("lower bound", 0.0);
  options.sublist("Nugget").sublist("Bounds").set("upper bound", 1.0);

  //! hyperparameter optimization options.
  options.set("num restarts", 20);

  gp_model.set_options(options);

}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::LoadParameters(const IntRealVectorMap &from,
                                      MatrixXd &into) const
{

  int var_dim = from.begin()->second.length();

  //! Dakota GP model wants parameters in an Eigen::MatrixXd
  //! with rows as the number of points and columns assert
  //! number of dimensions (or features).
  
  //! allocate size to MatrixXd
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

  //! Dakota GP model wants responses (i.e., true or observed values)
  //! in an Eigen::VectorXd with rows as the number of points.

  //! allocate size to VectorXd
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

bool
AdaptiveDecisionMaker::BuildGaussianProcessModel(const MatrixXd &samples,
                                                 const VectorXd &values,
                                                 double &loss,
                                                 double split_ratio,
                                                 const size_t seed)
{

  assert(samples.rows() == values.rows());

  if(split_ratio == 0 and verbose>0)  {
    Cout << "Adaptive SOGA Warning: Train/Test split ratio not provided. "
         << "Using default (20%).\n";
    split_ratio = 0.2;
  }

  //! split data set into training and testing sets.
  int var_dim     = samples.cols();
  int num_samples = samples.rows();
  int num_holdout = (int)(num_samples*split_ratio);
  int num_train   = num_samples - num_holdout;

  if(num_holdout == 0) {
    Cout << "Adaptive SOGA Warning: Insufficient data to create separate "
         << "training and test sets. Waiting till more samples are collected.\n";
    loss = 1e10; /* sending back a random large value. */
    return false;
  }

  MatrixXd training_samples(num_train, var_dim);
  VectorXd training_values(num_train);
  MatrixXd testing_samples(num_holdout, var_dim);
  VectorXd testing_values(num_holdout);

  VectorXi index_permutations = 
    VectorXi::LinSpaced(num_samples, 0, num_samples-1);

  //! randomize indices (using dakota's random_permutation)
  if(seed == 0) 
    random_permutation(num_samples, (unsigned int)std::time(0), 
                       index_permutations);
  else
    random_permutation(num_samples, (unsigned int)seed, 
                       index_permutations);

  int train_index = 0;
  int test_index = 0;
  for(int i=0; i<num_samples; ++i) {

    //! first "num_train" samples are used for training.
    if(i < num_train) {
      int samples_index = index_permutations(i);
      training_samples.row(train_index) = samples.row(samples_index);
      training_values(train_index, 0) = values(samples_index, 0);
      train_index++;
    }
    //! next populate test set.
    else {
      int samples_index = index_permutations(i);
      testing_samples.row(test_index) = samples.row(samples_index);
      testing_values(test_index, 0) = values(samples_index, 0);
      test_index++;
    }
  }
  
  //! build the model
  gp_model.build(training_samples, training_values);
  
  //! test the model
  VectorXd prediction_values = 
    gp_model.value(testing_samples);

  //! calculate loss 
  loss = (prediction_values - testing_values).norm();
  loss = std::sqrt(loss/num_holdout);

  return true;

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

  ++num_train_calls;

  //! Prepare to train 
  MatrixXd parameters; // continuous design variables
  VectorXd responses; // error values for the designs

  LoadParameters(id2var, parameters);
  LoadResponses(id2error, responses);

  //! Return when the database is too small (i.e., build failed).
  double loss = 0.0;
  if(!BuildGaussianProcessModel(parameters, responses, loss)) {
    ready_to_predict = false;
    return;
  }

  //! Switch the GP model on once loss is below a threshold.
  if(loss < 5e-2) { 
    ready_to_predict = true;

    if(verbose>0)
      Cout << "Adaptive SOGA: Gaussian Process Regression model is ready "
           << "to predict.\n";
  }

  //! Debug output
  if(verbose>1) {
    Cout << "Training " << num_train_calls 
         << " (loss = " << loss << ")" << ":\n";
    VectorXd variance = gp_model.variance(parameters);
    VectorXd p_responses = gp_model.value(parameters);

    int num_variables = parameters.cols();
    int num_points = parameters.rows();
    for(int i=0; i<num_points; ++i) {
      for(int j=0; j<num_variables; ++j)
        Cout << parameters(i,j) << "    ";
      Cout << responses(i) << "    "
           << p_responses(i) << "    "
           << variance(i) << "\n";
    }
    Cout << "\n";
  }

}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

