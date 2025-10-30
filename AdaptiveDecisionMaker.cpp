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

void AdaptiveDecisionMaker::GetEvalTypeAndMetaData(const RealVector &cont_vars, 
                                                   String &into_type,
                                                   IntVector &into_metadata, 
                                                   size_t num_points,
                                                   bool flag)
{

  assert(into_metadata.length() == num_points+1);

  //! First get evaluation type "ERROR", "TRUE", "APPROX"
  GetEvaluationType(cont_vars, into_type, flag);

  //! Avoid neighbor calculation if evaluation is "TRUE"
  if(into_type == "TRUE") {
    for(int i=0; i<into_metadata.length(); ++i) 
      into_metadata[i] = -1;
    return;
  }

  int target = -1; //! target evaluation id for "ERROR" runs.
  IntVector candidates(num_points);
  GetTargetAndNearestNeighbors(cont_vars, target, candidates, num_points, flag);

  //! Fill up neighbors for "ERROR" and "APPROX" types
  into_metadata[0] = (into_type == "APPROX") ? -1 : target;
  for(int i=0; i<num_points; ++i) {
    into_metadata[i+1] = candidates[i];
  }

}

//------------------------------------------------------------------------------

void AdaptiveDecisionMaker::GetTargetAndNearestNeighbors(const RealVector &cont_vars, 
                                                         int &target,
                                                         IntVector &candidates, 
                                                         size_t num_points,
                                                         bool flag)
{

  if(candidates.length() != num_points) {
    candidates.size(num_points);
  }

  //! Throw an error if force find is true but true-db is empty.
  if(flag and id2var.empty()) {
    Cerr << "Adaptive SOGA Error: Trying to find nearest neighbors "
         << "in an empty database.\n";
    abort_handler(METHOD_ERROR);
  }

  //! Throw an error true-db does not have enough points
  if(flag and id2var.size() < num_points) {
    Cerr << "Adaptive SOGA Error: Number of true evaluations stored "
         << "(" << id2var.size() << ") "
         << "is less than number of neighbors requested "
         << "(" << num_points << ").\n";
    abort_handler(METHOD_ERROR);
  }

  // AN: Currently doing naive search.
  // TODO: Could use KDTree.
  
  vector<pair<double, int>> dist2targ;
  IntRealVectorMap::const_iterator tr_it = id2var.begin();
  const IntRealVectorMap::const_iterator tr_e = id2var.end();
  for(; tr_it!=tr_e; ++tr_it) {
    //! NOTE: In asynchronous runs, an evaluation's ID and its associated
    //! variables are registered in our data containers as soon as the job
    //! is queued. Because the simulation might still be running, we skip
    //! any "true-evaluation" IDs that don't yet have an error entry, as
    //! the errors are logged only after the evaluation finishes.
    
    if(!flag and id2error.find(tr_it->first) == id2error.end()) {
      //! This check is bypassed when flag is set to true, which is the
      //! case for error evaluations.
      continue;
    }

    double d = EuclideanDistance(tr_it->second, cont_vars);
    dist2targ.push_back({d, tr_it->first});
  }

  sort(dist2targ.begin(), dist2targ.end());

  // fill up first "num_points" into candidates 
  int c=0;
  for(int i=0; i<dist2targ.size(); ++i) {
    if(dist2targ[i].first == 0) {
      target = dist2targ[i].second;
      continue;
    }

    if(c >= num_points) {
      break;
    }

    candidates[c] = dist2targ[i].second;
    ++c;
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
  if(id2var.empty()) {
    return;
  }
  
  //! Return when GP model is not ready.
  if(!ready_to_predict) {
    return;
  }

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

  VectorXd query_prediction = gp_model.value(query);
  VectorXd query_variance   = gp_model.variance(query);
  
  double standard_dev = std::sqrt(query_variance(0));

  if(3*standard_dev < 6e-5 /*1e-2*/) { // check if 99.7% CI is narrow 
    into = (query_prediction(0) < 5e-5 /*1.3e-2*/) ? "APPROX" : "TRUE";
  }

  if(verbose>1) { 
    //cont_vars.print(Cout) << 
    Cout << "\n";
    Cout << "variables  : " << query << "\n";
    Cout << "uncertainty: " << standard_dev << "\n";
    Cout << "prediction : " << query_prediction(0) << "\n"; 
    Cout << "decision   : " << into << "\n";
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
    if(it != id2var.end()) {
      if(verbose>1) {
        Cout << "Adaptive SOGA Warning: Duplicate evaluation ID " << eval_id
             << ", previously mapped to variables: " << it->second
             << ". Skipping ...\n";
      }
    }
    else {
      id2var[eval_id] = cont_vars;
    }
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
  IntStringMap::iterator it = id2type.find(eval_id);
  if(it != id2type.end() and verbose > 1) {
    Cout << "Adaptive SOGA Warning: Duplicate evaluation ID " << eval_id
         << ", previously mapped to evalutation type: " << it->second
         << ".Skipping ...\n";
    return;
  }

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
  options.sublist("Trend").sublist("Options").set("max degree", 0);
  options.sublist("Trend").sublist("Options").set("p-norm", 2.0);

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
  VectorXd prediction_values = gp_model.value(testing_samples);

  //! calculate loss (root mean squared loss)
  loss = (prediction_values - testing_values).squaredNorm();
  loss = std::sqrt(loss/num_holdout);

  //! calculate standardized loss
  //double test_mean   = testing_values.sum() / num_holdout; 
  //double numerator   = (testing_values - prediction_values).squaredNorm();
  //double denominator = (
  //  testing_values - test_mean * VectorXd::Ones(num_holdout)
  //).squaredNorm(); 

  //loss = std::sqrt(numerator/denominator); // root mean squared.

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
  if(loss < 5e-5/*5e-2*/) { 
    ready_to_predict = true;

    if(verbose>0)
      Cout << "Adaptive SOGA: Gaussian Process Regression model is ready "
           << "to predict.\n";
  }

  WriteGaussianProcessModel();
  WriteCurrentModelResults(parameters, responses, loss);

}

//------------------------------------------------------------------------------

void AdaptiveDecisionMaker::WriteGaussianProcessModel()
{
  //! Temporarily save the model -- will be added as an option later.
  char full_name[256];
  sprintf(full_name, "GaussianModel_%04d", num_train_calls);
  String fname(full_name);
  GaussianProcess::save(gp_model, fname, true);   
}

//------------------------------------------------------------------------------

void AdaptiveDecisionMaker::WriteCurrentModelResults(const MatrixXd &parameters,
                                                     const VectorXd &responses,
                                                     const double loss)
{
  //! Output data for plotting/debugging
  char full_name[256];
  sprintf(full_name, "GaussianModel_%04d_data.txt", num_train_calls);
  FILE *model_data_file = fopen(full_name, "w");
  if(!model_data_file) {
    Cerr << "Adaptive SOGA Error: Could not open the file "
         << "\"" << String(full_name) << "\".\n";
    abort_handler(METHOD_ERROR);
  }

  int num_variables = parameters.cols();
  int num_points    = parameters.rows();

  VectorXd variance = gp_model.variance(parameters);
  VectorXd p_responses = gp_model.value(parameters);

  fprintf(model_data_file, "## Training epoch %04d\n", num_train_calls);
  fprintf(model_data_file, "## Loss %16.8e\n", loss);
  fprintf(model_data_file, "## ");
  for(int i=0; i<num_variables; ++i) {
    fprintf(model_data_file, "Parameter %04d  |  ", i+1);
  }
  fprintf(model_data_file, "NRMSE  |  Prediction  |  Uncertainty\n");

  for(int i=0; i<num_points; ++i) {
    for(int j=0; j<num_variables; ++j) {
      fprintf(model_data_file, "%16.8e  ", parameters(i,j));
    }
    fprintf(model_data_file, "%16.8e  %16.8e  %16.8e\n", 
      responses(i), p_responses(i), std::sqrt(variance(i)));
  }
  fclose(model_data_file);
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

