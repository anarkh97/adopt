#include<DakotaEnvironmentWrapper.h>
#include<AdaptiveJegaOptimizer.h>

#include<DakotaModel.hpp>
#include<DakotaInterface.hpp>

using namespace std;
using namespace Dakota;

//-----------------------------------------------------------------------------

DakotaEnvironmentWrapper::DakotaEnvironmentWrapper(ProgramOptions options,
                                                   bool check_bcast_construct,
                                                   DbCallbackFunctionPtr callback,
                                                   void *callback_data)
                        : LibraryEnvironment(options, check_bcast_construct)
{
  
  SetupAdaptiveOptimizer();

}

//-----------------------------------------------------------------------------

DakotaEnvironmentWrapper::DakotaEnvironmentWrapper(MPI_Comm dakota_comm,
                                                   ProgramOptions options,
                                                   bool check_bcast_construct,
                                                   DbCallbackFunctionPtr callback,
                                                   void *callback_data)
                        : LibraryEnvironment(
                            dakota_comm, options, check_bcast_construct)
{

  SetupAdaptiveOptimizer();

}

//-----------------------------------------------------------------------------

void DakotaEnvironmentWrapper::SetupAdaptiveOptimizer()
{
  
  // reset the model to our model if top method is soga
  String top_method_name = topLevelIterator.method_string();

  // calling it here explicitly to make things clear and readable.
  ProblemDescDB &problem_db = problem_description_db();
  shared_ptr<Model> top_model = topLevelIterator.iterated_model();

  const String& top_model_type = top_model->model_type();

  if(top_method_name == "soga" and top_model_type == "simulation") {

    // ProblemDescDB::set does not allow tinkering of analysis
    // drivers. So, we simply create a new Model with the same
    // interface. We handle different kinds of evaluations 
    // in the driver through STATE variables.
    // check if a string set of state variables was defined.
/*
    StringSetArray string_sets = ModelUtils::discrete_set_string_values(
      *top_model, Dakota::MIXED_STATE);

    bool perform_adaptive_optimization = false;

    if(!string_sets.empty()) {

      // We are looking for a string set whose permissible values are
      // "TRUE", "APPROX", and "ERROR".
      int num_matches = 0;
      for(auto &sets : string_sets) {
        if(sets.size() > 3) continue;

	for(auto &val : sets)
          if(val == "TRUE" or val == "APPROX" or val == "ERROR")
            num_matches++;
      }

      // we found all three keyword, so we assume that the user
      // wants to perform our adaptive optimization.
      perform_adaptive_optimization = true;

    }
*/

    bool perform_adaptive_optimization = false;
    StringMultiArrayConstView disc_string_labels =
      ModelUtils::all_discrete_string_variable_labels(
        *top_model);

    for(const auto &label : disc_string_labels) {
      if(label == "SWITCH") 
        perform_adaptive_optimization = true;
    }

    if(perform_adaptive_optimization) {

      // get the communicator attached to environment
      ParallelLibrary &parallel_lib = parallel_library();
      ParLevLIter w_pl_iter = parallel_lib.w_parallel_level_iterator();

      Cout << "\n";
      Cout << "----------------------------\n";
      Cout << "--        SOGA-Ad.        --\n";
      Cout << "----------------------------\n";

      shared_ptr<Model> err_model =
        ModelUtils::get_model(problem_db);

      // assign our optimizer.
      topLevelIterator.assign_rep(
        make_shared<AdaptiveJegaOptimizer>(problem_db, top_model, err_model));
      topLevelIterator.init_communicators(w_pl_iter);

    }
  }

  abort_handler(OTHER_ERROR);

}

//-----------------------------------------------------------------------------

