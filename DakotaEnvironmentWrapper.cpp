#include<DakotaEnvironmentWrapper.h>
#include<AdaptiveJegaOptimizer.h>
#include<ForkApplicInterfaceWrapper.h>

#include<DakotaModel.hpp>
#include<DakotaInterface.hpp>

using namespace std;
using namespace Dakota;

//-----------------------------------------------------------------------------

DakotaEnvironmentWrapper::
DakotaEnvironmentWrapper(ProgramOptions options,
                         bool check_bcast_construct,
                         DbCallbackFunctionPtr callback,
                         void *callback_data)
: LibraryEnvironment(options, check_bcast_construct),
  adaptive_optimization(false)
{
  
  if(!check_bcast_construct)
    done_modifying_db();

  SetupAdaptiveOptimizationVariables();
  SetupAdaptiveOptimizer();

  if(adaptive_optimization) {
    Cout << "\n";
    Cout << "----------------------------\n";
    Cout << "--        SOGA-Ad.        --\n";
    Cout << "----------------------------\n";
  }
  
}

//-----------------------------------------------------------------------------

DakotaEnvironmentWrapper::
DakotaEnvironmentWrapper(MPI_Comm dakota_comm,
                         ProgramOptions options,
                         bool check_bcast_construct,
                         DbCallbackFunctionPtr callback,
                         void *callback_data)
: LibraryEnvironment(dakota_comm, options, check_bcast_construct),
  adaptive_optimization(false)
{

  if(!check_bcast_construct)
    done_modifying_db();

  SetupAdaptiveOptimizationVariables();
  SetupAdaptiveOptimizer();

  if(adaptive_optimization) {
    Cout << "\n";
    Cout << "----------------------------\n";
    Cout << "--        SOGA-Ad.        --\n";
    Cout << "----------------------------\n";
  }
  
}

//-----------------------------------------------------------------------------

void DakotaEnvironmentWrapper::SetupAdaptiveOptimizationVariables()
{

  //! Get world rank from Dakota
  ParallelLibrary &parallel_lib = parallel_library();
  int world_rank = parallel_lib.world_rank();

  //! Here we get the user inputs for discrete string state set variables.
  //! If the user specified the keyword "SWITCH" then we infer that the
  //! user wants to perform our adaptive optimization. We also infer that
  //! number of continuous state variables is the number of points for
  //! interpolation.
  
  ProblemDescDB &problem_db = problem_description_db();
  shared_ptr<Model> top_model = topLevelIterator.iterated_model();

  //! Get method and model information
  const String &top_method_name = topLevelIterator.method_string();
  const String &top_model_type = top_model->model_type();

  //! handle trivial case
  if(top_method_name != "soga" and top_model_type != "recast")
    return;

  //---------------------------------------------------------------------------
  // Optimization with Single Objective Genetic Algorithm
  //---------------------------------------------------------------------------

  //! Correct the inactive view to STATE
  Variables &top_model_vars = top_model->current_variables();
  top_model_vars.inactive_view(MIXED_STATE);

  //! Check for ADAPTIVE OPTIMIZATION
  StringMultiArrayConstView state_string_labels = 
    ModelUtils::inactive_discrete_string_variable_labels(*top_model);

  for(const auto &label : state_string_labels) {
    if(label == "SWITCH") {
      adaptive_optimization = true;
      break;
    }
  }

  //! Change continuous STATE variable labels for SOFICS
  size_t num_icv = ModelUtils::icv(*top_model);
  StringMultiArray new_state_var_labels(boost::extents[num_icv]);
  for(size_t i=0; i<num_icv; ++i)
    new_state_var_labels[i] = "NEIGHBOR_" + std::to_string(i+1);
  
  StringMultiArrayConstView state_var_view = 
    new_state_var_labels[boost::indices[idx_range(0,num_icv)]];
  ModelUtils::inactive_continuous_variable_labels(*top_model, 
    state_var_view);

}

//-----------------------------------------------------------------------------

void DakotaEnvironmentWrapper::SetupAdaptiveOptimizer()
{
  
  //! Get world rank from Dakota
  ParallelLibrary &parallel_lib = parallel_library();
  int world_rank = parallel_lib.world_rank();

  //! Get problem description database and iterated model 
  //! attached to the top method pointer.
  ProblemDescDB &problem_db = problem_description_db();
  shared_ptr<Model> top_model = topLevelIterator.iterated_model();

  //----------------------------------------------------------------------------
  // Setup top method iterator
  //---------------------------------------------------------------------------

  if(adaptive_optimization) {
  
    //! Reference for reseting communicator
    ParLevLIter w_pl_iter = parallel_lib.w_parallel_level_iterator();
  
    //! Create a new Dakota Model for error evaluations.
    shared_ptr<Model> err_model = ModelUtils::get_model(problem_db);
  
    //! Reset to our own interface. This is done to distinguish
    //! working directories of true models and error models.
    Interface &err_model_interface = err_model->derived_interface();   
    err_model_interface.assign_rep(
      std::make_shared<ForkApplicInterfaceWrapper>(problem_db));
  
    //! Assign our optimizer.
    topLevelIterator.assign_rep(
      make_shared<AdaptiveJegaOptimizer>(problem_db, top_model, err_model));

    //! Reset communicator
    topLevelIterator.init_communicators(w_pl_iter);
  
  }

  //abort_handler(OTHER_ERROR);

}

//-----------------------------------------------------------------------------

