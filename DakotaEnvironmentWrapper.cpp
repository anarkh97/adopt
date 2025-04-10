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
  
  if(check_bcast_construct) {
//    done_modifying_db();
    Cerr << "ADAPTIVE SOGA: Problem database has already been "
         << "broadcasted. Delay check and broadcast operations."
         << "\n";
    abort_handler(METHOD_ERROR);
  }


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

  if(check_bcast_construct) {
//    done_modifying_db();
    Cerr << "ADAPTIVE SOGA: Problem database has already been "
         << "broadcasted. Delay check and broadcast operations."
         << "\n";
    abort_handler(METHOD_ERROR);
  }

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
  problem_db.resolve_top_method(); // allow DB set/get operations.

  if(world_rank == 0) {

    //! Get the top (active) iterator and model pointers for restoration
    size_t method_index = problem_db.get_db_method_node();
    size_t model_index = problem_db.get_db_model_node();

    //! AN: We assume that there is only one method and one model
    //!     if user wants to do the adaptive optimization.
    problem_db.set_db_method_node(0);
    problem_db.set_db_model_nodes(0);

    unsigned short method_type = problem_db.get_ushort("method.algorithm");
    const String &model_type = problem_db.get_string("model.type");

    //-------------------------------------------------------------------------
    // Optimization with Single Objective Genetic Algorithm
    //-------------------------------------------------------------------------

    if(method_type == Dakota::SOGA and model_type == "simulation") {

      //! Check for ADAPTIVE OPTIMIZATION
      const StringArray &state_string_labels = 
        problem_db.get_sa("variables.discrete_state_set_string.labels");

      for(const auto &label : state_string_labels) {
        if(label == "SWITCH") {
          adaptive_optimization = true;
          break;
        }
      }

      if(adaptive_optimization) {

        //! Change continuous STATE variable labels for SOFICS
        size_t num_continuous_state = problem_db.get_sizet(
          "variables.continuous_state");

        if(num_continuous_state <= 0) {
          Cerr << "ADAPTIVE SOGA: Number of interpolation points "
               << "not specified. These need to be under continuous "
               << "state variables for adaptive optimization.\n";
          abort_handler(METHOD_ERROR);
        }

        StringArray updated_labels(num_continuous_state);

        //! These descriptors are recognized by SOFICS
        for(size_t i=0; i<num_continuous_state; ++i)
          updated_labels[i] = "NEIGHBOR_" + std::to_string(i+1);
        
        //! Set new labels
        problem_db.set("variables.continuous_state.labels", updated_labels);

      }

    }

    //! Reset pointers back to original
    problem_db.set_db_method_node(method_index);
    problem_db.set_db_model_nodes(model_index);

  }

  //! Check inputs and create Dakota objects
  done_modifying_db();

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

  //---------------------------------------------------------------------------
  // Setup top method iterator
  //---------------------------------------------------------------------------

  if(adaptive_optimization) {
  
    //! Make sure that STATE variables are set to inactive
    Variables &top_vars = top_model->current_variables();
    top_vars.inactive_view(MIXED_STATE);

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

