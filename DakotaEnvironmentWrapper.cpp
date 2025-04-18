#include<DakotaEnvironmentWrapper.h>
#include<AdaptiveJegaOptimizer.h>

#include<DakotaModel.hpp>
#include<SimulationModel.hpp>
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

      adaptive_optimization = false;

      //! Find the SWITCH discrete string set state variable
      for(const auto &label : state_string_labels) {
        if(label == "SWITCH") {
          adaptive_optimization = true;
          break;
        }
      }

      if(adaptive_optimization) {

        //---------------------------------------------------------------------
        // Additional data containers for ERROR model.
        //---------------------------------------------------------------------
        DataModel     data_model;
        DataInterface data_interf;
        DataResponses data_respns;

        //---------------------------------------------------------------------
        // Update the True simulation model
        //---------------------------------------------------------------------

        //! Change discrete STATE range variable labels for SOFICS
        size_t num_interp_points = problem_db.get_sizet(
          "variables.discrete_state_range");
        size_t max_func_evals = problem_db.get_sizet(
          "method.max_function_evaluations");

        if(num_interp_points <= 0) {
          Cerr << "ADAPTIVE SOGA: Number of interpolation points "
               << "not specified. These need to be under continuous "
               << "state variables for adaptive optimization.\n";
          abort_handler(METHOD_ERROR);
        }

        IntVector initial_state(num_interp_points);
        IntVector lower_bounds(num_interp_points);
        IntVector upper_bounds(num_interp_points);
        StringArray updated_labels(num_interp_points);
        BitArray updated_catgry(num_interp_points);

        //! These descriptors are recognized by SOFICS
        
        updated_labels[0] = "TARGET";
        updated_catgry[0] = true;
        initial_state[0]  = -1;
        lower_bounds[0]   = -1;
        upper_bounds[0]   = max_func_evals;

        for(size_t i=1; i<num_interp_points; ++i) {
          updated_labels[i] = "NEIGHBOR_" + std::to_string(i);
          updated_catgry[i] = true;

          initial_state[i] = -1;
          lower_bounds[i]  = -1;
          upper_bounds[i]  = max_func_evals;
        }
        
        //! Update the discrete STATE range variables
        problem_db.set("variables.discrete_state_range.labels", 
          updated_labels);
        problem_db.set("variables.discrete_state_range.categorical", 
          updated_catgry);
        problem_db.set("variables.discrete_state_range.initial_state",
          initial_state);
        problem_db.set("variables.discrete_state_range.lower_bounds",
          lower_bounds);
        problem_db.set("variables.discrete_state_range.upper_bounds",
          upper_bounds);

        //---------------------------------------------------------------------
        // Populate ERROR model (Dakota classes initialized later)
        //---------------------------------------------------------------------

        //! Create an error model
        data_model.data_rep()->idModel = "ERROR_MODEL";
        data_model.data_rep()->modelType = "simulation";
        data_model.data_rep()->interfacePointer = "ERROR_FORK_INTERF";
        data_model.data_rep()->responsesPointer = "MSE_RESPONSE";

        //! Create an interface for computing MSE
        data_interf.data_rep()->idInterface = "ERROR_FORK_INTERF";
        data_interf.data_rep()->interfaceType = Dakota::FORK_INTERFACE;
        data_interf.data_rep()->fileTagFlag = true;
        data_interf.data_rep()->fileSaveFlag = true;
        data_interf.data_rep()->useWorkdir = true;
        data_interf.data_rep()->dirTag = true;
        data_interf.data_rep()->dirSave = true;
        data_interf.data_rep()->failAction = "abort";
        data_interf.data_rep()->asynchLocalEvalConcurrency = 1;
        data_interf.data_rep()->asynchLocalEvalScheduling = Dakota::STATIC_SCHEDULING;

        //! Get info from True model
        const String &tr_work_dir = 
          problem_db.get_string("interface.workDir");
        const String &tr_param_file = 
          problem_db.get_string("interface.application.parameters_file");
        const String &tr_result_file = 
          problem_db.get_string("interface.application.results_file");
        const StringArray &tr_drivers = 
          problem_db.get_sa("interface.application.analysis_drivers");

        //! Set rest of the interface values.
        data_interf.data_rep()->workDir = "error_simulations/" + tr_work_dir;
        data_interf.data_rep()->parametersFile = tr_param_file;
        data_interf.data_rep()->resultsFile = tr_result_file;
        data_interf.data_rep()->analysisDrivers = tr_drivers;

        //! Create a set of responses for converying MSE to Optimizer
        data_respns.data_rep()->idResponses = "MSE_RESPONSE";
        data_respns.data_rep()->numResponseFunctions = 1;
        data_respns.data_rep()->responseLabels.push_back("MSE");

        //---------------------------------------------------------------------
        // Insert ERROR model into ProblemDescDB
        //---------------------------------------------------------------------
        problem_db.insert_node(data_model);
        problem_db.insert_node(data_interf);
        problem_db.insert_node(data_respns);

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
  shared_ptr<Model> top_model = topLevelIterator->iterated_model();

  //---------------------------------------------------------------------------
  // Setup top method iterator
  //---------------------------------------------------------------------------

  if(adaptive_optimization) {
  
    //! Make sure that STATE variables are set to inactive
    top_model->inactive_view(MIXED_STATE);

    //! Reference for reseting communicator
    ParLevLIter w_pl_iter = parallel_lib.w_parallel_level_iterator();
  
    //! For restoration
    size_t model_index = problem_db.get_db_model_node();
    problem_db.set_db_model_nodes("ERROR_MODEL");

    //! Create a new model
    std::shared_ptr<Model> err_model = problem_db.get_model();
    err_model->inactive_view(MIXED_STATE);

/*
    //! Reset to our own interface. This is done to distinguish
    //! working directories of true models and error models.
    err_model->derived_interface(
      std::make_shared<ForkApplicInterfaceWrapper>(problem_db));
*/

    //! Reset model pointer in problem db
    problem_db.set_db_model_nodes(model_index);

/*
    //! debug
    StringMultiArray string_var(boost::extents[1]);
    string_var[0] = "ERROR";
    const size_t idsv_len = string_var.num_elements();
    StringMultiArrayConstView idsv_view = string_var[
      boost::indices[idx_range(0,idsv_len)]];
    ModelUtils::inactive_discrete_string_variables(*err_model, idsv_view);

    StringMultiArrayConstView res = 
      ModelUtils::inactive_discrete_string_variables(*err_model);

    Cout << res[0] << "\n";

    abort_handler(OTHER_ERROR);
*/

    //! Assign our optimizer.
    topLevelIterator.reset(
       new AdaptiveJegaOptimizer(problem_db, top_model, err_model));

    //! Reset communicator
    topLevelIterator->init_communicators(w_pl_iter);
  
  }

}

//-----------------------------------------------------------------------------

