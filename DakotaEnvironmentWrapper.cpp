#include<DakotaEnvironmentWrapper.h>
#include<AdaptiveJegaOptimizer.h>
#include<ForkApplicInterfaceWrapper.h>

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
      const StringArray &metadata_labels = 
        problem_db.get_sa("responses.metadata_labels");

      bool found_switch_label = false;
      bool found_mse_label = false;

      //! Find the SWITCH discrete string set state variable
      for(const auto &label : state_string_labels) {
        if(label == "SWITCH") {
          found_switch_label = true;
          break;
        }
      }

      //! Find the MSE metadata response
      for(const auto &label : metadata_labels) {
        if(label == "MSE") {
          found_mse_label = true;
          break;
        }
      }

      //! Switch to adaptive optimization if both labels are found
      adaptive_optimization = (found_switch_label and found_mse_label);

      if(adaptive_optimization) {

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

  shared_ptr<Model> debug_err_model;

  if(adaptive_optimization) {
  
    //! Make sure that STATE variables are set to inactive
    top_model->inactive_view(MIXED_STATE);

    //! Reference for reseting communicator
    ParLevLIter w_pl_iter = parallel_lib.w_parallel_level_iterator();
  
    //! Create a new Dakota Model for error evaluations.
    shared_ptr<Model> err_model = make_shared<SimulationModel>(problem_db);
    //err_model->inactive_view(MIXED_STATE);
    ModelUtils::active_variables(*err_model, top_model->current_variables());
  
    //! Reset to our own interface. This is done to distinguish
    //! working directories of true models and error models.
    Interface &err_model_interface = err_model->derived_interface();   
    err_model_interface.assign_rep(
      std::make_shared<ForkApplicInterfaceWrapper>(problem_db));

    debug_err_model = err_model;

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
    topLevelIterator.assign_rep(
      make_shared<AdaptiveJegaOptimizer>(problem_db, top_model, err_model));

    //! Reset communicator
    topLevelIterator.init_communicators(w_pl_iter);
  
  }

/*
  //! debug
  RealVector cont_vars(3);
  IntVector int_vars(1);
  StringMultiArray string_var(boost::extents[1]);

  cont_vars[0] = cont_vars[1] = cont_vars[2] = 1.0;
  int_vars[0] = -1;
  string_var[0] = "TRUE";

  ModelUtils::continuous_variables(*top_model, cont_vars);
  ModelUtils::inactive_discrete_int_variables(*top_model, int_vars);
  StringMultiArrayConstView idsv_view = string_var[
    boost::indices[idx_range(0,1)]];
  ModelUtils::inactive_discrete_string_variables(*top_model, idsv_view);

  top_model->evaluate();

  string_var[0] = "ERROR";
  int_vars[0] = 1;
  ModelUtils::continuous_variables(*debug_err_model, cont_vars);
  ModelUtils::inactive_discrete_int_variables(*debug_err_model, int_vars);
  StringMultiArrayConstView idsv_view2 = string_var[
    boost::indices[idx_range(0,1)]];
  ModelUtils::inactive_discrete_string_variables(*debug_err_model, idsv_view2);

  debug_err_model->evaluate();

  abort_handler(OTHER_ERROR);
*/
}

//-----------------------------------------------------------------------------

