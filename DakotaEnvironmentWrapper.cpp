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

  if(top_method_name == "soga") {

    // calling it here explicitly to make things clear and readable.
    ProblemDescDB &problem_db = problem_description_db();

    // get all the available models
/*
    ModelList &models = problem_db.model_list();
    for(auto &m : models) {

      const String &m_id     = m->model_id();
      const String &m_type   = m->model_type();
      const String &m_source = m->method_source();

      Cout << "Found Model with\n"
           << "ID: " << m_id << std::endl
	   << "Type: " << m_type << std::endl
	   << "Source: " << m_source << std::endl;
    }

    abort_handler(-1);
*/

    shared_ptr<Model> top_model = 
      topLevelIterator.iterated_model();

    // ProblemDescDB::set does not allow tinkering of analysis
    // drivers. So, we simply create a new Model with the same
    // interface. We handle different kinds of evaluations 
    // in the driver through STATE variables.

    shared_ptr<Model> err_model =
      ModelUtils::get_model(problem_db);
/*
    // print drivers to verify --- they should be the same.
    Interface &true_interface = top_model->derived_interface();
    Interface &esti_interface = err_model->derived_interface();

    for(auto &driver : true_interface.analysis_drivers()) 
      Cout << "True driver: " << driver << std::endl;
    
    for(auto &driver : esti_interface.analysis_drivers())
      Cout << "Esti driver: " << driver << std::endl;

    abort_handler(-1);
*/

    // add continuous state set variables w/ permissible values
    // TRUE, INTERP, ERROR


    // assign our optimizer.
    topLevelIterator.assign_rep(
      make_shared<AdaptiveJegaOptimizer>(problem_db, top_model, err_model));

  }

}

//-----------------------------------------------------------------------------

