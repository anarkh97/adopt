#ifndef _DAKOTA_ENVIRONMENT_WRAPPER_H_
#define _DAKOTA_ENVIRONMENT_WRAPPER_H_

// Dakota includes
#include<ParallelLibrary.hpp>
#include<ProblemDescDB.hpp>
#include<LibraryEnvironment.hpp>

/************************************************
 * A lightweight wrapper for Dakota's 
 * LibraryEnvironment class. In the future, this 
 * class can be expanded to implement a "front-end"
 * using our custom "IoData", while integrating a 
 * "back-end" that invokes M2C and Aero-S core 
 * functions (main function converted to "call-back
 * functions"). By implementing a dedicated back-end,
 * FSI simulations can efficiently run on MPI 
 * sub-groups managed by Dakota's evaluation 
 * server, eliminating complicated node balancing 
 * typically performed in the analysis driver.
 ***********************************************/

class DakotaEnvironmentWrapper : public Dakota::LibraryEnvironment {

public:

  DakotaEnvironmentWrapper(Dakota::ProgramOptions options = Dakota::ProgramOptions(),
                           bool check_bcast_construct = true,
                           Dakota::DbCallbackFunctionPtr callback = nullptr,
                           void* callback_data = nullptr);

  DakotaEnvironmentWrapper(Dakota::MPI_Comm dakota_comm,
                           Dakota::ProgramOptions options = Dakota::ProgramOptions(),
			   bool check_bcast_construct = true,
			   Dakota::DbCallbackFunctionPtr callback = nullptr,
                           void* callback_data = nullptr);

  //! Object deletion left to LibraryEnvironment for now
  ~DakotaEnvironmentWrapper() { }

private:

  void SetupAdaptiveOptimizer();

};

#endif
