#include<ForkApplicInterfaceWrapper.h>

#include<DakotaResponse.hpp>
#include<ParamResponsePair.hpp>
#include<ForkApplicInterface.hpp>
#include<ProblemDescDB.hpp>
#include<ParallelLibrary.hpp>
#include<WorkdirHelper.hpp>
#include<sys/wait.h> // for wait and waitpid
#include<unistd.h>   // for fork, execvp, setgpid
#include<algorithm>
#include<thread>


//-----------------------------------------------------------------------------

ForkApplicInterfaceWrapper::ForkApplicInterfaceWrapper(
                              const Dakota::ProblemDescDB &problem_db)
                          : Dakota::ForkApplicInterface(problem_db)
{
  // error simulations are stored in a separate sub-directory
  workDirName = "error_simulations/" + workDirName;
}

//-----------------------------------------------------------------------------

