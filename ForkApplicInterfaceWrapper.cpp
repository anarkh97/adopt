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
  // concatenate _error to work directory
  workDirName = workDirName + "_error";
}

//-----------------------------------------------------------------------------

