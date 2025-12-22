#include <cassert>
#include <string>

#include <DakotaEnvironmentWrapper.h>

#ifndef DAKOTA_HAVE_MPI
#define MPI_COMM_WORLD 0
#endif

//-----------------------------------------------------------------------------

int main(int argc, char *argv[])
{

  //! whether dakota running in parallel (e.g. mpiexec adopt -i adopt.in)
  bool parallel = Dakota::MPIManager::detect_parallel_launch(argc, argv);

  //! Define MPI_DEBUG in dakota_global_defs.cpp to cause a hold here
  Dakota::mpi_debug_hold();

#ifdef DAKOTA_HAVE_MPI
  //! AN: Stop execution in case someone tries mpiexec adopt.
  //! AN: Currently, we are assuming Dakota's operations are cheap
  //! compared to the function evaluations (i.e., simulations).
  //! Consequently, Dakota can be on one processor.
  if (parallel)
  {
    //MPI_Init(&argc, &argv); // initialize MPI
    Cerr << "***Error: ADOPT is not configured to run with "
         << "Dakota's coarse/fine grained parallelism." << std::endl;
    Dakota::abort_handler(Dakota::OTHER_ERROR);
  }
#endif // DAKOTA_HAVE_MPI

  //! Parse command line through Dakota
  //! Note that we are assuming that we are on the root rank.
  Dakota::ProgramOptions opts(argc, argv, 0);

  //! Create the Library environment which sets the Data containers
  DakotaEnvironmentWrapper env(opts, false);

  //! Perform whatever study was requested.
  env.execute();

#ifdef DAKOTA_HAVE_MPI
  //if (parallel)
  //MPI_Finalize(); // finalize MPI
#endif // DAKOTA_HAVE_MPI

  return EXIT_SUCCESS;
}

//-----------------------------------------------------------------------------
//
