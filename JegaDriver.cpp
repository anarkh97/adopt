#include<JegaDriver.h>

using namespace std;
using namespace JEGA::Logging;
using namespace JEGA::FrontEnd;
using namespace eddy::utilities;
using namespace JEGA::Utilities;
using namespace JEGA::Algorithms;

//-----------------------------------------------------------------------------

JegaDriver::JegaDriver(const ProblemConfig &config)
          : JEGA::FrontEnd::Driver(config)
{
  EDDY_FUNC_DEBUGSCOPE
}

//-----------------------------------------------------------------------------

//! Refactored from JEGA::FrontEnd::Driver::PerformIterations
DesignOFSortSet
JegaDriver::PerformIterations(GeneticAlgorithm *ga)
{

  EDDY_FUNC_DEBUGSCOPE

  try{

    // prepare to time this run.
    JEGA_LOGGING_IF_ON(clock_t start = clock();)

    // initialize the populations
    ga->AlgorithmInitialize();

    while(true) {

      // This will call JegaEvaluator class through JEGA.
      if(!PerformNextIteration(ga)) break;

      // Place Adaptive decision model here.

    }

    // Clean-ups on JEGA side.
    ga->AlgorithmFinalize();

    JEGA_LOGGING_IF_ON(double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;)

    JEGALOG_II(ga->GetLogger(), lquiet(), this, ostream_entry(lquiet(), 
               "JEGA Front End: " + ga->GetName() + " execution took ") 
               << elapsed << " seconds.")

    JEGALOG_II_G(lquiet(), this, ostream_entry(lquiet(), "JEGA Front End: "
                 "Execution took ") << elapsed << " seconds.")

    return DeepDuplicate(ga->GetCurrentSolution());

  } catch(const exception &e) {
    JEGALOG_II_G_F(this, text_entry(lfatal(), "JEGA Front End Error: "
                   "Exception caught at application level reading: \"") 
                   << e.what() << "\".")
  }

  return DesignOFSortSet();

}

//-----------------------------------------------------------------------------

