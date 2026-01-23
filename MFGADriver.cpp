#include <MFGADriver.h>

using namespace std;
using namespace JEGA::Logging;
using namespace JEGA::FrontEnd;
using namespace eddy::utilities;
using namespace JEGA::Utilities;
using namespace JEGA::Algorithms;

//-----------------------------------------------------------------------------

namespace MultiFidelityOptimizer {

namespace detail {

MFGADriver::MFGADriver(const ProblemConfig &config)
    : JEGA::FrontEnd::Driver(config)
{
  EDDY_FUNC_DEBUGSCOPE
}

}

}

//-----------------------------------------------------------------------------

