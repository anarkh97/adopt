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

