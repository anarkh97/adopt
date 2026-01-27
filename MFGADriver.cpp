/************************************************************************
 * MODIFICATIONS IN THIS DERIVED FILE
 *
 * Derived from:
 *   - Dakota's JEGAOptimizer::Driver
 *   - Upstream commit: 0daaaa2237bd79445349c47d157a4a2b73db7452
 *
 * Modifications by Aditya Narkhede, 2025:
 *
 * Licensing:
 *   This file continues to be licensed under the GNU Lesser General Public
 *   License v2.1 or (at your option) any later version. See the bundled
 *   license file `LICENSE.LGPL-2.1` or
 *   <https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html>.
 ************************************************************************/

#include <MFGADriver.h>

using namespace std;
using namespace JEGA::Logging;
using namespace JEGA::FrontEnd;
using namespace eddy::utilities;
using namespace JEGA::Utilities;
using namespace JEGA::Algorithms;

//-----------------------------------------------------------------------------

namespace MultiFidelityOptimizer
{

namespace detail
{

MFGADriver::MFGADriver(const ProblemConfig &config)
    : JEGA::FrontEnd::Driver(config)
{
  EDDY_FUNC_DEBUGSCOPE
}

} // namespace detail

} // namespace MultiFidelityOptimizer

//-----------------------------------------------------------------------------

