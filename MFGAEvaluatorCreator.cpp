/************************************************************************
 * MODIFICATIONS IN THIS DERIVED FILE
 *
 * Derived from:
 *   - Dakota's JEGAOptimizer::EvaluatorCreator
 *   - Upstream commit: 0daaaa2237bd79445349c47d157a4a2b73db7452
 *
 * Modifications by Aditya Narkhede, 2025:
 *   - Included a seperate Dakota::Model for error calculations.
 *   - Included a AdaptiveDecisionMaker
 *
 * Licensing:
 *   This file continues to be licensed under the GNU Lesser General Public
 *   License v2.1 or (at your option) any later version. See the bundled
 *   license file `LICENSE.LGPL-2.1` or
 *   <https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html>.
 ************************************************************************/

#include <MFGAEvaluatorCreator.h>
#include <MFGAEvaluator.h>
#include <utilities/include/EDDY_DebugScope.hpp>

using namespace std;
using namespace Dakota;
using namespace JEGA::FrontEnd;
using namespace JEGA::Algorithms;

namespace MultiFidelityOptimizer
{

namespace detail
{

//-----------------------------------------------------------------------------

MFGAEvaluatorCreator::MFGAEvaluatorCreator(Model &tmodel_, Model &emodel_,
                                           AdaptiveDecisionMaker &dmaker_)
    : sim_model(tmodel_), error_model(emodel_),
      decision_maker(dmaker_)
{
  EDDY_FUNC_DEBUGSCOPE
}

//-----------------------------------------------------------------------------

//! Return our Evaluator model.
GeneticAlgorithmEvaluator* 
MFGAEvaluatorCreator::CreateEvaluator(GeneticAlgorithm &alg)
{
  EDDY_FUNC_DEBUGSCOPE
  return new MFGAEvaluator(alg, sim_model, error_model, decision_maker);
}

} // detail

} // MultiFidelityOptimizer
