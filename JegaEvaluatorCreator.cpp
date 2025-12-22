#include <JegaEvaluatorCreator.h>
#include <JegaEvaluator.h>
#include <utilities/include/EDDY_DebugScope.hpp>

using namespace std;
using namespace Dakota;
using namespace JEGA::FrontEnd;
using namespace JEGA::Algorithms;

//-----------------------------------------------------------------------------

JegaEvaluatorCreator::JegaEvaluatorCreator(Model &tmodel_, Model &emodel_,
                                           AdaptiveDecisionMaker &dmaker_)
    : sim_model(tmodel_), error_model(emodel_),
      decision_maker(dmaker_){EDDY_FUNC_DEBUGSCOPE}

        //-----------------------------------------------------------------------------

        //! Return our Evaluator model.
        GeneticAlgorithmEvaluator
        * JegaEvaluatorCreator::CreateEvaluator(GeneticAlgorithm & alg)
{
  EDDY_FUNC_DEBUGSCOPE
  return new JegaEvaluator(alg, sim_model, error_model, decision_maker);
}
