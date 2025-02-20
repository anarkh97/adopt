#include<JegaEvaluatorCreator.h>
#include<JegaEvaluator.h>
#include<utilities/include/EDDY_DebugScope.hpp>

using namespace std;
using namespace Dakota;
using namespace JEGA::FrontEnd;
using namespace JEGA::Algorithms;

//-----------------------------------------------------------------------------

JegaEvaluatorCreator::JegaEvaluatorCreator(Model &tmodel_, Model &emodel_)
                    : true_model(tmodel_), error_model(emodel_)
{
  EDDY_FUNC_DEBUGSCOPE 
}

//-----------------------------------------------------------------------------

//! Return our Evaluator model.
GeneticAlgorithmEvaluator*
JegaEvaluatorCreator::CreateEvaluator(GeneticAlgorithm &alg)
{

  EDDY_FUNC_DEBUGSCOPE
  return new JegaEvaluator(alg, true_model, error_model);

}
