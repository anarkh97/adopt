#include<JegaEvaluator.h>

using namespace std;
using namespace Dakota;
using namespace JEGA::Logging;
using namespace JEGA::FrontEnd;
using namespace eddy::utilities;
using namespace JEGA::Utilities;
using namespace JEGA::Algorithms;

//-----------------------------------------------------------------------------

JegaEvaluator::JegaEvaluator(GeneticAlgorithm &algorithm, Model &model_)
             : GeneticAlgorithmEvaluator(algorithm), model(model_)
{
  EDDY_FUNC_DEBUGSCOPE
}

//-----------------------------------------------------------------------------

JegaEvaluator::JegaEvaluator(const JegaEvaluator &copy)
             : GeneticAlgorithmEvaluator(copy), model(copy.model)
{
  EDDY_FUNC_DEBUGSCOPE
}

//-----------------------------------------------------------------------------

JegaEvaluator::JegaEvaluator(const JegaEvaluator &copy, GeneticAlgorithm &algo,
                             Model &model_)
             : GeneticAlgorithmEvaluator(copy, algorithm), model(model_)
{
  EDDY_FUNC_DEBUGSCOPE
}

//-----------------------------------------------------------------------------

void
JegaEvaluator::SeparateVariables(const Design &from, RealVector &into_cont) const
{

}

//-----------------------------------------------------------------------------

void
JegaEvaluator::RecordResponses(const RealVector &from, Design &into) const
{

}

//-----------------------------------------------------------------------------

//! Performs design evaluations using the model interface.
void
JegaEvaluator::Evaluate(DesignGroup &group)
{

}

//-----------------------------------------------------------------------------


