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

JegaEvaluator::JegaEvaluator(GeneticAlgorithm &alg)
{

}

//-----------------------------------------------------------------------------

void
JegaEvaluator::SeparateVariables(const Design &from, RealVector &into_cont) const
{

  EDDY_FUNC_DEBUGSCOPE

  size_t num_cv  = model.cv();
  size_t num_div = model.div();
  size_t num_drv = model.drv();
  size_t num_dsv = model.dsv();

  // Currently we only support continuous real variables.
  if(num_div > 0) {
    JEGALOG_II_G_F(this, text_entry(lfatal(), "Adaptive JEGA Error: "
                   "Currently, Discrete Integer design variables are not "
                   "supported.\n"))
  }
  if(num_drv > 0) {
    JEGALOG_II_G_F(this, text_entry(lfatal(), "Adaptive JEGA Error: "
                   "Currently, Discrete Real design variables are not "
                   "supported.\n"))
  }
  if(num_dsv > 0) {
    JEGALOG_II_G_F(this, text_entry(lfatal(), "Adaptive JEGA Error: "
                   "Currently, Discrete String design variables are not "
                   "supported.\n"))
  }

  if(into_cont.length() != num_cv) into_cont.size(num_cv);
 
  // Get the target from JEGA 
  const DesignTarget &target = from.GetDesignTarget();
  const DesignVariableInfoVector &dv_info = target.GetDesignVariableInfos();

  for(int i=0; i<num_cv; ++i) {
    EDDY_ASSERT(dv_info[i]->IsContinuum());
    into_cont[i] = dv_info[i]->WhichValue(from);
  }

}

//-----------------------------------------------------------------------------

void
JegaEvaluator::RecordResponses(const RealVector &from, Design &into) const
{

  EDDY_FUNC_DEBUGSCOPE

  // get the target for information.
  const DesignTarget& target = GetDesignTarget();

  // get the information about the constraints.
  const ConstraintInfoVector& cnis = target.GetConstraintInfos();

  // prepare to store the location in the responses vector.
  RealVector::ordinalType loc = 0;

  // find out how many objective and constraint functions there are.
  const size_t nof = target.GetNOF();
  const size_t ncn = target.GetNCN();

  // record the objective functions first.
  for(size_t i=0; i<nof; ++i, ++loc)
    into.SetObjective(i, from[loc]);

  // now record the nonlinear constraints.  To do this,
  // we will need to know how many there are.  They will be the
  // first of the constraints in the design.
  const size_t num_nonlin_cn = GetNumberNonLinearConstraints();
  for(size_t cn=0; cn<num_nonlin_cn && cn<ncn; ++cn, ++loc) {
    into.SetConstraint(cn, from[loc]);
    cnis[cn]->RecordViolation(into);
  }
  
}

//-----------------------------------------------------------------------------

//! Performs design evaluations using the model interface.
void
JegaEvaluator::Evaluate(DesignGroup &group)
{
// TODO: This is where Dakota Model is called which launches simulations.
}

//-----------------------------------------------------------------------------


