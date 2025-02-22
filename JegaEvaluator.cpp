#include<JegaEvaluator.h>

// JEGAConfig.hpp should be the first include in all JEGA files.
#include <../Utilities/include/JEGAConfig.hpp>
#include <../Utilities/include/Logging.hpp>

// JEGA utility includes.
#include <../Utilities/include/DesignGroup.hpp>
#include <../Utilities/include/ConstraintInfo.hpp>

using namespace std;
using namespace Dakota;
using namespace JEGA::Logging;
using namespace eddy::utilities;
using namespace JEGA::Utilities;
using namespace JEGA::Algorithms;

//-----------------------------------------------------------------------------

JegaEvaluator::JegaEvaluator(GeneticAlgorithm &algorithm, Model &tmodel_, 
                             Model &emodel_)
             : GeneticAlgorithmEvaluator(algorithm), true_model(tmodel_), 
               error_model(emodel_)
{
  EDDY_FUNC_DEBUGSCOPE
}

//-----------------------------------------------------------------------------

JegaEvaluator::JegaEvaluator(const JegaEvaluator &copy)
             : GeneticAlgorithmEvaluator(copy), true_model(copy.true_model),
               error_model(copy.error_model)
{
  EDDY_FUNC_DEBUGSCOPE
}

//-----------------------------------------------------------------------------

JegaEvaluator::JegaEvaluator(const JegaEvaluator &copy, GeneticAlgorithm &algo,
                             Model &tmodel_, Model &emodel_)
             : GeneticAlgorithmEvaluator(copy, algo), true_model(tmodel_),
               error_model(emodel_)
{
  EDDY_FUNC_DEBUGSCOPE
}

//-----------------------------------------------------------------------------

bool JegaEvaluator::Evaluate(Design &des)
{
  EDDY_FUNC_DEBUGSCOPE
  JEGALOG_II_F(GetLogger(), this, text_entry(lfatal(), GetName() + 
               ": You cannot use Evaluate(Design&) with this "
               "evaluator...ever."))
  return false;
}

//-----------------------------------------------------------------------------

string JegaEvaluator::GetName() const
{
  EDDY_FUNC_DEBUGSCOPE
  return JegaEvaluator::Name();
}

//-----------------------------------------------------------------------------

string JegaEvaluator::GetDescription() const
{
  EDDY_FUNC_DEBUGSCOPE
  return JegaEvaluator::Description();
}

//-----------------------------------------------------------------------------

GeneticAlgorithmOperator*
JegaEvaluator::Clone(GeneticAlgorithm &algorithm) const
{
  EDDY_FUNC_DEBUGSCOPE
  return new JegaEvaluator(*this, algorithm, true_model, error_model);
}

//-----------------------------------------------------------------------------

size_t JegaEvaluator::GetNumberNonLinearConstraints() const
{
  EDDY_FUNC_DEBUGSCOPE
  return ModelUtils::num_nonlinear_eq_constraints(true_model) +
         ModelUtils::num_nonlinear_ineq_constraints(true_model);
}

//-----------------------------------------------------------------------------

size_t JegaEvaluator::GetNumberLinearConstraints() const
{
  EDDY_FUNC_DEBUGSCOPE
  return ModelUtils::num_linear_eq_constraints(true_model) +
         ModelUtils::num_linear_ineq_constraints(true_model);
}

//-----------------------------------------------------------------------------

void
JegaEvaluator::SeparateVariables(const Design &from, RealVector &into_cont) const
{

  EDDY_FUNC_DEBUGSCOPE

  size_t num_cv  = ModelUtils::cv(true_model);
  size_t num_div = ModelUtils::div(true_model);
  size_t num_drv = ModelUtils::drv(true_model);
  size_t num_dsv = ModelUtils::dsv(true_model);

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
//! Refactored from JEGAOptimizer::Evaluator::Evaluate
bool
JegaEvaluator::Evaluate(DesignGroup &group)
{
  EDDY_FUNC_DEBUGSCOPE

  JEGALOG_II(GetLogger(), ldebug(), this,
    text_entry(ldebug(), GetName() + ": Performing group evaluation."))

  // check for trivial abort conditions
  if(group.IsEmpty()) return true;

  // first, let's see if we can avoid any evaluations.
  ResolveClones(group);

  // we'll prepare containers for repeated use without re-construction
  RealVector       continuous_variables;
/*
  IntVector        disc_int_vars;
  RealVector       disc_real_vars;
  StringMultiArray disc_string_vars;
*/

  // prepare to iterate over the group
  DesignDVSortSet::const_iterator it(group.BeginDV());
  const DesignDVSortSet::const_iterator e(group.EndDV());

  // these quantities will be used below
  const DesignTarget& target = GetDesignTarget();

  // Find out the counts on the different types of constraints
  const size_t num_nonlin_cn = GetNumberNonLinearConstraints();
  // const size_t num_lin_cn = GetNumberLinearConstraints();

  // Get the information about the constraints.
  const ConstraintInfoVector& cninfos = target.GetConstraintInfos();

  // Prepare to store the number of requests in order to avoid overshooting
  // the limit.
  const eddy::utilities::uint64_t prior_reqs = GetNumberEvaluations();
  eddy::utilities::uint64_t num_eval_reqs = 0;

  // store an iterator to the first linear constraint so that
  // we can evaluate it using the cinfo objects.

  // prepare to iterate
  ConstraintInfoVector::const_iterator flincn(cninfos.begin() + num_nonlin_cn);
  ConstraintInfoVector::const_iterator cit;

  // prepare to return the success of this.  Success occurs only if all
  // designs wind up evaluated and non-illconditioned.
  bool ret = true;

  for(; it!=e; ++it) {
    // If this Design is evaluated, let's skip it.
    if((*it)->IsEvaluated()) continue;

    // If we've reached our maximum allowable number of
    // evaluations, tag remaining as evaluated and illconditioned.
    // By doing so, they will be flushed from the algorithm.
    if((prior_reqs + num_eval_reqs) >= GetMaxEvaluations()) {
      (*it)->SetEvaluated(true);
      (*it)->SetIllconditioned(true);
      ret = false;
      continue;
    }

    // extract the real and continuous variables
    // from the current Design
    SeparateVariables(**it, continuous_variables /*, disc_int_vars, disc_real_vars,
        disc_string_vars*/);

    // send this guy out for evaluation using the "true_model"

    // first, set the current values of the variables in the model
    ModelUtils::continuous_variables(true_model, continuous_variables);
    //ModelUtils::discrete_int_variables(true_model, disc_int_vars);
    //ModelUtils::discrete_real_variables(true_model, disc_real_vars);
    //// Strings set by calling single value setter for each
    ////for (size_t i=0; i<disc_string_vars.num_elements(); ++i)
    ////  ModelUtils::discrete_string_variable(true_model, disc_string_vars[i],i);
    //// Could use discrete_string_varables to avoid overhead of repeated 
    //// function calls, but it takes a StringMultiArrayConstView, which
    //// must be created from disc_string_vars. Maybe there's a simpler way,
    //// but...
    //// const size_t &dsv_len = disc_string_vars.num_elements();
    //// StringMultiArrayConstView dsv_view = disc_string_vars[ 
    ////   boost::indices[idx_range(0,dsv_len)]];
    //// ModelUtils::discrete_string_variables(this->_model, dsv_view);
    
    // now request the evaluation in synchronous or asyncronous mode.
    if(true_model.asynch_flag()) {
      // The following method call will use the default
      // Active set vector which is to just compute
      // function values, no gradients or hessians.
      true_model.evaluate_nowait();
    }
    else {
      // The following method call will use the default
      // Active set vector which is to just compute
      // function values, no gradients or hessians.
      true_model.evaluate();

      // increment the number of performed evaluations by 1
      IncrementNumberEvaluations();

      // Record the responses back into the Design
      const RealVector& ftn_vals =
          true_model.current_response().function_values();

      RecordResponses(ftn_vals, **it);

      // Label this guy as now being evaluated.
      (*it)->SetEvaluated(true);

      // now check the feasibility of this design
      target.CheckFeasibility(**it);
    }

    // no matter what, we want to increment the number of evaluations
    // or evaluation requests.
    ++num_eval_reqs;

    // The responses do not (or will not) include the linear
    // constraint values. We have to compute them ourselves.
    // we can do it using the "EvaluateConstraint" method which
    // will only succeed for linear constraints for which
    // coefficients have been supplied.
    for(cit=flincn; cit!=cninfos.end(); ++cit) {
      (*cit)->EvaluateConstraint(**it);
      (*cit)->RecordViolation(**it);
    }
  }

  // If we did our evaluations asynchronously, we did not yet record
  // the results (because they were not available).  We need to do so
  // now.  The call to _model.synchronize causes the program to block
  // until all the results are available.  We can then record the
  // responses in the same fashion as above.  Note that the linear
  // constraints have already been computed!!
  if(true_model.asynch_flag()) {
    // Wait for the responses.
    const IntResponseMap& response_map = true_model.synchronize();
    size_t num_resp = response_map.size();

    // increment the number of evaluations by the number of responses
    IncrementNumberEvaluations(num_resp);

    EDDY_ASSERT(num_resp == num_eval_reqs);

    // prepare to access the elements of the response_map by iterator.
    IntRespMCIter r_cit = response_map.begin();

    // Record the set of responses in the DesignGroup
    for(it=group.BeginDV(); it!=e; ++it) {
      // we didn't send already-evaluated Designs out for evaluation
      // so skip them here as well.
      if((*it)->IsEvaluated()) continue;

      // Put the responses into the Design properly.
      RecordResponses(r_cit->second.function_values(), **it);

      // Label this guy as now being evaluated.
      (*it)->SetEvaluated(true);

      // now check the feasibility of this design
      target.CheckFeasibility(**it);

      // only increment for newly evaluated points contained in response_map
      ++r_cit;
    }
  }

  // Now get the error estimate at each true evaluation using the error_model
  //if(ret) {

  //}

  return ret;
}


//-----------------------------------------------------------------------------


