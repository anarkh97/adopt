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

JegaEvaluator::JegaEvaluator(GeneticAlgorithm &algorithm, 
                             Model &tmodel_, 
                             Model &emodel_,
                             AdaptiveDecisionMaker &decision_maker_)
             : GeneticAlgorithmEvaluator(algorithm), sim_model(tmodel_), 
               error_model(emodel_), decision_maker(decision_maker_)
{
  EDDY_FUNC_DEBUGSCOPE
}

//-----------------------------------------------------------------------------

JegaEvaluator::JegaEvaluator(const JegaEvaluator &copy)
             : GeneticAlgorithmEvaluator(copy), sim_model(copy.sim_model),
               error_model(copy.error_model), decision_maker(copy.decision_maker)
{
  EDDY_FUNC_DEBUGSCOPE
}

//-----------------------------------------------------------------------------

JegaEvaluator::JegaEvaluator(const JegaEvaluator &copy, 
                             GeneticAlgorithm &algo,
                             Model &tmodel_, 
                             Model &emodel_,
                             AdaptiveDecisionMaker &decision_maker_)
             : GeneticAlgorithmEvaluator(copy, algo), sim_model(tmodel_),
               error_model(emodel_), decision_maker(decision_maker_)
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
  return new JegaEvaluator(*this, algorithm, sim_model, error_model);
}

//-----------------------------------------------------------------------------

size_t JegaEvaluator::GetNumberNonLinearConstraints() const
{
  EDDY_FUNC_DEBUGSCOPE
  return ModelUtils::num_nonlinear_eq_constraints(sim_model) +
         ModelUtils::num_nonlinear_ineq_constraints(sim_model);
}

//-----------------------------------------------------------------------------

size_t JegaEvaluator::GetNumberLinearConstraints() const
{
  EDDY_FUNC_DEBUGSCOPE
  return ModelUtils::num_linear_eq_constraints(sim_model) +
         ModelUtils::num_linear_ineq_constraints(sim_model);
}

//-----------------------------------------------------------------------------

void
JegaEvaluator::SeparateVariables(const Design &from, RealVector &into_cont) 
const
{

  EDDY_FUNC_DEBUGSCOPE

  size_t num_cv  = ModelUtils::cv(sim_model);
  size_t num_div = ModelUtils::div(sim_model);
  size_t num_drv = ModelUtils::drv(sim_model);
  size_t num_dsv = ModelUtils::dsv(sim_model);

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
JegaEvaluator::SetStateVariables(const RealVector& cont_vars, 
                                 IntVector &into_disc_int,
                                 StringMultiArray &into_disc_string,
                                 const bool error_flag)
{
    EDDY_FUNC_DEBUGSCOPE

/*
  size_t num_dsv  = ModelUtils::dsv(sim_model);
  size_t num_idsv = ModelUtils::idsv(sim_model);
  size_t num_adsv = ModelUtils::adsv(sim_model);

  Cout << "[Active]Discrete String variables = " << num_dsv << "\n"
       << "[Inactive]Discrete String variables = " << num_idsv << "\n"
       << "[Total]Discrete String variables = " << num_adsv << "\n";
  abort_handler(OTHER_ERROR);
*/

  size_t num_cv    = ModelUtils::cv(sim_model);
  size_t num_idiv  = ModelUtils::idiv(sim_model);
  size_t num_idsv  = ModelUtils::idsv(sim_model);

  EDDY_ASSERT(cont_vars.length() == num_cv);

  // allocate memory for arrays handled by this function.
  if(into_disc_int.length() != num_idiv) 
    into_disc_int.size(num_idiv);

  if(into_disc_string.num_elements() != num_idsv) {
    StringMultiArray::extent_gen extents;
    into_disc_string.resize(extents[num_idsv]);
  }
  
  // Prepare to set inactive (state) discrete string set variable
  StringMultiArrayConstView disc_string_labels = 
    ModelUtils::inactive_discrete_string_variable_labels(sim_model);
  bool found_label = false;

  // Set inactive discrete string set variable -> Evaluation Type
  // Set inactive discrete integer variable -> Neighbor evaluation IDs.  
  if(!error_flag) {

    decision_maker.GetNearestNeighbors(cont_vars, into_disc_int, num_idiv);
    bool eval_decision = 
      decision_maker.GetEvaluationDecision(cont_vars);

    // Pass evaluation decision to the analysis driver 
    // through discrete string set
    for(size_t i=0; i<num_idsv; ++i) {

      const auto &label = disc_string_labels[i];
      if(label == "SWITCH") {
        into_disc_string[i] = (eval_decision) ? "TRUE" : "APPROX";
        found_label         = true;
	break;
      }

    }

  }
  else {
    // Function was called for error_model
    decision_maker.GetNearestNeighbors(cont_vars,
      into_disc_int, num_idiv, true /*force-find*/);

    // Pass Error flag to the analysis driver 
    // through discrete string set
    for(size_t i=0; i<num_idsv; ++i) {

      const auto &label = disc_string_labels[i];
      if(label == "SWITCH") {
        into_disc_string[i] = "ERROR";
        found_label         = true;
	break;
      }

    }

  }
  
  if(!found_label) {
    JEGALOG_II_G_F(this, text_entry(lfatal(),
                 "Adaptive JEGA Error: Optimizer requires a "
                 "discrete string set state variable with "
                 "label \"SWITCH\".\n"))
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

void
JegaEvaluator::RecordEvaluationInDecisionMaker(const int id, 
                                               const RealVector &cont_vars,
                                               const StringMultiArray &disc_strings)
{

  EDDY_FUNC_DEBUGSCOPE

  // Identify the decision for current evaluation ID 

  bool eval_decision = true;
  size_t num_idsv    = ModelUtils::idsv(sim_model);

  // Prepare to set inactive (state) discrete string set variable
  StringMultiArrayConstView disc_string_labels = 
    ModelUtils::inactive_discrete_string_variable_labels(sim_model);
  bool found_label = false;

  for(size_t i=0; i<num_idsv; ++i) {

    const auto &label = disc_string_labels[i];
    if(label == "SWITCH") {
      eval_decision = (disc_strings[i] == "TRUE");
      found_label         = true;
      break;
    }

  }

  if(!found_label) {
    JEGALOG_II_G_F(this, text_entry(lfatal(),
                 "Adaptive JEGA Error: Optimizer requires a "
                 "discrete string set state variable with "
                 "label \"SWITCH\".\n"))
  }

  // Update the decision maker w/ evaluation id and variables.
  decision_maker.RecordEvaluationDecision(id, cont_vars, eval_decision);
  
}

//-----------------------------------------------------------------------------

void
JegaEvaluator::RecordErrorInDecisionMaker(const int id, /* evaluation id from true model */
                                          const std::vector<RespMetadataT> &vals, 
                                          const StringArray &labels, 
                                          const RealVector &cont_vars)
{

  EDDY_FUNC_DEBUGSCOPE

  size_t num_cv = ModelUtils::cv(sim_model);

  EDDY_ASSERT(cont_vars.length() == num_cv);

  size_t loc = 0;
  bool found_label = false;
  for(const auto &lbl : labels) {
    if(lbl == "MSE" and !found_label) {
      decision_maker.RecordEvaluationError(id, cont_vars, vals[loc]);
      found_label = true;
    }
    ++loc;
  }

  if(!found_label) {
    JEGALOG_II_G_F(this, text_entry(lfatal(),
                 "Adaptive JEGA Error: Optimizer requires a "
                 "discrete string set state variable with "
                 "label \"SWITCH\".\n"))
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

  // arrays/vectors for adaptive decision model
  IntVector        state_int_vars;
  StringMultiArray state_string_vars; // "ERROR", "APPROX", "TRUE"

  // prepare to iterate over the group
  DesignDVSortSet::const_iterator it(group.BeginDV());
  const DesignDVSortSet::const_iterator e(group.EndDV());

  // these quantities will be used below
  const DesignTarget& target = GetDesignTarget();

  // Find out the counts on the different types of constraints
  const size_t num_nonlin_cn = GetNumberNonLinearConstraints();
  // NOT COMMENTED BY AN
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

    // send this guy out for evaluation using the "sim_model"

    // first, set the current values of the variables in the model
    ModelUtils::continuous_variables(sim_model, continuous_variables);
    //ModelUtils::discrete_int_variables(sim_model, disc_int_vars);
    //ModelUtils::discrete_real_variables(sim_model, disc_real_vars);
    //// Strings set by calling single value setter for each
    ////for (size_t i=0; i<disc_string_vars.num_elements(); ++i)
    ////  ModelUtils::discrete_string_variable(sim_model, disc_string_vars[i],i);
    //// Could use discrete_string_varables to avoid overhead of repeated 
    //// function calls, but it takes a StringMultiArrayConstView, which
    //// must be created from disc_string_vars. Maybe there's a simpler way,
    //// but...
    //// const size_t &dsv_len = disc_string_vars.num_elements();
    //// StringMultiArrayConstView dsv_view = disc_string_vars[ 
    ////   boost::indices[idx_range(0,dsv_len)]];
    //// ModelUtils::discrete_string_variables(this->_model, dsv_view);
    
    //=========================================================================
    // Pass decision to analysis driver.
    //=========================================================================
    SetStateVariables(continuous_variables, state_int_vars, state_string_vars);
    ModelUtils::inactive_discrete_int_variables(sim_model, state_int_vars);
    const size_t idsv_len = state_string_vars.num_elements();
    StringMultiArrayConstView idsv_view = state_string_vars[
      boost::indices[idx_range(0,idsv_len)]];
    ModelUtils::inactive_discrete_string_variables(sim_model, idsv_view);

    // now request the evaluation in synchronous or asyncronous mode.
    if(sim_model.asynch_flag()) {
      // The following method call will use the default
      // Active set vector which is to just compute
      // function values, no gradients or hessians.
      sim_model.evaluate_nowait();

      // Always call this after the derived interface 
      // has mapped parameters with responses, i.e., once
      // evaluation id has been updated.
      int eval_id = sim_model.evaluation_id();
      RecordEvaluationInDecisionMaker(eval_id, 
        continuous_variables, state_string_vars);
    }
    else {
      // The following method call will use the default
      // Active set vector which is to just compute
      // function values, no gradients or hessians.
      sim_model.evaluate();

      // Always call this after the derived interface 
      // has mapped parameters with responses, i.e., once
      // evaluation id has been updated.
      int eval_id = sim_model.evaluation_id();
      RecordEvaluationInDecisionMaker(eval_id, 
        continuous_variables, state_string_vars);

      // increment the number of performed evaluations by 1
      IncrementNumberEvaluations();

      // Record the responses back into the Design
      const RealVector& ftn_vals =
          sim_model.current_response().function_values();

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
  // now.  The call to sim_model.synchronize causes the program to block
  // until all the results are available.  We can then record the
  // responses in the same fashion as above.  Note that the linear
  // constraints have already been computed!!
  if(sim_model.asynch_flag()) {
    // Wait for the responses.
    const IntResponseMap& response_map = sim_model.synchronize();
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
  // Only run error_model when sim_model was successfull. Otherwise, we let
  // the code to fall through back to Dakota for printing errors messages.
  if(ret) {

    // Now we go over all true designs and map errors.
    IntRealVectorMap::const_iterator tr_it(
      decision_maker.GetBeginForTrueDatabase());
    const IntRealVectorMap::const_iterator tr_e(
      decision_maker.GetEndForTrueDatabase());
    for(; tr_it!=tr_e; ++tr_it) {

      // Note: We assume that the designs have already been evaluated. Furthermore, 
      // we do not check whether the function evaluation limit has been exceeded, 
      // since error evaluations are considered auxiliary and should not contribute
      //  to the user-specified maximum function evaluation count.
    
      // send this guy out for evaluation using the "error_model"

      // first, set the current values of the variables in the model
      ModelUtils::continuous_variables(error_model, tr_it->second);
      //ModelUtils::discrete_int_variables(error_model, disc_int_vars);
      //ModelUtils::discrete_real_variables(error_model, disc_real_vars);
      //// Strings set by calling single value setter for each
      ////for (size_t i=0; i<disc_string_vars.num_elements(); ++i)
      ////  ModelUtils::discrete_string_variable(sim_model, disc_string_vars[i],i);
      //// Could use discrete_string_varables to avoid overhead of repeated 
      //// function calls, but it takes a StringMultiArrayConstView, which
      //// must be created from disc_string_vars. Maybe there's a simpler way,
      //// but...
      //// const size_t &dsv_len = disc_string_vars.num_elements();
      //// StringMultiArrayConstView dsv_view = disc_string_vars[ 
      ////   boost::indices[idx_range(0,dsv_len)]];
      //// ModelUtils::discrete_string_variables(this->_model, dsv_view);
      
      //========================================================================
      // Pass error flag to analysis driver.
      //========================================================================
      SetStateVariables(tr_it->second, state_int_vars, state_string_vars, true); 
      ModelUtils::inactive_discrete_int_variables(error_model, state_int_vars);
      const size_t idsv_len = state_string_vars.num_elements();
      StringMultiArrayConstView idsv_view = state_string_vars[
        boost::indices[idx_range(0,idsv_len)]];
      ModelUtils::inactive_discrete_string_variables(error_model, idsv_view);

      // now request the evaluation in synchronous or asyncronous mode.
      if(error_model.asynch_flag()) 
        error_model.evaluate_nowait();
      else {
        // The following method call will use the default
        // Active set vector which is to just compute
        // function values, no gradients or hessians.
        error_model.evaluate();

        // Record the error responses from metadata
        const std::vector<RespMetadataT>& mtd_vals =
          error_model.current_response().metadata();
        const StringArray& mtd_labels = 
          error_model.current_response().shared_data().metadata_labels();

	      // Error response should not be sent to JEGA.
        RecordErrorInDecisionMaker(tr_it->first, mtd_vals, mtd_labels, tr_it->second);
      }

    }

    // If we did our evaluations asynchronously, we did not yet record
    // the results (because they were not available).  We need to do so
    // now.  The call to error_model.synchronize causes the program to block
    // until all the results are available.  We can then record the
    // responses in the same fashion as above.  Note that the linear
    // constraints have already been computed!!
    if(error_model.asynch_flag()) {
      // Wait for the responses.
      const IntResponseMap& response_map = error_model.synchronize();
      size_t num_resp = response_map.size();

      EDDY_ASSERT(num_resp == num_eval_reqs);

      // prepare to access the elements of the response_map by iterator.
      IntRespMCIter r_cit = response_map.begin();

      // Record the set of responses in the DesignGroup
      for(tr_it=decision_maker.GetBeginForTrueDatabase(); 
          tr_it!=tr_e; ++tr_it) {
        // Error response should not be sent to JEGA
        const std::vector<RespMetadataT>& mtd_vals =
          r_cit->second.metadata();
        const StringArray& mtd_labels = 
          r_cit->second.shared_data().metadata_labels();
        RecordErrorInDecisionMaker(tr_it->first, mtd_vals, mtd_labels, 
                                   tr_it->second);

        //increment
        ++r_cit;
      }
    }

  }

  TrainDecisionMaker();

  return ret;
}


//-----------------------------------------------------------------------------


void JegaEvaluator::TrainDecisionMaker()
{

}

//-----------------------------------------------------------------------------

