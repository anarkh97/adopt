/************************************************************************
 * MODIFICATIONS IN THIS DERIVED FILE
 *
 * Derived from:
 *   - Dakota's JEGAOptimizer
 *   - Upstream commit: 0daaaa2237bd79445349c47d157a4a2b73db7452
 *
 * Modifications by Aditya Narkhede, 2025:
 *   - Renamed JEGAOptimizer to AdaptiveJegaOptimizer
 *   - Removed MOGA capability.
 *   - Included a seperate Dakota::Model for error calculations.
 *   - Pulled out EvaluatorCreator, Evaluator, and Driver
 *     sub-classes.
 *
 * Licensing:
 *   This file continues to be licensed under the GNU Lesser General Public
 *   License v2.1 or (at your option) any later version. See the bundled
 *   license file `LICENSE.LGPL-2.1` or
 *   <https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html>.
 ************************************************************************/

#include<AdaptiveJegaOptimizer.h>

// JEGAConfig.hpp should be the first include in all JEGA files.
#include<../Utilities/include/JEGAConfig.hpp>
#include<../Utilities/include/Logging.hpp>

// Standard includes.

// JEGA core includes.
#include<GeneticAlgorithm.hpp>
#include<GeneticAlgorithmEvaluator.hpp>
#include<GeneticAlgorithmInitializer.hpp>
#include<OperatorGroups/AllOperators.hpp>

// SOGA-specific API
#include<../SOGA/include/SOGA.hpp>

// JEGA front end includes.
#include<../FrontEnd/Core/include/Driver.hpp>
#include<../FrontEnd/Core/include/ProblemConfig.hpp>
#include<../FrontEnd/Core/include/AlgorithmConfig.hpp>
#include<../FrontEnd/Core/include/EvaluatorCreator.hpp>

// Dakota includes.
#include<ProblemDescDB.hpp>
#include<MarginalsCorrDistribution.hpp>
#include<model_utils.hpp>

// Eddy utility includes.
#include<utilities/include/EDDY_DebugScope.hpp>

// JEGA utility includes.
#include<../Utilities/include/DesignGroup.hpp>
#include<../Utilities/include/ConstraintInfo.hpp>
#include<../Utilities/include/ParameterExtractor.hpp>
#include<../Utilities/include/BasicParameterDatabaseImpl.hpp>
#include<../Utilities/include/MultiObjectiveStatistician.hpp>
#include<../Utilities/include/SingleObjectiveStatistician.hpp>

#include<algorithm>
#include<sstream>

using namespace std;
using namespace Dakota;
using namespace JEGA::Logging;
using namespace JEGA::FrontEnd;
using namespace eddy::utilities;
using namespace JEGA::Utilities;
using namespace JEGA::Algorithms;

//-----------------------------------------------------------------------------

template <typename T>
string asstring(const T& val)
{
  EDDY_FUNC_DEBUGSCOPE
  ostringstream ostr;
  ostr << val;
  return ostr.str();
}

//-----------------------------------------------------------------------------


AdaptiveJegaOptimizer::AdaptiveJegaOptimizer(ProblemDescDB &prob_db, 
                                             shared_ptr<Model> sim_model,
                                             shared_ptr<Model> err_model)
                     : Optimizer(prob_db, 
                                 sim_model,
                                 std::shared_ptr<AdaptiveJegaTraits>(
                                 new AdaptiveJegaTraits())),
                       param_db(nullptr), eval_creator(nullptr), 
                       decision_maker(nullptr)
{

  EDDY_FUNC_DEBUGSCOPE

  //! Make sure that the inactive view is set to STATE variables.
  //! NOTE: Here, iterated model can be a Scaling Model when
  //        user calls for scaling in the input file. Consequently,
  //        sim_model and iteratedModel do not point at the same thing.
  iteratedModel->inactive_view(MIXED_STATE);
  err_model->inactive_view(MIXED_STATE);

  //! Exit if AdaptiveJega was called with moga.
  if(methodName == Dakota::MOGA) {
    JEGALOG_II_G(lfatal(), this, text_entry(lfatal(), 
                 "Adaptive JEGA Error: Multi-objective "
                 "optimization currently not supported."))
  }
  else if(methodName != Dakota::SOGA) {
    JEGALOG_II_G(lfatal(), this, text_entry(lfatal(), 
                 "Adaptive JEGA Error: \"" + 
                 method_enum_to_string(methodName) + 
                 "\" is an invalid method specification."))
  }

  //! Initialize JEGA using the Front End Driver.
  if(!Driver::IsJEGAInitialized()) {

    int seed = prob_db.get_int("method.random_seed");

    // if user did not specify a seed pass "0" to JEGA
    size_t rand_seed = (seed < 0) ? 0 : (size_t)seed;

    // user output level
    short dakota_level = prob_db.get_short("method.output");

    LogLevel jega_level;
    switch(dakota_level) {

      case SILENT_OUTPUT : jega_level = lsilent(); break;
      case NORMAL_OUTPUT: jega_level = lnormal(); break;
      case DEBUG_OUTPUT: jega_level = ldebug(); break;
      case QUIET_OUTPUT: jega_level = lquiet(); break;
      case VERBOSE_OUTPUT: jega_level = lverbose(); break;
      default: jega_level = ldefault();

    }

    const bool jega_register_signals = false;
    Driver::InitializeJEGA(
        "JEGAGlobal.log", jega_level, rand_seed, Logger::ABORT,
        jega_register_signals
        );
  }

  // exit if initialization failed
  JEGAIFLOG_II_G_F(!Driver::IsJEGAInitialized(), this,
      text_entry(lfatal(), "JEGAOptimizer Error: Unable to initialize JEGA")
      );

  LoadParameterDatabase();

  // population_size is extracted by JEGA in
  // GeneticAlgorithmInitializer::PollForParameters(), but it is
  // also needed here to specify the algorithmic concurrency.  Note
  // that the JEGA population size may grow or shrink during its
  // iterations, so this is only an initial estimate.
  int pop_size = prob_db.get_int("method.population_size");
  maxEvalConcurrency *= pop_size;

  // Create the adaptive decision maker object
  // Chose to instantiate here so that we can pass
  // user I/O objects here.
  decision_maker = make_shared<AdaptiveDecisionMaker>(prob_db);

  // iteratedModel is the sim_model set internally 
  // by DakotaOptimizer
  eval_creator = make_shared<JegaEvaluatorCreator>(*iteratedModel, 
                                                   *err_model,
                                                   *decision_maker);
  
}

//-----------------------------------------------------------------------------

//! Smart pointers deleted automatically
AdaptiveJegaOptimizer::~AdaptiveJegaOptimizer()
{

}

//-----------------------------------------------------------------------------

//! Repurposed code from JEGAOptimizer::core_run()
void
AdaptiveJegaOptimizer::core_run()
{
  EDDY_FUNC_DEBUGSCOPE

  // Load up an algorithm config and a problem config.
  ProblemConfig prob_config;
  LoadProblemConfig(prob_config);

  AlgorithmConfig alg_config(*eval_creator, *param_db);
  LoadAlgorithmConfig(alg_config);

  // retrieve parameter database for repeated use (is actaully param_db)
  ParameterDatabase& pdb = alg_config.GetParameterDB();

  // Create a new driver for AdaptiveJega.
  JegaDriver driver(prob_config);

  // Get the algorithm separately (rather than simply running the current
  // configuration) in case we need to change the initializer.
  GeneticAlgorithm* ga_algorithm = driver.ExtractAllData(alg_config);

  // Get the name of the GA for repeated use below.  We need this regardless
  // of whether or not logging b/c it is used in a fatal error.
  const string& name = ga_algorithm->GetName();

  JEGALOG_II_G(lverbose(), this, text_entry(lverbose(), name + 
               ": About to perform algorithm execution."))

  // For SOGA this multiset contains designs with same "fitness".
  DesignOFSortSet bests(driver.PerformIterations(ga_algorithm)); 

  JEGALOG_II_G(lverbose(), this, ostream_entry(lverbose(), name + 
               ": algorithm execution completed. ") << bests.size() 
               << " solutions found. Passing them back to DAKOTA.")

  // Return up to numBest solutions to DAKOTA, sorted first by L2
  // constraint violation, then (utopia distance or weighted sum
  // objective value).  So the single "best" will be at the front.
  
  // For SOGA, we have to return solution with the same best "fitness".
  // MOGA is currently not supported.  
	
  // populate the sorted map of best solutions (fairly lightweight
  // map) key is pair<constraint_violation, fitness>, where fitness
  // is either utopia distance (MOGA) or objective value (SOGA)
  std::multimap<RealRealPair, Design*> designSortMap;
  GetBestSOSolutions(bests, *ga_algorithm, designSortMap);

  JEGAIFLOG_II_G(designSortMap.size() == 0, lquiet(), this,
      text_entry(lquiet(), name + ": was unable to identify at least one "
          "best solution.  The Dakota best variables and best responses "
          "objects will be empty.\n\n")
      )

  // load the map into the DAKOTA vectors
  resize_best_resp_array(designSortMap.size());
  resize_best_vars_array(designSortMap.size());
  
  std::multimap<RealRealPair, Design*>::const_iterator best_it = 
      designSortMap.begin(); 
  const std::multimap<RealRealPair, Design*>::const_iterator best_end = 
      designSortMap.end(); 
  ResponseArray::size_type index = 0;
  for( ; best_it != best_end; ++best_it, ++index) {
    LoadDakotaResponses(*(best_it->second), bestVariablesArray[index],
                        bestResponseArray[index]);
  }

  // now we are done with our solution set so we can flush it
  // per Driver rules.
  bests.flush();

  JEGALOG_II_G(lquiet(), this, text_entry(lquiet(), name + 
               ": find optimum completed and all results have been passed "
               "back to DAKOTA.\n\n"))

  // We can not destroy our GA --- Let JEGA handle this
  driver.DestroyAlgorithm(ga_algorithm);

}

//-----------------------------------------------------------------------------

void AdaptiveJegaOptimizer::LoadParameterDatabase()
{
  EDDY_FUNC_DEBUGSCOPE

  ProblemDescDB &prob_db = probDescDB;

  // create a new parameter database
  if(param_db) param_db.reset();
  param_db = shared_ptr<BasicParameterDatabaseImpl>(
    new BasicParameterDatabaseImpl()		  
  );

  // Duplicate in all the integral parameters.
  const int& random_seed = prob_db.get_int("method.random_seed");
  if (random_seed != 0)
    param_db->AddIntegralParam("method.random_seed", random_seed);
  else
    param_db->AddIntegralParam("method.random_seed", -1);

  // now all the reals
  const Real& constraint_penalty = prob_db.get_real("method.constraint_penalty");
  // HARD CODED as there is some bug in DAKOTA.
  if (constraint_penalty >= 0.)
    //param_db->AddDoubleParam("method.constraint_penalty", constraint_penalty);
    param_db->AddDoubleParam("method.constraint_penalty", 10.0);
  else
    param_db->AddDoubleParam("method.constraint_penalty", 10.0);

  const Real& crossover_rate = prob_db.get_real("method.crossover_rate");
  if (crossover_rate >= 0.)
    param_db->AddDoubleParam("method.crossover_rate", crossover_rate);
  else
    param_db->AddDoubleParam("method.crossover_rate", 0.75);

  if (prob_db.get_string("method.mutation_type") != "")
    param_db->AddDoubleParam("method.mutation_rate", prob_db.get_real("method.mutation_rate"));
  else
    param_db->AddDoubleParam("method.mutation_rate", 0.1);

  param_db->AddDoubleParam("method.mutation_scale", prob_db.get_real("method.mutation_scale"));

  double perc_change = prob_db.get_real("method.jega.percent_change");
  param_db->AddDoubleParam("method.jega.percent_change",
                           (perc_change < 0) ? 1.0e-4 : perc_change);

  double conv_tol = prob_db.get_real("method.convergence_tolerance");
  param_db->AddDoubleParam("method.convergence_tolerance",
                           (conv_tol < 0) ? 1.0e-4 : conv_tol);

  param_db->AddDoubleParam("method.jega.shrinkage_percentage",
                           prob_db.get_real("method.jega.shrinkage_percentage"));

  param_db->AddDoubleParam("method.jega.fitness_limit",
                           prob_db.get_real("method.jega.fitness_limit"));

  // now get all the size_t's
  param_db->AddSizeTypeParam("method.jega.num_cross_points",
                             prob_db.get_sizet("method.jega.num_cross_points"));

  param_db->AddSizeTypeParam("method.jega.num_parents",
                             prob_db.get_sizet("method.jega.num_parents"));

  param_db->AddSizeTypeParam("method.jega.num_offspring",
                             prob_db.get_sizet("method.jega.num_offspring"));

  param_db->AddSizeTypeParam("method.jega.num_generations",
                             prob_db.get_sizet("method.jega.num_generations"));

  param_db->AddSizeTypeParam("method.jega.max_designs",
                             prob_db.get_sizet("method.jega.num_designs"));

  // Note that the population size, max evals, and max gens are in as ints.
  // Do a conversion for each here.
  param_db->AddSizeTypeParam("method.population_size",
                             (size_t)(prob_db.get_int("method.population_size")));

  param_db->AddSizeTypeParam("method.max_iterations",
                             (size_t)(prob_db.get_sizet("method.max_iterations")));

  param_db->AddSizeTypeParam("method.max_function_evaluations",
                             (size_t)(
                             prob_db.get_sizet("method.max_function_evaluations")));

  // Dakota does not currently expose the input that indicates that there is
  // a minimum allowable population size for the below limit selector.  Add
  // a default explicitly here to prevent any problems.
  param_db->AddSizeTypeParam("method.jega.minimum_selections", 2);

  // Dakota does not currently expose the evaluation concurrency flag
  // through the interface nor should it b/c Dakota handles evaluations.
  param_db->AddSizeTypeParam("method.jega.eval_concurrency", 1);

  // Now get all the booleans
  param_db->AddBooleanParam("method.print_each_pop",
                            prob_db.get_bool("method.print_each_pop"));

  // Dakota does not currently expose the flag to instruct the GA whether or
  // not to write the final data file and discards file.  Put those in here
  // to avoid warnings about them missing.
  param_db->AddBooleanParam("method.print_final_data", true);
  param_db->AddBooleanParam("method.print_discards", true);

  // now get all the strings.
  
  // Dakota does not expose the ability to specify the weighted sum
  // only fitness assessor b/c it is only for use with the favor
  // feasible selector.  Likewise, the favor feasible selector can
  // only be used with the weighted sum fitness asessor.  Because of
  // this, we will detect use of the favor feasible and enforce
  // the use of the weighted sum only.  We will write a log message
  // about it.
  const string& selector = prob_db.get_string("method.replacement_type");
  if (selector != "")
    param_db->AddStringParam("method.replacement_type", selector);
  else if (methodName == Dakota::SOGA)
    param_db->AddStringParam("method.replacement_type", "elitist");
  else if (methodName == Dakota::MOGA) {
    JEGALOG_II_G(lfatal(), this, text_entry(lfatal(), 
                 "Adaptive JEGA Error: Multi-objective "
                 "optimization currently not supported."))
  }
  
  const string& fitness = prob_db.get_string("method.fitness_type");
  if(selector == "favor_feasible") {
    JEGALOG_II_G(lquiet(), this, text_entry(lquiet(),
                 "Use of the favor_feasible selector has been detected.  Use of "
                 "this selector type requires use of the \"weighted_sum_only\" "
                 "fitness assessor.  Therefore, use of the \"" + fitness +
                 "\" will be changed to use of the \"weighted_sum_only\"."))

    param_db->AddStringParam("method.fitness_type", "weighted_sum_only");
  }
  else {
    if (fitness != "")
      param_db->AddStringParam("method.fitness_type", fitness);
    else if (methodName == Dakota::SOGA)
      param_db->AddStringParam("method.fitness_type", "merit_function");
    else if (methodName == Dakota::MOGA) {
      JEGALOG_II_G(lfatal(), this, text_entry(lfatal(), 
                   "Adaptive JEGA Error: Multi-objective "
                   "optimization currently not supported."))
    }
  }

  const string& crossover_operator = prob_db.get_string("method.crossover_type");
  if (crossover_operator != "")
    param_db->AddStringParam("method.crossover_type", crossover_operator);
  else
    param_db->AddStringParam("method.crossover_type", "shuffle_random");

  const string& mutation_operator = prob_db.get_string("method.mutation_type");
  if (mutation_operator != "")
    param_db->AddStringParam("method.mutation_type", mutation_operator);
  else
    param_db->AddStringParam("method.mutation_type", "replace_uniform");

  param_db->AddIntegralParam("method.output", prob_db.get_short("method.output"));

  param_db->AddStringParam("method.initialization_type",
                           prob_db.get_string("method.initialization_type"));

  param_db->AddStringParam("method.flat_file",
                           prob_db.get_string("method.flat_file"));

  // Dakota does not currently expose the input that allows one to specify
  // the location of any data files written by JEGA.  Use the default but
  // specify it explicitly here to prevent any problems.
  param_db->AddStringParam("method.jega.data_directory", "./");

  // Dakota does not currently expose the final data file name pattern
  // through the interface.
  param_db->AddStringParam("method.jega.final_data_filename", "finaldata#.dat");

  // Dakota does not currently expose the main loop operator selection
  // through the interface.
  param_db->AddStringParam("method.jega.mainloop_type", "duplicate_free");

  // The log file gets special attention.  If it is the default global log
  // file name, we will replace it with an empty string b/c we don't want the
  // created GA to think it owns the global log file.
  string log_file = prob_db.get_string("method.log_file");
  param_db->AddStringParam("method.log_file", 
                           log_file == "JEGAGlobal.log" ? "" : log_file);
  const string& convergence_operator = prob_db.get_string("method.jega.convergence_type");

  if (convergence_operator != "")
    param_db->AddStringParam("method.jega.convergence_type", convergence_operator);
  else if (this->methodName == Dakota::SOGA)
    param_db->AddStringParam("method.jega.convergence_type", "average_fitness_tracker");
  else if (this->methodName == Dakota::MOGA) {
    JEGALOG_II_G(lfatal(), this, text_entry(lfatal(), 
                 "Adaptive JEGA Error: Multi-objective "
                 "optimization currently not supported."))
  }

  //!AN: This is required by JEGA
  param_db->AddStringParam("method.jega.niching_type",
                           prob_db.get_string("method.jega.niching_type"));

  param_db->AddStringParam("method.jega.postprocessor_type",
                           prob_db.get_string("method.jega.postprocessor_type"));

  // Dakota does not expose a flat file delimiter for the case where we
  // are using the flat file initializer but the initializer is going to
  // be looking for one if it is in use.
  param_db->AddStringParam("method.jega.initializer_delimiter", "");

  // now get all vector of doubles.
  const RealVector *dak_rv = &prob_db.get_rv("method.jega.niche_vector");

  JEGA::DoubleVector niche_vector(dak_rv->values(), 
                                  dak_rv->values() + dak_rv->length());

  param_db->AddDoubleVectorParam("method.jega.niche_vector",
                                 niche_vector);

  dak_rv = &prob_db.get_rv("method.jega.distance_vector");
  JEGA::DoubleVector distance_vector(dak_rv->values(),
                                     dak_rv->values() + dak_rv->length());

  param_db->AddDoubleVectorParam("method.jega.distance_vector", distance_vector);
 
  // when recasting is active, the weights may be transformed; get off Model
  dak_rv = &iteratedModel->primary_response_fn_weights();
  JEGA::DoubleVector mow_vector(dak_rv->values(), 
                                dak_rv->values() + dak_rv->length());

  param_db->AddDoubleVectorParam("responses.multi_objective_weights",
                                 mow_vector);

  // Dakota does not expose a capability to enter multiple flat file names
  // as a vector (can delimit them in the single string and JEGA will parse
  // them).  We will add an empty vector to prevent the warning from JEGA.
  param_db->AddStringVectorParam("method.flat_files", 
                                 std::vector<std::string>());



  // now get all int vector params.
  // nothing to do here.

  // now get all the double matrices.
  // nothing to do here.

  // now get all integer list params.
  // nothing to do here.

  // now get all string vectors
  // nothing to do here.

  // now get all string lists
  // nothing to do here.


}

//-----------------------------------------------------------------------------

void
AdaptiveJegaOptimizer::LoadAlgorithmConfig(AlgorithmConfig &config)
{
  EDDY_FUNC_DEBUGSCOPE

  ParameterDatabase& pdb = config.GetParameterDB();
  
  // Determine what kind of algorithm we are creating (MOGA or SOGA)
  // based on the methodName base class variable.
  AlgorithmConfig::AlgType alg_type;
  
  if(methodName == Dakota::MOGA)
    JEGALOG_II_G(lfatal(), this, text_entry(lfatal(), 
                 "Adaptive JEGA Error: Multi-objective "
                 "optimization currently not supported."))
  else if(methodName == Dakota::SOGA)
    alg_type = AlgorithmConfig::SOGA;
  else
      JEGALOG_II_G_F(this, text_entry(lfatal(), "Adaptive JEGA Error: \"" +
                     method_enum_to_string(methodName) + "\" is an invalid "
                     "method specification."))
  
  config.SetAlgorithmType(alg_type);
  
  // We will use the method id as the algorithm name if it is non-empty and
  // we will otherwise use the method name.
  config.SetAlgorithmName(method_id().empty() ?
                          method_enum_to_string(methodName) : 
                          method_id());
  
}

//-----------------------------------------------------------------------------

void
AdaptiveJegaOptimizer::LoadProblemConfig(ProblemConfig &config)
{
  EDDY_FUNC_DEBUGSCOPE

  LoadDesignVariables(config);
  LoadObjectiveFunctions(config);
  LoadConstraints(config);
}

//-----------------------------------------------------------------------------

void
AdaptiveJegaOptimizer::LoadDesignVariables(ProblemConfig &config)
{
  EDDY_FUNC_DEBUGSCOPE

  Model &iterated_model = *iteratedModel;

  // Currently, we only deal with continuous real design variables.
  // To simply the implementation we only load continuous real variables
  // into JEGA and raise errors when other variables are encountered.

  if(numDiscreteIntVars > 0) {
    JEGALOG_II_G_F(this, text_entry(lfatal(), "Adaptive JEGA Error: "
                   "Currently, Discrete Integer design variables are not "
                   "supported.\n"))
  }
  if(numDiscreteRealVars > 0) {
    JEGALOG_II_G_F(this, text_entry(lfatal(), "Adaptive JEGA Error: "
                   "Currently, Discrete Real design variables are not "
                   "supported.\n"))
  }
  if(numDiscreteStringVars > 0) {
    JEGALOG_II_G_F(this, text_entry(lfatal(), "Adaptive JEGA Error: "
                   "Currently, Discrete String design variables are not "
                   "supported.\n"))
  }

  const RealVector& clbs = ModelUtils::continuous_lower_bounds(iterated_model);
  const RealVector& cubs = ModelUtils::continuous_upper_bounds(iterated_model);
  StringMultiArrayConstView clabels = 
    iterated_model.current_variables().continuous_variable_labels();
  for(int i=0; i<numContinuousVars; ++i)
    config.AddContinuumRealVariable(clabels[i], clbs[i], cubs[i], 6);
  
  // Now make sure that an info was created for each variable.
  EDDY_ASSERT(config.GetDesignTarget().GetNDV() == numContinuousVars);

}

//-----------------------------------------------------------------------------

void
AdaptiveJegaOptimizer::LoadObjectiveFunctions(ProblemConfig &config)
{
  EDDY_FUNC_DEBUGSCOPE
  
  Model &iterated_model = *iteratedModel;

  if(numObjectiveFns > 1) {
    JEGALOG_II_G_F(this, text_entry(lfatal(), "Adaptive JEGA Error: "
                   "Currently, Multi-objective Design optimization "
                   "studies are not supported.\n"))
  }

  // For now, all objectives will be of type minimize.  Hopefully,
  // Dakota will soon support mixed extremization schemes.
  // Dakota does not support labeling objectives.  Until it does,
  // we will create a label that looks like "Nature Type Index".
  const StringArray&  labels = ModelUtils::response_labels(iterated_model);
  const BoolDeque& max_sense = iterated_model.primary_response_fn_sense();
  bool use_sense = !max_sense.empty();
  for(size_t i=0; i<numObjectiveFns; ++i)
    if (use_sense && max_sense[i])
      config.AddNonlinearMaximizeObjective("Non-Linear Maximize " + labels[i]);
    else
      config.AddNonlinearMinimizeObjective("Non-Linear Minimize " + labels[i]);
  
  // see to it that the numbers match up.
  EDDY_ASSERT(config.GetDesignTarget().GetNOF() == numObjectiveFns);
}

//-----------------------------------------------------------------------------

void
AdaptiveJegaOptimizer::LoadConstraints(ProblemConfig &config)
{
  EDDY_FUNC_DEBUGSCOPE
  
  // The information needed to create the constraint infos
  // is contained in the data structures of the base classes.
  // In particular, the Model (iteratedModel) has most of the
  // info. 
  const Model &m = *iteratedModel;
  
  /**************************************************************************
  
  Note the order in which these are created.  Do not change this order.  It
  is this way because of the order in which responses are returned out of the
  Model.  Technically, it only involves the first two blocks which create the
  non-linear constraints.  But don't mess with any of it anyway.
  
  **************************************************************************/
  
  
  // start with non-linear (2-sided) inequality constraints.
  // The information we need for these is in
  // nonlinear_ineq_constraint_lower_bounds and
  // nonlinear_ineq_constraint_upper_bounds.  As with objective
  // functions, Dakota does not allow labeling of constraints.
  // we will create a label that looks like "Nature Type Index".
  const RealVector& nln_ineq_lwr_bnds = 
    ModelUtils::nonlinear_ineq_constraint_lower_bounds(m);
  const RealVector& nln_ineq_upr_bnds =
    ModelUtils::nonlinear_ineq_constraint_upper_bounds(m);

//PDH: Dakota nonlinear constraints to JEGA nonlinear constraints.
//     Don't know what the JEGA data structure is.  These are all
//     mapped one entry at a time.
//     Looks like we don't have to worry about (JEGA) order.  Need to
//     determine if they have to be two-sided for JEGA.

  // Loop over all two sided non linear inequality constraitns and add an
  // info object for each.
  for(size_t i=0; i<numNonlinearIneqConstraints; ++i)
    config.AddNonlinearTwoSidedInequalityConstraint("Non-Linear Two-Sided "
                                                    "Inequality " + asstring(i),
                                                    nln_ineq_lwr_bnds[i], 
                                                    nln_ineq_upr_bnds[i]);
  
  // now do non-linear equality constraints.  The information we need for
  // these is in nonlinear_eq_constraint_targets.
  const RealVector& nln_eq_targets = 
    ModelUtils::nonlinear_eq_constraint_targets(m);
  for(size_t i=0; i<numNonlinearEqConstraints; ++i)
    config.AddNonlinearEqualityConstraint("Non-Linear Equality " + asstring(i),
                                          nln_eq_targets[i]);

//PDH: Dakota linear constraints to JEGA linear constraints.
//     Don't know what the JEGA data structure is.  These are all
//     mapped one entry at a time.
//     Looks like we don't have to worry about (JEGA) order.  Need to
//     determine if they have to be two-sided for JEGA.

  // now do linear (2-sided) inequality constraints  The information we need
  // for these is in linear_ineq_constraint_lower_bounds and
  // linear_ineq_constraint_upper_bounds.
  // In addition to bounds, these accept coefficients for possible shortcut
  // evaluation.  That information is in linear_ineq_constraint_coeffs.
  const RealVector& lin_ineq_lwr_bnds
      = ModelUtils::linear_ineq_constraint_lower_bounds(m);
  const RealVector& lin_ineq_upr_bnds
      = ModelUtils::linear_ineq_constraint_upper_bounds(m);
  const RealMatrix& lin_ineq_coeffs
      = ModelUtils::linear_ineq_constraint_coeffs(m);
  
  JEGA::DoubleVector lin_ineq_coeffs_row(lin_ineq_coeffs.numCols());

//PDH: RealMatrix -> set of std::vector
//     Just need the individual rows.  Check copy_row_vector to see if
//     transpose is also needed.

  for(size_t i=0; i<numLinearIneqConstraints; ++i) {
    copy_row_vector(lin_ineq_coeffs, i, lin_ineq_coeffs_row);
    
    config.AddLinearTwoSidedInequalityConstraint("Linear Two-Sided Inequality " 
                                                 + asstring(i), 
                                                 lin_ineq_lwr_bnds[i],
                                                 lin_ineq_upr_bnds[i],
                                                 lin_ineq_coeffs_row);
  }

  // now do linear equality constraints.  The information we need for these
  // is in lin_eq_targets. In addition to targets, these accept coefficients
  // for possible shortcut evaluation.  That information is in
  // linear_eq_constraint_coeffs.
  const RealVector& lin_eq_targets = ModelUtils::linear_eq_constraint_targets(m);
  const RealMatrix& lin_eq_coeffs = ModelUtils::linear_eq_constraint_coeffs(m);

  JEGA::DoubleVector lin_eq_coeffs_row(lin_eq_coeffs.numCols());

//PDH: RealMatrix -> set of std::vector
//     Just need the individual rows.  Check copy_row_vector to see if
//     transpose is also needed.

  for(size_t i=0; i<numLinearEqConstraints; ++i) {
    copy_row_vector(lin_eq_coeffs, i, lin_eq_coeffs_row);

    config.AddLinearEqualityConstraint("Linear Equality " + asstring(i),
                                       lin_eq_targets[i], 0.0, 
				       lin_eq_coeffs_row);
  }

  // see to it that the numbers match up.
  EDDY_ASSERT(pConfig.GetDesignTarget().GetNCN() == (numNonlinearIneqConstraints
              + numLinearIneqConstraints + numNonlinearEqConstraints 
              + numLinearEqConstraints));

}

//-----------------------------------------------------------------------------

//! Copied from JEGAOptimizer::GetBestSOSolutions()
void
AdaptiveJegaOptimizer::GetBestSOSolutions(const DesignOFSortSet &from, 
                                          const GeneticAlgorithm &ga,
                                          multimap<RealRealPair, Design*> &designSortMap)
{

  EDDY_FUNC_DEBUGSCOPE

  if(from.empty()) return;

  const DesignTarget& target = from.front()->GetDesignTarget();
  const ConstraintInfoVector& constraint_info = target.GetConstraintInfos();

  // get number of objective functions
  const eddy::utilities::uint64_t nof = target.GetNOF();

  EDDY_ASSERT(nof == 1)

  // get total number of constraints (nonlinear and linear)
  const eddy::utilities::uint64_t noc = target.GetNCN();

  // in order to order the points, need the weights; get them from
  // the GA to ensure solver/final results consistency


  //AN: Weights are manually fixed here.
  //    Was getting a segmentation fault here. Don't know why.

  JEGA::DoubleVector weights(nof, 1.0);
/*  
  try {
    const JEGA::Algorithms::SOGA& ga = 
      dynamic_cast<const JEGA::Algorithms::SOGA&>(ga);
    weights = ga.GetWeights();
  } catch(const std::bad_cast& bc_except) {
    Cerr << "\nError: could not cast GeneticAlgorithm to SOGA; exception:\n" 
         << bc_except.what() << std::endl;
    abort_handler(-1);
  }
*/

  // iterate the designs and sort first by constraint violation,
  // then objective function
  DesignOFSortSet::const_iterator design_it(from.begin());
  const DesignOFSortSet::const_iterator design_end(from.end());

  for(; design_it != design_end; ++design_it) {
    // L2 constraint violation for this design
    double constraint_violation = 0.0;

    for(size_t i=0; i<noc; ++i)
      constraint_violation +=
          Math::Pow(constraint_info[i]->GetViolationAmount(**design_it), 2);

    // Multi-objective sum for this Design over objective functions.
    // In the single objective case we can store
    // objective even if there's a constraint violation.
    double objectiveFunction = SingleObjectiveStatistician::ComputeWeightedSum(
      **design_it, weights);

    // insert the design into the map, keeping only numBest
    RealRealPair metrics(constraint_violation, objectiveFunction);

    if(designSortMap.size() < numFinalSolutions)
      designSortMap.insert(std::make_pair(metrics, *design_it));
    else {
      // if this Design is better than the worst, remove worst
      // and insert this one
      std::multimap<RealRealPair, Design*>::iterator worst_it =
          --designSortMap.end();

      if(metrics < worst_it->first) {
        designSortMap.erase(worst_it);
        designSortMap.insert(std::make_pair(metrics, *design_it));
      }
    }
  }

}

//-----------------------------------------------------------------------------

//! Copied from JEGAOptimizer::LoadDakotaResponses()
void
AdaptiveJegaOptimizer::LoadDakotaResponses(const Design &from, Variables &vars, 
                                           Response &resp) const
{

  size_t num_cv = vars.cv();
  size_t num_idiv = vars.idiv();
  size_t num_idsv = vars.idsv();

  RealVector c_vars(num_cv);
  IntVector d_i_vars(num_idiv);
  StringMultiArray d_s_vars(boost::extents[num_idsv]);

//PDH: JEGA variables to Dakota variables.
//     Don't know what the JEGA data structure is.  These are all
//     mapped one entry at a time.

    // The first numContinuousVars of a design will be all the continuous
    // variables of the problem (see LoadTheDesignVariables).
  for(size_t i=0; i<numContinuousVars; ++i)
    c_vars[i] = from.GetVariableValue(i);

  //! Get information from the decision maker
  
  StringMultiArrayConstView d_s_labels = 
    vars.inactive_discrete_string_variable_labels();

  size_t index=-1;
  size_t counter=0;
  for(const auto &label : d_s_labels) {
    if(label == "SWITCH") {
      index = counter;
      break;
    }
    counter++;
  }

  decision_maker->GetEvalTypeAndMetaData(c_vars, d_s_vars[index], 
    d_i_vars, num_idiv-1);

  size_t idsv_len = d_s_vars.num_elements();
  StringMultiArrayConstView idsv_view = d_s_vars[
    boost::indices[idx_range(0,idsv_len)]];

  vars.continuous_variables(c_vars);
  vars.inactive_discrete_int_variables(d_i_vars);
  vars.inactive_discrete_string_variables(idsv_view);

//PDH: JEGA responses to Dakota responses.
//     Don't know what the JEGA data structure is.  These are all
//     mapped one entry at a time.
//     Need to respect constraint ordering.

  // BMA TODO: Could always populate constraints and just get
  // primary responses from the DB, as in SNLL
  RealVector fn_vals(resp.num_functions());
  if (!localObjectiveRecast) {  // else local_recast_retrieve
    for(size_t i=0; i<numObjectiveFns; i++)
      fn_vals[i]= from.GetObjective(i);
  }

  // JEGA constraint ordering is nonlinear inequality, nonlinear equality,
  // linear inequality, linear equality
  // (see JEGAOptimizer::LoadTheConstraints()).
  for(size_t i=0; i<numNonlinearConstraints; ++i)
    fn_vals[i+numUserPrimaryFns] = from.GetConstraint(i);

  resp.function_values(fn_vals);

}

//-----------------------------------------------------------------------------

