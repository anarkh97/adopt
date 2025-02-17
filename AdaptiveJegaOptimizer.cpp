#include<AdaptiveJegaOptimizer.h>

using namespace std;
using namespace Dakota;
using namespace JEGA::Logging;
using namespace JEGA::FrontEnd;
using namespace eddy::utilities;
using namespace JEGA::Utilities;
using namespace JEGA::Algorithms;

//-----------------------------------------------------------------------------

AdaptiveJegaOptimizer::AdaptiveJegaOptimizer(ProblemDescDB &prob_db, Model &model)
                     : Optimizer(prob_db, model, 
                                 std::shared_ptr<AdaptiveJegaTraits>(
                                 new AdaptiveJegaTraits())),
                       param_db(nullptr), eval_creator(nullptr)
{

  EDDY_FUNC_DEBUGSCOPE

  //! Exit if AdaptiveJega was called with moga.
  if(methodName == MOGA) {
    JEGALOG_II_G(this, text_entry(lfatal(), "Adaptive JEGA Error: Multi-objective "
                                            "optimization currently not supported."))
  }
  else if(methodName != SOGA) {
    JEGALOG_II_G(this, text_entry(lfatal(), "Adaptive JEGA Error: \"" + 
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

  LoadParameters();

  // population_size is extracted by JEGA in
  // GeneticAlgorithmInitializer::PollForParameters(), but it is
  // also needed here to specify the algorithmic concurrency.  Note
  // that the JEGA population size may grow or shrink during its
  // iterations, so this is only an initial estimate.
  int pop_size = this->prob_db.get_int("method.population_size");
  maxEvalConcurrency = *= pop_size;



}

//-----------------------------------------------------------------------------

//! Repurposed code from JEGAOptimizer::core_run()
void
AdaptiveJegaOptimizer::core_run()
{
    EDDY_FUNC_DEBUGSCOPE

    // Load up an algorithm config and a problem config.
    ProblemConfig pConfig;
    LoadProblemConfig(pConfig);

    // TODO: AlgorithmConfig aConfig(*this->_theEvalCreator, *this->_theParamDB);
    LoadAlgorithmConfig(aConfig);

    // retrieve parameter database for repeated use (is actaully _theParamDB)
    ParameterDatabase& pdb = aConfig.GetParameterDB();

    // Create a new driver for JEGA.
    // TODO: JEGAOptimizer::Driver driver(pConfig);

    // Get the algorithm separately (rather than simply running the current
    // configuration) in case we need to change the initializer.
    // TODO: GeneticAlgorithm* theGA = driver.ExtractAllData(aConfig);

    // Get the name of the GA for repeated use below.  We need this regardless
    // of whether or not logging b/c it is used in a fatal error.
    const string& name = theGA->GetName();

    // The initializer requires some additional logic to account for the
    // possibility that JEGA is being used in a Dakota strategy.  If that is
    // the case, the _initPts array will be non-empty and we will use them
    // instead of whatever initialization has been specified by the user.
    if(!this->_initPts.empty())
    {
        const GeneticAlgorithmInitializer& oldInit =
            theGA->GetOperatorSet().GetInitializer();

        JEGALOG_II_G(lquiet(), this,
            text_entry(lquiet(), name + ": discovered multiple initial "
                "points presumably supplied by a previous iterator in a "
                "strategy.  The \"" + oldInit.GetName() + "\" initializer "
                "will not be used and instead will be replaced with the "
                "double_matrix initializer which will read the supplied "
                "initial points."
                )
            )

        pdb.AddIntegralParam(
            "method.population_size", static_cast<int>(oldInit.GetSize())
            );

        pdb.AddDoubleMatrixParam(
            "method.jega.design_matrix", ToDoubleMatrix(initial_points())
            );

        GeneticAlgorithmInitializer* newInit =
            AllOperators::FullInstance().GetInitializer(
                "double_matrix", *theGA
                );

        JEGAIFLOG_II_G_F(newInit == 0x0, this,
            text_entry(lfatal(), name + ": Unable to resolve "
                "Initializer \"double_matrix\".")
            );

        JEGAIFLOG_II_F(!theGA->SetInitializer(newInit),
            theGA->GetLogger(), this,
            text_entry(lfatal(), name + ": Unable to set the initializer to "
                "double_matrix because it is incompatible with the other "
                "operators."
                )
            )

        JEGAIFLOG_II_F(
            !newInit->ExtractParameters(pdb), theGA->GetLogger(), this,
            text_entry(lfatal(),
                name + ": Failed to retrieve the parameters for \"" +
                newInit->GetName() + "\".")
            );

    }

    JEGALOG_II_G(lverbose(), this,
        text_entry(lverbose(),
            name + ": About to perform algorithm execution.")
            )

    // TODO: DesignOFSortSet bests(driver.PerformIterations(theGA));

    JEGALOG_II_G(lverbose(), this,
        ostream_entry(lverbose(), name + ": algorithm execution completed. ")
            << bests.size() << " solutions found. Passing them back to DAKOTA."
        )

    // Return up to numBest solutions to DAKOTA, sorted first by L2
    // constraint violation, then (utopia distance or weighted sum
    // objective value).  So the single "best" will be at the front.
    //
    // Load up to numBest solutions into the arrays of best responses
    // and variables.  If this is MOGA, the array will then contain
    // the Pareto set.  If it is SOGA, it will contain all the
    // solutions with the same best "fitness".
      
    // populate the sorted map of best solutions (fairly lightweight
    // map) key is pair<constraintViolation, fitness>, where fitness
    // is either utopia distance (MOGA) or objective value (SOGA)
    std::multimap<RealRealPair, Design*> designSortMap;
    this->GetBestSOSolutions(bests, *theGA, designSortMap);

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
    for( ; best_it != best_end; ++best_it, ++index)
    {
        this->LoadDakotaResponses(
            *(best_it->second),
            this->bestVariablesArray[index],
            this->bestResponseArray[index]
            );
    }

    // now we are done with our solution set so we can flush it
    // per Driver rules.
    bests.flush();

    JEGALOG_II_G(lquiet(), this,
        text_entry(lquiet(), name + ": find optimum completed and all "
            "results have been passed back to DAKOTA.\n\n")
        )

    // We can not destroy our GA.
    // TODO: driver.DestroyAlgorithm(theGA);
}

//-----------------------------------------------------------------------------

