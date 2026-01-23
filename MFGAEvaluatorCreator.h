#ifndef _MFGA_EVALUATOR_CREATOR_H_
#define _MFGA_EVALUATOR_CREATOR_H_

#include <../Utilities/include/JEGAConfig.hpp>
#include <../Utilities/include/Logging.hpp>

#include <GeneticAlgorithm.hpp>
#include <GeneticAlgorithmEvaluator.hpp>
#include <../FrontEnd/Core/include/EvaluatorCreator.hpp>

// Include the main optimizer class to get Model and other dakota things.
#include <DakotaOptimizer.hpp>

#include <AdaptiveDecisionMaker.h>

/************************************************
 * MFGAEvaluatorCreator is a utility class that 
 * recieves the sim_model and error_model from 
 * the AdaptiveJegaOptimizer and passes them to 
 * the Evaluator class. This is the standard
 * procedure required by JEGA::FrontEnd.
 ***********************************************/

namespace MultiFidelityOptimizer {

namespace detail {

class MFGAEvaluatorCreator : public JEGA::FrontEnd::EvaluatorCreator
{

  Dakota::Model         &sim_model;
  Dakota::Model         &error_model;
  AdaptiveDecisionMaker &decision_maker;

public:
  //! Constructor
  MFGAEvaluatorCreator(Dakota::Model &sim_model, Dakota::Model &error_model,
                       AdaptiveDecisionMaker &decision_maker);

  //! Destructor
  ~MFGAEvaluatorCreator() {};

  //! Called internally by JEGA libarary. Overriden to pass Dakota Model
  JEGA::Algorithms::GeneticAlgorithmEvaluator *
  CreateEvaluator(JEGA::Algorithms::GeneticAlgorithm &alg) override;
};

} // detail

} // MultiFidelityOptimizer

#endif
