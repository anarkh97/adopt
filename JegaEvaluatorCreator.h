#ifndef _JEGA_EVALUATOR_CREATOR_H_
#define _JEGA_EVALUATOR_CREATOR_H_

#include<../Utilities/include/JEGAConfig.hpp>
#include<../Utilities/include/Logging.hpp>

#include<GeneticAlgorithm.hpp>
#include<GeneticAlgorithmEvaluator.hpp>
#include<../FrontEnd/Core/include/EvaluatorCreator.hpp>

// Include the main optimizer class to get Model and other dakota things.
#include<DakotaOptimizer.hpp>

class JegaEvaluatorCreator : public JEGA::FrontEnd::EvaluatorCreator {

  Dakota::Model &true_model;
  Dakota::Model &error_model;

public:

  //! Constructor
  JegaEvaluatorCreator(Dakota::Model& true_model, Dakota::Model& error_model); 

  //! Destructor
  ~JegaEvaluatorCreator() { };

  //! Called internally by JEGA libarary. Overriden to pass Dakota Model
  JEGA::Algorithms::GeneticAlgorithmEvaluator*
  CreateEvaluator(JEGA::Algorithms::GeneticAlgorithm& alg) override;

};

#endif
