#ifndef _JEGA_EVALUATOR_H_
#define _JEGA_EVALUATOR_H_

#include<../Utilities/include/JEGAConfig.hpp>
#include<../Utilities/include/JEGATypes.hpp>

// JEGA Core
#include<GeneticAlgorithm.hpp>
#include<GeneticAlgorithmEvaluator.hpp>

// JEGA utility includes.
#include<../Utilities/include/DesignGroup.hpp>

// Eddy utility includes.
#include<utilities/include/EDDY_DebugScope.hpp>

// Include the main optimizer class to get Model and other dakota things.
#include<DakotaOptimizer.hpp>

#include<AdaptiveDecisionMaker.h>

/*****************************************************
 * Custom implementation of the 
 * GeneticAlgorithmEvaluator class in JEGA::Algorithms.
 * This class mimics JEGAOptimizer::Evaluator class
 * but adds the Gaussian process model.
 ****************************************************/

class JegaEvaluator : public JEGA::Algorithms::GeneticAlgorithmEvaluator {

private:

  size_t switch_label_idx;

  Dakota::Model& sim_model;
  Dakota::Model& error_model;
  AdaptiveDecisionMaker& decision_maker;

public:

  static const std::string& Name() {
    EDDY_FUNC_DEBUGSCOPE
    static const std::string ret("Adaptive JEGA Evaluator");
    return ret;
  } 

  static const std::string& Description() {
    EDDY_FUNC_DEBUGSCOPE
    static const std::string ret(
        "This evaluator uses Sandia's DAKOTA optimization software to "
        "evaluate the passed in Designs.  This makes it possible to "
        "take advantage of the fact that DAKOTA is designed to run on "
        "massively parallel machines."
        );
    return ret;
  }

  std::string GetName() const override;

  std::string GetDescription() const override;

  JEGA::Algorithms::GeneticAlgorithmOperator* 
  Clone(JEGA::Algorithms::GeneticAlgorithm& algorithm) const override;

  //! This constructor should be used.
  JegaEvaluator(JEGA::Algorithms::GeneticAlgorithm& algorithm_,
                Dakota::Model& sim_model_, 
                Dakota::Model& error_model_,
                AdaptiveDecisionMaker& decision_maker_);

  //! Copy constructor
  JegaEvaluator(const JegaEvaluator& copy);

  //! Copy constructor w/ algorithm and model
  JegaEvaluator(const JegaEvaluator& copy, 
                JEGA::Algorithms::GeneticAlgorithm& algorithm_,
                Dakota::Model& sim_model_, 
                Dakota::Model& error_model_,
                AdaptiveDecisionMaker& decision_maker_);

  bool Evaluate(JEGA::Utilities::DesignGroup& group) override;

  //! This Evaluate is overriden so that an error is raised.
  //! According to JEGAOptimizer::Evaluator this function
  //! is not compatible w/ DAKOTA's multi-level parallelism.
  bool Evaluate(JEGA::Utilities::Design& des) override;

  //! Custom utility to run error model
  void ErrorEvaluationLoop();
  bool EvaluationLoop(JEGA::Utilities::DesignGroup& group, Dakota::String type);

protected:

  //! Here we only deal w/ continuous design variables.
  void SeparateVariables(const JEGA::Utilities::Design& from, 
                         Dakota::RealVector& into_cont) const;

  Dakota::String GetEvaluationType(const Dakota::RealVector& cont_vars);

  void SetStateVariables(const Dakota::RealVector& cont_vars,
                         Dakota::IntVector& into_disc_int,
			                   Dakota::StringMultiArray& into_disc_string,
                         Dakota::String type);

  void RecordResponses(const Dakota::RealVector& from, 
                       JEGA::Utilities::Design& into) const;

  //! Functions to update decision maker.
  void RecordEvaluationInDecisionMaker(const int id,
                                       const Dakota::RealVector& cont_vars,
                                       const Dakota::StringMultiArray& disc_strings);
  void RecordErrorInDecisionMaker(const int id,
                                  const Dakota::RealVector& from,
                                  const Dakota::RealVector& cont_vars);

  std::size_t GetNumberNonLinearConstraints() const;
  std::size_t GetNumberLinearConstraints() const;

private:

  void SetupLabelIndex();

  //! Should not be used.
  JegaEvaluator(JEGA::Algorithms::GeneticAlgorithm& algorithm);

};	

#endif
