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

#include<DakotaOptimizer.hpp>

/*****************************************************
 * Custom implementation of the 
 * GeneticAlgorithmEvaluator class in JEGA::Algorithms.
 * This class mimics JEGAOptimizer::Evaluator class
 * but adds the Gaussian process model.
 ****************************************************/

class JegaEvaluator : public JEGA::Algorithms::GeneticAlgorithmEvaluator {

private:

  Dakota::Model& model;

public:

  static const std::string& Name() {
    EDDY_FUNC_DEBUGSCOPE
    static const string ret("ADAPTIVE JEGA Evaluator");
    return ret;
  }

  static const std::string& Description() {
    EDDY_FUNC_DEBUGSCOPE
    static const string ret(
        "This evaluator uses Sandia's DAKOTA optimization software to "
        "evaluate the passed in Designs.  This makes it possible to "
        "take advantage of the fact that DAKOTA is designed to run on "
        "massively parallel machines."
        );
    return ret;
  }

  //! This Evaluate is overriden so that an error is raised.
  //! According to JEGAOptimizer::Evaluator this function
  //! is not compatible w/ DAKOTA's multi-level parallelism.
  bool Evaluate(JEGA::Utilities::Design& des) override {
    EDDY_FUNC_DEBUGSCOPE
  
    JEGALOG_II_F(GetLogger(), this, text_entry(lfatal(), GetName() + 
                 ": You cannot use Evaluate(Design&) with this "
                 "evaluator...ever."))
    return false;
  }

  std::string GetName() const override {
    EDDY_FUNC_DEBUGSCOPE
    return JegaEvaluator::Name();
  }

  std::string GetDescription() const override {
    EDDY_FUNC_DEBUGSCOPE
    return JegaEvaluator::Description();
  }

  JEGA::Algorithms::GeneticAlgorithmOperator* 
  Clone(JEGA::Algorithms::GeneticAlgorithm& algorithm) const override {
    EDDY_FUNC_DEBUGSCOPE
    return new JegaEvaluator(*this, algorithm, _model);
  }

protected:

  //! Here we only deal w/ continuous design variables.
  void SeparateVariables(const Design& from, RealVector& intoCont) const;

  void RecordResponses(const RealVector& from, Design& into) const;

  std::size_t GetNumberNonLinearConstraints() const {
    EDDY_FUNC_DEBUGSCOPE
    return _model.num_nonlinear_eq_constraints() +
           _model.num_nonlinear_ineq_constraints();
  }

  std::size_t GetNumberLinearConstraints() const {
    EDDY_FUNC_DEBUGSCOPE
    return _model.num_linear_eq_constraints() +
           _model.num_linear_ineq_constraints();
  }

public:

  //! This constructor should be used.
  JegaEvaluator(JEGA::Algorithms::GeneticAlgorithm& algorithm, 
                Dakota::Model& model);

  //! Copy constructor
  JegaEvaluator(const JegaEvaluator& copy);

  //! Copy constructor w/ algorithm and model
  JegaEvaluator(const JegaEvaluator& copy, 
                JEGA::Algorithms::GeneticAlgorithm& algorithm,
                Dakota::Model& model);

  bool Evaluate(JEGA::Utilities::DesignGroup& group) override;

private:

  //! Should not be used.
  JegaEvaluator(JEGA::Algorithms::GeneticAlgorithm& algorithm);

};	

#endif
