#ifndef _ADAPTIVE_JEGA_OPTIMIZER_H_
#define _ADAPTIVE_JEGA_OPTIMIZER_H_

#include<JEGAOptimizer.hpp>

/*****************************************************
 *  This class extends JEGAOptimizer, Dakota's 
 *  implementation of John Eddy's Genetic 
 *  Algorithm (JEGA) library. We modify its optimization 
 *  loop to integrate an Adaptive Decision Model based 
 *  on Gaussian Process. This model predicts whether 
 *  high-fidelity fluid dynamics simulations (M2C) 
 *  are necessary or if simpler interpolation estimates 
 *  (FEST) will suffice. Note that JEGAOptimizer's 
 *  Evaluate, Driver, and EvaluatorCreator classes are 
 *  private and inaccessible here.
 *****************************************************/

class AdaptiveJegaOptimizer : public Dakota::JEGAOptimizer {

public:

  //! Default constructor
  AdaptiveJegaOptimizer(ProblemDescDB& prob_db, Model& model);

  //! Default Destructor
  ~AdaptiveJegaOptimizer() { };

  //! Initializes population
  void initialize_run() override;

  //! Interface to Dakota's core_run
  void core_run() override;

};


/*****************************************************
 *  This class derives from JEGATraits. Dakota 
 *  uses this class to verify the optimizer's 
 *  capabilities, such as support for 
 *  non-linear constraints and its types. Overrides 
 *  are generally unnecessary here as we want to
 *  mimic JEGA as much as possible. That being said,
 *  the functions are included here for clarity.
 *****************************************************/

class AdaptiveJegaTraits : public JEGATraits {

public:

  //! Default constructor
  AdaptiveJegaTraits() { };

  //! Destructor
  virtual ~AdaptiveJegaTraits() { };

  //! This is needed to handle constraints
  inline static double noValue() { return std::numeric_limits<Real>::max(); }

  //! Return the flag indicating whether method supports continuous variables
  bool supports_continuous_variables() override { return true; }

  //! Return the flag indicating whether method supports nonlinear equality constrinats
  bool supports_nonlinear_equality() override { return true; }

  //! Return the flag indicating whether method supports nonlinear inequality constrinats
  bool supports_nonlinear_inequality() override { return true; }

  //! Return format for nonlinear inequality constraints
  NONLINEAR_EQUALITY_FORMAT nonlinear_equality_format() override { 
    return NONLINEAR_EQUALITY_FORMAT::TRUE_EQUALITY; 
  }

  //! Return format for nonlinear inequality constraints
  NONLINEAR_INEQUALITY_FORMAT nonlinear_inequality_format() override { 
    return NONLINEAR_INEQUALITY_FORMAT::TWO_SIDED;
    //return NONLINEAR_INEQUALITY_FORMAT::ONE_SIDED_UPPER; 
  }

};

#endif
