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


#ifndef _ADAPTIVE_JEGA_OPTIMIZER_H_
#define _ADAPTIVE_JEGA_OPTIMIZER_H_

//! Optimization related headers
#include<JEGAOptimizer.hpp>
#include<JegaEvaluatorCreator.h>
#include<JegaEvaluator.h>
#include<JegaDriver.h>
//! Surrogate model related headers
#include<AdaptiveDecisionMaker.h>

/*****************************************************
 *  This class extends JEGAOptimizer, Dakota's 
 *  implementation of John Eddy's Genetic Algorithm 
 *  (JEGA) library. We modify its optimization loop to 
 *  integrate an Adaptive Decision Model based on 
 *  Gaussian Process. This model predicts whether 
 *  high-fidelity fluid dynamics simulations (M2C) are 
 *  necessary or if simpler interpolation estimates 
 *  (FEST) will suffice. Note that Dakota::JEGAOptimizer's 
 *  Evaluator, Driver, and EvaluatorCreator classes 
 *  along with ParameterDatabase variables are private 
 *  and inaccessible here. Hence, we do not directly 
 *  derive from JEGAOptimizer but refactor most of its 
 *  functionality.
 *****************************************************/

class AdaptiveJegaOptimizer : public Dakota::Optimizer {

  std::shared_ptr<JEGA::Utilities::ParameterDatabase> param_db; 
  std::shared_ptr<JegaEvaluatorCreator> eval_creator;
  std::shared_ptr<AdaptiveDecisionMaker> decision_maker;

public:

  //! Default constructor
  AdaptiveJegaOptimizer(Dakota::ProblemDescDB& prob_db, 
                        std::shared_ptr<Dakota::Model> sim_model,
                        std::shared_ptr<Dakota::Model> err_model);

  //! Default Destructor -- Calls ~Optimizer automatically
  ~AdaptiveJegaOptimizer() override;

  //! Interface to Dakota's core_run
  void core_run() override;

  //! Note we can return multiple best solutions even in the case of SOGA.
  bool returns_multiple_points() const override { return true; };

private:

  //! JEGA related methods
  void LoadParameterDatabase();
  void LoadAlgorithmConfig(JEGA::FrontEnd::AlgorithmConfig& config);
  void LoadProblemConfig(JEGA::FrontEnd::ProblemConfig& config);
  void LoadDesignVariables(JEGA::FrontEnd::ProblemConfig& config);
  void LoadObjectiveFunctions(JEGA::FrontEnd::ProblemConfig& config);
  void LoadConstraints(JEGA::FrontEnd::ProblemConfig& config);

  //! Dakota related methods
  //! MOGA functionality has been removed for now.
  void GetBestSOSolutions(const JEGA::Utilities::DesignOFSortSet& from,
                          const JEGA::Algorithms::GeneticAlgorithm& ga,
                          std::multimap<Dakota::RealRealPair, 
                                        JEGA::Utilities::Design*>&);

  void LoadDakotaResponses(const JEGA::Utilities::Design& from, 
                           Dakota::Variables& vars, 
                           Dakota::Response& resp) const;

};


/*****************************************************
 *  This class derives from TraitsBase. Dakota 
 *  uses this class to verify the optimizer's 
 *  capabilities, such as support for 
 *  non-linear constraints and its types. Hre as we 
 *  want to mimic JEGATraits as much as possible. 
 *****************************************************/

class AdaptiveJegaTraits : public Dakota::TraitsBase {

public:

  //! Default constructor
  AdaptiveJegaTraits() { };

  //! Destructor
  virtual ~AdaptiveJegaTraits() { };

  //! This is needed to handle constraints
  inline static double noValue() { return std::numeric_limits<Dakota::Real>::max(); }

  //! Return the flag indicating whether method supports continuous variables
  bool supports_continuous_variables() override { return true; }

  //! Return the flag indicating whether method supports nonlinear equality constrinats
  bool supports_nonlinear_equality() override { return true; }

  //! Return the flag indicating whether method supports nonlinear inequality constrinats
  bool supports_nonlinear_inequality() override { return true; }

  //! Return format for nonlinear inequality constraints
  Dakota::NONLINEAR_EQUALITY_FORMAT nonlinear_equality_format() override { 
    return Dakota::NONLINEAR_EQUALITY_FORMAT::TRUE_EQUALITY; 
  }

  //! Return format for nonlinear inequality constraints
  Dakota::NONLINEAR_INEQUALITY_FORMAT nonlinear_inequality_format() override { 
    return Dakota::NONLINEAR_INEQUALITY_FORMAT::TWO_SIDED;
    //return NONLINEAR_INEQUALITY_FORMAT::ONE_SIDED_UPPER; 
  }

};

#endif
