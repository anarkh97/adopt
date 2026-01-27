/************************************************************************
 * MODIFICATIONS IN THIS DERIVED FILE
 *
 * Derived from:
 *   - Dakota's JEGAOptimizer::EvaluatorCreator
 *   - Upstream commit: 0daaaa2237bd79445349c47d157a4a2b73db7452
 *
 * Modifications by Aditya Narkhede, 2025:
 *   - Included a seperate Dakota::Model for error calculations.
 *   - Included a AdaptiveDecisionMaker
 *
 * Licensing:
 *   This file continues to be licensed under the GNU Lesser General Public
 *   License v2.1 or (at your option) any later version. See the bundled
 *   license file `LICENSE.LGPL-2.1` or
 *   <https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html>.
 ************************************************************************/

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

namespace MultiFidelityOptimizer
{

namespace detail
{

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
  ~MFGAEvaluatorCreator(){};

  //! Called internally by JEGA libarary. Overriden to pass Dakota Model
  JEGA::Algorithms::GeneticAlgorithmEvaluator *
  CreateEvaluator(JEGA::Algorithms::GeneticAlgorithm &alg) override;
};

} // namespace detail

} // namespace MultiFidelityOptimizer

#endif
