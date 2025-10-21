/************************************************************************
 * MODIFICATIONS IN THIS DERIVED FILE
 *
 * Derived from:
 *   - Dakota's JEGAOptimizer::Driver
 *   - Upstream commit: 0daaaa2237bd79445349c47d157a4a2b73db7452
 *
 * Modifications by Aditya Narkhede, 2025:
 *
 * Licensing:
 *   This file continues to be licensed under the GNU Lesser General Public
 *   License v2.1 or (at your option) any later version. See the bundled
 *   license file `LICENSE.LGPL-2.1` or
 *   <https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html>.
 ************************************************************************/

#ifndef _JEGA_DRIVER_H_
#define _JEGA_DRIVER_H_

#include<../Utilities/include/JEGAConfig.hpp>
#include<../Utilities/include/Logging.hpp>

// JEGA core includes.
#include<GeneticAlgorithm.hpp>
#include<GeneticAlgorithmEvaluator.hpp>

// JEGA front end includes.
#include<../FrontEnd/Core/include/Driver.hpp>
#include<../FrontEnd/Core/include/ProblemConfig.hpp>
#include<../FrontEnd/Core/include/AlgorithmConfig.hpp>
#include<../FrontEnd/Core/include/EvaluatorCreator.hpp>

// Eddy utility includes.
#include<utilities/include/EDDY_DebugScope.hpp>

// JEGA utility includes.
#include<../Utilities/include/DesignGroup.hpp>

// Include the main optimizer class to get Model and other dakota things.
#include<DakotaOptimizer.hpp>

/************************************************
 * JegaDriver is the core GA optimization driver 
 * class, managing the execution of optimization 
 * iterations. Here we implement a simple wrapper
 * for JEGA::FrontEnd::Driver and add our 
 * adaptive decision maker.
 ***********************************************/

class JegaDriver : public JEGA::FrontEnd::Driver {

public:

  //! Constructor
  JegaDriver(const JEGA::FrontEnd::ProblemConfig &config);

  //! Wrappers
  JEGA::Algorithms::GeneticAlgorithm* 
  ExtractAllData(const JEGA::FrontEnd::AlgorithmConfig &config) {
    return JEGA::FrontEnd::Driver::ExtractAllData(config);
  }

  void DestroyAlgorithm(JEGA::Algorithms::GeneticAlgorithm *ga) {
    JEGA::FrontEnd::Driver::DestroyAlgorithm(ga);
  }

  //! Custom implementations of JEGA::FrontEnd::Driver functions
  JEGA::Utilities::DesignOFSortSet 
  PerformIterations(JEGA::Algorithms::GeneticAlgorithm* ga);

};

#endif
