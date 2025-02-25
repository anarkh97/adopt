#ifndef _ADAPTIVE_DECISION_MAKER_H_
#define _ADAPTIVE_DECISION_MAKER_H_

#include<KDTree.h>

// Dakota includes
#include<dakota_data_types.hpp>
#include<SurrogatesGaussianProcess.hpp>

class AdaptiveDecisionMaker {

  dakota::surrogates::GaussianProcess gp_model;

public:

  AdaptiveDecisionMaker() { };
  ~AdaptiveDecisionMaker() { };

  //! Functions that interface with Dakota's Iterator
  void GetNearestNeighbors(Dakota::IntVector& into, size_t num_neighbors);
  void GetAndUpdateEvaluationDecision(const Dakota::RealVector& variables, 
                                      size_t num_vars, bool& decision);
  void RecordErrorForVariables(const Dakota::RealVector& variables, 
                               size_t num_vars, 
                               const double& error);
  bool IsEvaluationApprox(const Dakota::RealVector& variables);
  bool AddToDatabase(const Dakota::RealVector& variables);

  //! Functions related to training and model assessment
  void Train();

};

#endif
