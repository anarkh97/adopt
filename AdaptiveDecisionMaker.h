#ifndef _ADAPTIVE_DECISION_MAKER_H_
#define _ADAPTIVE_DECISION_MAKER_H_

#include<KDTree.h>

// Dakota includes
#include<dakota_data_types.hpp>
#include<SurrogatesGaussianProcess.hpp>

class AdaptiveDecisionMaker {

  //! Database containers
  Dakota::IntRealVectorMap true_evals;
  Dakota::IntRealVectorMap approx_evals; // internal var.

  //! Gaussian regression model --- operates on approximation errors.
  dakota::surrogates::GaussianProcess gp_model;

public:

  AdaptiveDecisionMaker() { };
  ~AdaptiveDecisionMaker() { };

  //! Iterating functions
  Dakota::IntRealVectorMap::iterator GetBeginForTrueDatabase();
  Dakota::IntRealVectorMap::const_iterator GetBeginForTrueDatabase() const;

  Dakota::IntRealVectorMap::iterator GetEndForTrueDatabase();
  Dakota::IntRealVectorMap::const_iterator GetEndForTrueDatabase() const;

  //! Functions that interface with Dakota's Iterator (Optimizer)
  void GetNearestNeighbors(const Dakota::RealVector& variables,
                           Dakota::IntVector& into, size_t num_neighbors);
  void GetEvaluationDecision(const Dakota::RealVector& variables, 
                             bool& decision);
  void RecordErrorForVariables(const Dakota::RealVector& variables, 
                               const double& error);
  bool IsEvaluationApprox(const Dakota::RealVector& variables);
  bool UpdateEvaluationDecision(int id, const Dakota::RealVector& variables,
                                const Dakota::String& decision);

  //! Functions related to training and model assessment
  void Train();

};

// Inline definitions

inline Dakota::IntRealVectorMap::iterator 
AdaptiveDecisionMaker::GetBeginForTrueDatabase()
{
  return true_evals.begin();
}

inline Dakota::IntRealVectorMap::const_iterator
AdaptiveDecisionMaker::GetBeginForTrueDatabase() const
{
  return true_evals.begin();
}

inline Dakota::IntRealVectorMap::iterator 
AdaptiveDecisionMaker::GetEndForTrueDatabase()
{
  return true_evals.end();
}

inline Dakota::IntRealVectorMap::const_iterator
AdaptiveDecisionMaker::GetEndForTrueDatabase() const
{
  return true_evals.end();
}

#endif
