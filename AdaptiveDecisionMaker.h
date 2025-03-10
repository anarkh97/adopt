#ifndef _ADAPTIVE_DECISION_MAKER_H_
#define _ADAPTIVE_DECISION_MAKER_H_

#include<array>
#include<KDTree.h>

// Dakota includes
#include<dakota_data_types.hpp>
#include<SurrogatesGaussianProcess.hpp>

class AdaptiveDecisionMaker {

  //! Database containers
  Dakota::IntRealVectorMap true_evals;
  Dakota::IntRealVectorMap approx_evals; // internal var.
  Dakota::IntRealMap       error_vals; // internal var. 

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
                           Dakota::IntVector& into, size_t num_neighbors,
                           bool force=false);
  bool GetEvaluationDecision(const Dakota::RealVector& variables);
  void RecordErrorForVariables(const Dakota::RealVector& variables, 
                               const double& error);
  bool IsEvaluationApprox(const Dakota::RealVector& variables);
  void UpdateEvaluationDecision(const int id, 
                                const Dakota::RealVector& variables,
                                const bool decision);

  //! Functions related to training and model assessment
  void Train();

};

//-----------------------------------------------------------------------------

// An instantiation of "Obj" in KDTree.h
template <int dim>
class PointInND {
public:
  int id;
  std::array<double, dim> x;
public:
  PointInND() : id(0), x({}) { };
  PointInND(int i, std::array<double,dim> &xin) : id(i), x(xin) { }
  double val(int i) const { return x[i]; }
  double width(int i) const { return 0.0; }
  int pid() const { return id; }
};

//-----------------------------------------------------------------------------
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
