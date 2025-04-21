#ifndef _ADAPTIVE_DECISION_MAKER_H_
#define _ADAPTIVE_DECISION_MAKER_H_

// Dakota includes
#include<dakota_global_defs.hpp>
#include<dakota_data_types.hpp>
#include<ProblemDescDB.hpp>
#include<SurrogatesGaussianProcess.hpp>

// Eigen includes.
#include<Eigen/Dense>

/************************************************
 *
 ***********************************************/

class AdaptiveDecisionMaker {

  //! Verbosity
  int verbose;

  //! flag to switch on the GP model.
  bool ready_to_predict;

  //! Number of times train is called.
  int num_train_calls;

  //! Problem description
  const Dakota::ProblemDescDB& problem_db;

  //! Database containers
  Dakota::IntRealVectorMap true_evals; // internal var.
  Dakota::IntRealVectorMap approx_evals; // internal var.
  Dakota::IntRealMap       error_vals; // internal var. 

  //! Gaussian regression model --- operates on approximation errors.
  dakota::surrogates::GaussianProcess gp_model;

  //! Neighbor search model --- (AN) currently using naive search

public:

  AdaptiveDecisionMaker(const Dakota::ProblemDescDB& problem_db);
  ~AdaptiveDecisionMaker() { };

  //! Iterating functions
  size_t GetTrueDatabaseSize() const { return true_evals.size(); }
  Dakota::IntRealVectorMap::iterator GetBeginForTrueDatabase();
  Dakota::IntRealVectorMap::const_iterator GetBeginForTrueDatabase() const;

  Dakota::IntRealVectorMap::iterator GetEndForTrueDatabase();
  Dakota::IntRealVectorMap::const_iterator GetEndForTrueDatabase() const;

  //! Functions that interface with Dakota's Iterator (Optimizer)
  //! Getters
  void GetNearestNeighbors(const Dakota::RealVector& variables,
                           Dakota::IntVector& into,
                           size_t num_int_points,
                           bool force=false);
  void GetEvaluationType(const Dakota::RealVector& variables,
                         Dakota::String& into,
                         const bool error_sim=false);

  //! Update functions
  void RecordEvaluationError(const int id, 
                             const Dakota::RealVector& variables, 
                             const double& error);
  void RecordEvaluationDecision(const int id, 
                                const Dakota::RealVector& variables,
                                const Dakota::String& eval_type);

  //! Functions related to training and model assessment
  void Train();

private:

  void ReadOptionsFile(const Dakota::String filename);

  void GetEvaluationDecision(const Dakota::RealVector& variables,
                             Dakota::String& into);

  //! Helper functions to interface with dakota::surrogates
  void LoadGaussianProcesssOptions();
  //! (re-)builds the GaussianProcess model
  bool BuildGaussianProcessModel(const Eigen::MatrixXd& samples,
                                 const Eigen::VectorXd& values,
                                 double& loss,
                                 double split_ratio=0.0,
                                 /* default value of dakota::surrogates*/
                                 const size_t seed=42); 
  void CrossValidateGausssianModel() const;

  void LoadParameters(const Dakota::IntRealVectorMap& from, 
                      Eigen::MatrixXd& into) const;
  void LoadResponses(const Dakota::IntRealMap& from, 
                     Eigen::VectorXd& into) const;

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
