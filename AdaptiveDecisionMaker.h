#ifndef _ADAPTIVE_DECISION_MAKER_H_
#define _ADAPTIVE_DECISION_MAKER_H_

// Dakota includes
#include<dakota_global_defs.hpp>
#include<dakota_data_types.hpp>
#include<ProblemDescDB.hpp>
#include<SurrogatesGaussianProcess.hpp>

// Eigen includes.
#include<Eigen/Dense>

typedef std::map<int, Dakota::String> IntStringMap;

/************************************************
 *
 ***********************************************/

class AdaptiveDecisionMaker {

  int verbose; //! Verbosity
  int num_train_calls; //! Number of times model is trained.
  bool ready_to_predict; //! flag to switch on the GP model.

  //! Problem description
  const Dakota::ProblemDescDB& problem_db;

  //! Database containers
  Dakota::IntRealVectorMap id2var; // internal var.
  Dakota::IntRealMap       id2error; // internal var. 
  IntStringMap             id2type; // internal var.

  //! Gaussian regression model --- operates on approximation errors.
  dakota::surrogates::GaussianProcess gp_model;

  //! Neighbor search model --- (AN) currently using naive search

public:

  AdaptiveDecisionMaker(const Dakota::ProblemDescDB& problem_db);
  ~AdaptiveDecisionMaker() { };

  //! Iterating functions
  size_t GetTrueDatabaseSize() const { return id2var.size(); }
  Dakota::IntRealVectorMap::iterator GetBeginForTrueDatabase();
  Dakota::IntRealVectorMap::const_iterator GetBeginForTrueDatabase() const;

  Dakota::IntRealVectorMap::iterator GetEndForTrueDatabase();
  Dakota::IntRealVectorMap::const_iterator GetEndForTrueDatabase() const;

  //! Functions that interface with Dakota's Iterator (Optimizer)
  //! Getters
  void GetEvalTypeAndMetaData(const Dakota::RealVector& variables,
                              Dakota::String& into_type,
                              Dakota::IntVector& into_metadata,
                              size_t num_points,
                              bool flag=false);


  //! Update functions
  bool NeedToComputeErrors();
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

  void GetTargetAndNearestNeighbors(const Dakota::RealVector& variables,
                                    int& target,
                                    Dakota::IntVector& candidates,
                                    size_t num_int_points,
                                    bool flag=false);
  void GetEvaluationType(const Dakota::RealVector& variables,
                         Dakota::String& into,
                         bool flag=false);

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
  return id2var.begin();
}

inline Dakota::IntRealVectorMap::const_iterator
AdaptiveDecisionMaker::GetBeginForTrueDatabase() const
{
  return id2var.begin();
}

inline Dakota::IntRealVectorMap::iterator 
AdaptiveDecisionMaker::GetEndForTrueDatabase()
{
  return id2var.end();
}

inline Dakota::IntRealVectorMap::const_iterator
AdaptiveDecisionMaker::GetEndForTrueDatabase() const
{
  return id2var.end();
}

#endif
