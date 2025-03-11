#include<AdaptiveDecisionMaker.h>

using namespace std;
using namespace Dakota;
using namespace dakota::surrogates;

//------------------------------------------------------------------------------

void AdaptiveDecisionMaker::GetNearestNeighbors(const RealVector &cont_vars, 
                                                IntVector &into, 
                                                size_t num_neighbors,
                                                bool force)
{

  // AN: Throw an error if force find is true but true-db is empty.
  // This distinction is made so that when decision maker is called
  // for error calculations one can check if database with true 
  // evaluations is empty.
  if(true_evals.empty() and force) {
    Cerr << "ADAPTIVE SOGA: Trying to find nearest neighbors in an empty "
         << "database.\n";
    abort_handler(METHOD_ERROR);
  }

  // Send garbage values when true-db is empty.
  if(true_evals.empty()) {
    if(into.length() != num_neighbors) into.size(num_neighbors);
    for(int i=0; i<num_neighbors; ++i) into[i] = -1;
  }

}

//------------------------------------------------------------------------------

bool
AdaptiveDecisionMaker::GetEvaluationDecision(const RealVector &cont_vars)
{
  
  // Send "TRUE" decision when true-db is empty.
  if(true_evals.empty()) return true;

}	

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::RecordErrorForVariables(const int eval_id,
                                               const RealVector &cont_vars,
                                               const double &error)
{

  const RealVector &stored_vars = true_evals[eval_id];

  // Check if variables match.  This is done to ensure error values stored
  // in error_vals are mapped to expected continuous variables stored in 
  // true_evals.
  if(cont_vars != stored_vars) {
    Cerr << "ADAPTIVE SOGA: Trying to map error values for unknown design "
         << "variables.\n";
    abort_handler(METHOD_ERROR);
  }

  error_vals[eval_id] = error;

}

//------------------------------------------------------------------------------

bool
AdaptiveDecisionMaker::IsEvaluationApprox(const RealVector &cont_vars)
{

}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::UpdateEvaluationDecision(int id, 
                                                const RealVector &cont_vars, 
                                                const bool decision)
{

}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::Train()
{

}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

