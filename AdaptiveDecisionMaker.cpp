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

  // Send garbage values when true-db is empty.
  // AN: Throw an error if force find is true but true-db is empty.

}

//------------------------------------------------------------------------------

bool
AdaptiveDecisionMaker::GetEvaluationDecision(const RealVector &cont_vars)
{
  
  // Send "TRUE" decision when true-db is empty.

}	

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::RecordErrorForVariables(const RealVector &cont_vars,
                                               const double &error)
{

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

