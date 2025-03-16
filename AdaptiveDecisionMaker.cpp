#include<algorithm>
#include<AdaptiveDecisionMaker.h>

// Dakota header for read/write
// Note: This header also has functions for 
// read/write from/to file.
#include<dakota_data_io.hpp>

using namespace std;
using namespace Dakota;
using namespace dakota::surrogates;

//------------------------------------------------------------------------------
//! Helper function
double EuclideanDistance(const RealVector &v1, const RealVector &v2)
{
  assert(v1.length() == v2.length());
  double dist=0;
  int dim=v1.length();

  for(int i=0; i<dim; ++i)
    dist += pow(v1[i]-v2[i],2);

  return sqrt(dist);

}

//------------------------------------------------------------------------------

AdaptiveDecisionMaker::AdaptiveDecisionMaker()
{

  // Initialize GP model

}

//------------------------------------------------------------------------------

void AdaptiveDecisionMaker::GetNearestNeighbors(const RealVector &cont_vars, 
                                                IntVector &into, 
                                                size_t num_neighbors,
                                                bool force)
{

  if(into.length() != num_neighbors) 
    into.size(num_neighbors);

  // AN: Throw an error if force find is true but true-db is empty.
  // This distinction is made so that when decision maker is called
  // for error calculations one can check if database with true 
  // evaluations is empty.
  if(true_evals.empty() and force) {
    Cerr << "ADAPTIVE SOGA: Trying to find nearest neighbors "
         << "in an empty database.\n";
    abort_handler(METHOD_ERROR);
  }

  // Send garbage values when true-db is empty.
  if(true_evals.empty()) {
    for(int i=0; i<num_neighbors; ++i) into[i] = -1;
    return;
  }

  if(true_evals.size() < num_neighbors) {
    Cerr << "ADAPTIVE SOGA: Number of true evaluations stored "
         << "(" << true_evals.size() << ") "
         << "is less than number of neighbors requested "
         << "(" << num_neighbors << ").\n";
    abort_handler(METHOD_ERROR);
  }

  // AN: Currently doing naive search.
  // TODO: Could use KDTree.
  
  vector<pair<double, int>> dist2targ;
  IntRealVectorMap::const_iterator tr_it = true_evals.begin();
  const IntRealVectorMap::const_iterator tr_e = true_evals.end();
  for(; tr_it!=tr_e; ++tr_it) {
    double d = EuclideanDistance(tr_it->second, cont_vars);
    dist2targ.push_back({d, tr_it->first});
  }

  sort(dist2targ.begin(), dist2targ.end());

  // fill up first "num_neighbors" into IntVector
  for(int i=0; i<num_neighbors; ++i) {
    into[i] = dist2targ[i].second;
  }

}

//------------------------------------------------------------------------------

bool
AdaptiveDecisionMaker::GetEvaluationDecision(const RealVector &cont_vars)
{
  
  // Send "TRUE" decision when true-db is empty.
  if(true_evals.empty()) return true;
  return false;
}

//------------------------------------------------------------------------------

void
AdaptiveDecisionMaker::RecordEvaluationError(const int eval_id,
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

void
AdaptiveDecisionMaker::RecordEvaluationDecision(int eval_id, 
                                                const RealVector &cont_vars, 
                                                const bool eval_type)
{

  if(eval_type == true) {
    // First check if evaluation id has already been mapped or not.
    IntRealVectorMap::iterator it = true_evals.find(eval_id);
    if(it != true_evals.end()) {
      Cerr << "ADAPTIVE SOGA: Duplicate evaluation ID " << eval_id
           << ", previously mapped to variables: " << it->second;
      abort_handler(METHOD_ERROR);
    }

    true_evals[eval_id] = cont_vars;
  }
  else {
    // Firt check if evaluation id has already been mapped or not.
    IntRealVectorMap::iterator it = approx_evals.find(eval_id);
    if(it != approx_evals.end()) {
      Cerr << "ADAPTIVE SOGA: Duplicate evaluation ID " << eval_id
           << ", previously mapped to variables: " << it->second;
      abort_handler(METHOD_ERROR);
    }

    approx_evals[eval_id] = cont_vars;
  }

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

