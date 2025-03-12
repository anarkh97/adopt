#ifndef _SIMPLIFIED_KD_TREE_
#define _SIMPLIFIED_KD_TREE_

#include<vector>
#include<memory>
#include<cassert>
#include<algorithm>

/******************************************************************************
 * Simplified implementation of M2C's KDTree.
 * The data type for this tree is fixed to double arrays, while the size
 * of the array and dimensions are inferred from the inputs.
 *****************************************************************************/

class SimplifiedKDTree {

protected:

  int num_points; //! total number of points
  int dimensions; //! total number of dimensions
  int direction; //! current split direction

  //! points stored in contiguous array of size [num_points x dim]
  std::shared_ptr<double> stored_points;

public:

  SimplifiedKDTree();
  SimplifiedKDTree(int dim, int num_points, double* points);
  ~SimplifiedKDTree() { }; //! smart pointers are automatically deleted.
  

  //! Give the maximum width of object in this tree.
  //! note this is not the same as the width of the bounding box
  //! of the tree.
  double GetMaxWidth(int dir);

  int GetMaxDepth();
  int GetMaxLeaf();

  int UpdateKDTree(int dim, int num_points, double* new_points);

  //! Finding functions
  int FindCandidates(int dim, double* query, double* output, 
                     int max_candidates, int depth=0);
  int FindCloseCandidates(int dim, double* query, double* output,
                          int max_candidates, double distance);
  int FindCandidatesWithin(int dim, double* query, double* output,
                           int max_candidates, double distance);
  int FindCandidatesInBox(int dim, double* qmin, double* qmax, double* output,
                          int max_candidates);

protected:

  std::shared_ptr<SimplifiedKDTree> BuildKDTree(int dim, 
                                                int num_points, 
                                                double* points,
                                                int depth=0);

  int GetBestSplit(int num_points, double* objects, int dir);

};

//-----------------------------------------------------------------------------
//! Implementation
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

SimplifiedKDTree::SimplifiedKDTree()
                : dimensions(-1), num_points(-1), stored_points(nullptr),
                  direction(-1)
{
  // empty c_tor
}

//-----------------------------------------------------------------------------

SimplifiedKDTree::SimplifiedKDTree(int dim_, int n_points_, double *points_,
                                   int depth)
                : dimensions(dim_), num_points(n_points_), 
                  stored_points(points_), direction(-1)
{
  direction = depth % dimensions;
}

//-----------------------------------------------------------------------------

int
SimplifiedKDTree::BuildKDTree()
{

  assert(stored_points); // can not be null

  int first_dir = direction;
  int best_dir  = 0;
  int split;
  double best_quality = 0.0;
  std::vector<double> quality(dimensions, 0.0);

  // Any large max iterations value should be fine.
  int iter=0;
  for(; iter<500; ++iter) {
    if(num_points <= 4) {
      split = num_points;
      quality[direction] = 0.0;
      best_dir = direction;
      break;
    }

    split = GetBestSplit(direction);

    // calculate the quality of this split
    quality[direction] = std::min(split, num_points-split);
    
    // check if this is the best value seen till now.
    if(direction == first_dir || quality[direction]>best_quality) {
      best_quality = quality[direction];
      best_dir     = direction;
    }

    // check if we need to improve.
    if(quality[direction] > 0.25*num_points)
      break;
    else {
      int new_dir = (direction+1) % dimensions;
      if(new_dir == first_dir) break;
      direction = new_dir;
    }

  }

  assert(iter<500); // we should converge, if not, something went wrong.

  if(quality[best_dir] < 0.25*num_points) {
    
  }

}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

#endif
