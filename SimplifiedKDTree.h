#ifndef _SIMPLIFIED_KD_TREE_
#define _SIMPLIFIED_KD_TREE_

#include<vector>
#include<memory>
#include<cfloat>
#include<string>
#include<cassert>
#include<iterator>
#include<iostream>
#include<algorithm>
#include<type_traits>

//-----------------------------------------------------------------------------
// Helpers
//-----------------------------------------------------------------------------

template<typename CoordinateValueType>
class DefaultComparator;

// deduce things at compile time
template <typename CoordinateType>
concept IsBidrectional = requires(CoordinateType v) {
  v.second;
};

template<typename CoordinateType>
constexpr decltype(auto)
GetCoordinate(const CoordinateType& value)
{
  if constexpr (IsBidrectional<CoordinateType>)
    return value.second; // bidirectional iterators
  else
    return value; // for random access iterators
};

/******************************************************************************
 * Simplified implementation of M2C's KDTree, it now takes the dimension as
 * inputs. Additionally the code uses iterators instead of raw pointers.
 *****************************************************************************/


template<
typename CoordinateType, 
typename Compare,
typename Iterator
>
class SimplifiedKDTree {

protected:

  int num_points; //! total number of points
  int num_dims; //! total number of dimensions
  int split_dir; //! current split direction
  double split_value; //! value at which tree was split
  bool has_child; //! flag for sub-trees

  //! points at this (sub-)tree.
  Iterator points_b;
  Iterator points_e;

  //! bounds for current (sub-) tree
  std::vector<CoordinateType> low;
  std::vector<CoordinateType> upp;

  //! sub trees
  std::shared_ptr<SimplifiedKDTree<CoordinateType,Compare,Iterator>> 
  left, right;

public:

  SimplifiedKDTree() = default;
  SimplifiedKDTree(int k,
                   const Iterator& input_b,
                   const Iterator& input_e,
                   int depth=0);
  ~SimplifiedKDTree() { }; //! smart pointers are automatically deleted.
  

/*
  //! Give the maximum width of object in this tree.
  //! note this is not the same as the width of the bounding box
  //! of the tree.
  double GetMaxWidth(int dir);
*/

  int GetMaxDepth() {
    if(!has_child) 
      return 1;
    else return std::max(left->GetMaxDepth(), right->GetMaxDepth());
  };
  int GetMaxLeaf() {
    if(!has_child) 
      return num_points;
    else return std::max(left->GetMaxLeaf(), right->GetMaxLeaf());
  };

  //! Finding functions
  int FindCandidates(const CoordinateType& query, 
                     std::vector<CoordinateType>& output,
                     int max_candidates,
                     int depth=0);
  int FindCloseCandidates(const CoordinateType& query, 
                          std::vector<CoordinateType>& output,
                          int max_candidates,
                          double distance);
  int FindCandidatesWithin(const CoordinateType& query, 
                           std::vector<CoordinateType>& output,
                           int max_candidates,
                           double distance);
  int FindCandidatesInBox(const CoordinateType& qmin,/* query bounding-box min */
                          const CoordinateType& qmax,/* query bounding-box max*/
                          std::vector<CoordinateType>& output,
                          int max_candidates);

  void Print(std::ostream& os, int depth=0);

private:

/*
  std::shared_ptr<SimplifiedKDTree> BuildKDTree(int k, 
                                                int num_points, 
                                                double* ipoints,
                                                int depth=0);
*/

  int GetBestSplit(const Iterator input_b,
                   const Iterator input_e,
                   int dir);
  void SetBoundingBox();


};

//-----------------------------------------------------------------------------
//! Implementation
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

/*
template<typename CoordinateType, typename Compare, typename Iterator>
SimplifiedKDTree<CoordinateType,Compare, Iterator>::
SimplifiedKDTree() : num_dims(-1), num_points(-1), 
                     points_b(), points_e(), 
                     split_dir(-1), split_value(-DBL_MAX),
                     has_child(false)
{
  // empty c_tor
}
*/

//-----------------------------------------------------------------------------

template<typename CoordinateType,typename Compare, typename Iterator>
SimplifiedKDTree<CoordinateType,Compare,Iterator>::
SimplifiedKDTree(int k,
                 const Iterator& input_b,
                 const Iterator& input_e,
                 int depth/*=0*/) : points_b(), points_e(), 
                                    left(nullptr), right(nullptr)
{

  // empty tree
  if(input_b == input_e) return;

  num_dims   = k;
  num_points = std::distance(input_b, input_e);
  split_dir  = depth % k;

  assert(num_points>0);

  int first_dir = split_dir;
  int best_dir  = 0;
  int split;
  double best_quality = 0.0;
  std::vector<double> quality(num_dims, 0.0);

  // Any large max iterations value should be fine.
  int iter=0;
  for(; iter<500; ++iter) {
    if(num_points <= 4) {
      split = num_points;
      quality[split_dir] = 0.0;
      best_dir = split_dir;
      break;
    }

    split = GetBestSplit(input_b, input_e, split_dir);

    // calculate the quality of this split
    quality[split_dir] = std::min(split, num_points-split);
    
    // check if this is the best value seen till now.
    if(split_dir == first_dir || quality[split_dir]>best_quality) {
      best_quality = quality[split_dir];
      best_dir     = split_dir;
    }

    // check if we need to improve.
    if(quality[split_dir] > 0.25*num_points)
      break;
    else {
      int new_dir = (split_dir+1) % num_dims;
      if(new_dir == first_dir) break;
      split_dir = new_dir;
    }

  }

  assert(iter<500); // we should converge, if not, something went wrong.

  // base case --- store and return
  if(quality[best_dir] < 0.25*num_points) {
    // make copies of begin and end iterators.
    has_child = false;
    points_b  = input_b;
    points_e  = input_e;
    // set bounding box
    SetBoundingBox();
    return;
  }

  // we have subtrees
  if(split_dir != best_dir) {
    split     = GetBestSplit(input_b, input_e, best_dir);
    split_dir = best_dir;
  }

  const Iterator &sp  =
    std::next(input_b, split);
  const Iterator &spm = 
    std::next(input_b, split-1);

  split_value = (GetCoordinate(*sp)[split_dir] + GetCoordinate(*spm)[split_dir])/2.0;

  has_child = true;
  left = std::make_shared<SimplifiedKDTree<CoordinateType,Compare,Iterator>>(
    num_dims, 
    input_b, 
    sp,
    split_dir+1);
  right = std::make_shared<SimplifiedKDTree<CoordinateType,Compare,Iterator>>(
    num_dims, 
    std::next(sp),
    input_e,
    split_dir+1);

}

//-----------------------------------------------------------------------------

/*
std::shared_ptr<SimplifiedKDTree>
SimplifiedKDTree::BuildKDTree(int k, int input_size, double *input, int depth)
{

  direction = depth % k;
  int first_dir = dir;
  int best_dir  = 0;
  int split;
  double best_quality = 0.0;
  std::vector<double> quality(k, 0.0);

  // Any large max iterations value should be fine.
  int iter=0;
  for(; iter<500; ++iter) {
    if(input_size <= 4) {
      split = input_size;
      quality[dir] = 0.0;
      best_dir = dir;
      break;
    }

    split = GetBestSplit(dir);

    // calculate the quality of this split
    quality[dir] = std::min(split, input_size-split);
    
    // check if this is the best value seen till now.
    if(dir == first_dir || quality[dir]>best_quality) {
      best_quality = quality[dir];
      best_dir     = dir;
    }

    // check if we need to improve.
    if(quality[dir] > 0.25*input_size)
      break;
    else {
      int new_dir = (dir+1) % k;
      if(new_dir == first_dir) break;
      dir = new_dir;
    }

  }

  assert(iter<500); // we should converge, if not, something went wrong.

  // base case --- update current instance and return
  if(quality[best_dir] < 0.25*input_size) {
    points      = std::make_shared<double>(input);
    num_dims      = k;
    num_points  = input_size;
    split_value = (input[split*dim+dir] + input[(split-1)*dim+dir])/2.0;
    return std::make_shared<SimplifiedKDTree>(this);
  }

  // we have subtrees
  if(dir != best_dir) {
    split = GetBestSplit(input_size, input, best_dir);
    dir   = best_dir;
  }

  split_value = (input[split*dim+dir] + input[(split-1)*dim+dir])/2.0;
  points = nullptr;
  left = BuildKDTree(k, split, input, dir+1);
  right = BuildKDTree(k, input_size-split, input, dir+1);

}

*/

//-----------------------------------------------------------------------------

template<typename CoordinateType, typename Compare, typename Iterator>
int
SimplifiedKDTree<CoordinateType,Compare,Iterator>::
GetBestSplit(const Iterator input_b,
             const Iterator input_e,
             int dir)
{

  assert(num_dims > 0);
  assert(input_b != input_e);

  int input_size = std::distance(input_b, input_e);

  assert(input_size>0);

  if(input_size<=4) return input_size;

  // input could be an bidirectional contianer type (e.g., map, set, etc.)
  using CoordinateValueType = decltype(GetCoordinate(*input_b));
  std::vector<CoordinateValueType> points;

  for(auto it=input_b; it != input_e; ++it)
    points.push_back(GetCoordinate(*it));

  // sort based on co-ordinate.
  Compare compare(dir);
  std::sort(points.begin(), points.end(), compare);

  int left_split  = input_size/2;
  int right_split = left_split;

  // Find the best median. Make sure that the split points
  // is such that the object at the median and the one before it
  // are not collocated.
  while(left_split > 0) {

/*
    const Iterator &lsp  = 
      std::next(input_b, left_split);
    const Iterator &lspm = 
      std::next(input_b, left_split-1);
    if(GetCoordinate(*lsp)[dir] == GetCoordinate(*lspm)[dir])
      left_split--;
    else
      break;
*/

    const CoordinateValueType &lsp  = points[left_split];
    const CoordinateValueType &lspm = points[left_split-1];
    if(lsp[dir] == lspm[dir])
      left_split--;
    else
      break;

  }

  // Note that since input_size>4 we are certain that right_split
  // can not be 0.
  while(right_split < input_size) {

/*
    const Iterator &rsp  = 
      std::next(input_b, right_split);
    const Iterator &rspm = 
      std::next(input_b, right_split-1);

    if(GetCoordinate(*rsp)[dir] == GetCoordinate(*rspm)[dir])
      right_split++;
    else
      break;
*/

    const CoordinateValueType &rsp  = points[right_split];
    const CoordinateValueType &rspm = points[right_split-1];
    if(rsp[dir] == rspm[dir])
      right_split++;
    else
      break;

  }

  int split = 
    (input_size-left_split < right_split) ? left_split : right_split;
  return split;

}

//-----------------------------------------------------------------------------

template<typename CoordinateType, typename Compare, typename Iterator>
int
SimplifiedKDTree<CoordinateType,Compare,Iterator>::
FindCandidates(const CoordinateType& query, 
               std::vector<CoordinateType>& output,
               int max_candidates,
               int depth)
{
/*
  // base case
  if(!has_child) {

    int num_candidates=0;
    Iterator pit = points_b; // points stored
    for(; pit!=points_e; ++pit) {
      for(int k=0; k<num_dims; ++k) {
        double val = std::abs(GetCoordinate(query)[k]-GetCoordinate(*pit)[k]);
        local_dist = std::max(local_dist, val);
      }

      if(local_dist > distance)
        continue;
      if(num_candidates >= max_candidates)
        return max_candidates+1; // found maximum requested candidates.

      // update iterators
      output.push_back(*pit);
      num_candidates++;
    }

    assert(num_candidates == output.size());

    return num_candidates;

  }

  int num_found = 0;

  // go to left tree if current co-ordinate value is less than
  // split value.
  if(query[split_dir] <= split_value and left)
    num_found = left->FindCandidates(query, output, max_candidates, depth+1);

  // go to the right tree if the current co-ordinate value is 
  // less than split value.
  if(query[split_dir] >= split_value)
    num_found += right->FindCandidates(query, output, max_candidates, depth+1);

  return num_found;
*/
}

//-----------------------------------------------------------------------------

template<typename CoordinateType, typename Compare, typename Iterator>
int
SimplifiedKDTree<CoordinateType,Compare,Iterator>::
FindCloseCandidates(const CoordinateType& query, 
                    std::vector<CoordinateType>& output,
                    int max_candidates,
                    double distance)
{

}

//-----------------------------------------------------------------------------

template<typename CoordinateType, typename Compare, typename Iterator>
int
SimplifiedKDTree<CoordinateType,Compare,Iterator>::
FindCandidatesWithin(const CoordinateType& query, 
                     std::vector<CoordinateType>& output,
                     int max_candidates,
                     double distance)
{

  // base case
  if(!has_child) {
    
    int num_candidates=0;
    Iterator pit = points_b; // stored points
    for(; pit!=points_e; ++pit) {
      double local_dist = 0;
      for(int k=0; k<num_dims; ++k) {
        double val = std::abs(GetCoordinate(query)[k]-GetCoordinate(*pit)[k]);
        local_dist = std::max(local_dist, val);
      }

      if(local_dist > distance)
        continue;
      if(num_candidates >= max_candidates)
        return max_candidates+1; // found maximum requested candidates.

      // update output
      output.push_back(*pit);
      num_candidates++;
    }

    assert(num_candidates == output.size());

    return num_candidates;

  }

  int num_found = 0;

  // go to left tree if current co-ordinate value is less than
  // split value.
  if(GetCoordinate(query)[split_dir] <= split_value and left)
    num_found = left->FindCandidatesWithin(query, output, max_candidates,
                                           distance);

  // go to the right tree if the current co-ordinate value is 
  // less than split value.
  if(GetCoordinate(query)[split_dir] >= split_value and right)
    num_found += right->FindCandidatesWithin(query, output, max_candidates,
                                             distance);

  return num_found;

}


//-----------------------------------------------------------------------------

template<typename CoordinateType, typename Compare, typename Iterator>
int
SimplifiedKDTree<CoordinateType,Compare,Iterator>::
FindCandidatesInBox(const CoordinateType& qmin,/* query bounding-box min */
                    const CoordinateType& qmax,/* query bounding-box max*/
                    std::vector<CoordinateType>& output,
                    int max_candidates)
{

}


//-----------------------------------------------------------------------------

template<typename CoordinateType, typename Compare, typename Iterator>
void
SimplifiedKDTree<CoordinateType,Compare,Iterator>::
SetBoundingBox()
{

}

//-----------------------------------------------------------------------------

template<typename CoordinateType, typename Compare, typename Iterator>
void
SimplifiedKDTree<CoordinateType,Compare,Iterator>::
Print(std::ostream &os, int depth)
{

  std::string indent(depth, '  ');

  // current tree stats.
  os << indent << "Node:\n"
     << indent << " num_points:" << num_points << "\n"
     << indent << " split_dir:" << split_dir << "\n"
     << indent << " split_value:" << split_value << "\n"
     << indent << " sub-trees?:" << (has_child ? "yes" : "no") << "\n";

  // base case
  if(!has_child) {
    os << indent << "  Leaf with " << num_points << " point(s).\n";
    return;
  }

  if(left) {
    os << indent << "Left tree (depth: " << depth << "):\n";
    left->Print(os, depth+1);
  }

  if(right) {
    os << indent << "Right tree (depth: " << depth <<  "):\n";
    right->Print(os, depth+1);
  }

}

//-----------------------------------------------------------------------------

template<typename CoordinateValueType>
class DefaultComparator {
  int dir;
public:
  DefaultComparator(int d) : dir(d) { }
  bool operator()(const CoordinateValueType &lhs, 
                  const CoordinateValueType &rhs) const {
    return lhs[dir] < rhs[dir];
  }
};

//-----------------------------------------------------------------------------

// CTAD (Class Template Argument Deduction) for compiler refer 
// 1. https://en.cppreference.com/w/cpp/language/class_template_argument_deduction
// 2. https://stackoverflow.com/questions/40951697/what-are-template-deduction-guides-and-when-should-we-use-them

// Deduce everything from the iterator type.
template<typename Iterator> 
SimplifiedKDTree(int, const Iterator&, const Iterator&, int) -> 
  SimplifiedKDTree<
     typename std::iterator_traits<Iterator>::value_type,
     /*DefaultComparator<typename std::iterator_traits<Iterator>::value_type>,*/
     DefaultComparator<
      std::decay_t<decltype(GetCoordinate(*std::declval<Iterator>()))> 
     >, 
     Iterator
    >;

template<typename Iterator> 
SimplifiedKDTree(int, const Iterator&, const Iterator&) -> 
  SimplifiedKDTree<
     typename std::iterator_traits<Iterator>::value_type,
     /*DefaultComparator<typename std::iterator_traits<Iterator>::value_type>,*/
     DefaultComparator<
      std::decay_t<decltype(GetCoordinate(*std::declval<Iterator>()))> 
     >, 
     Iterator
    >;

//-----------------------------------------------------------------------------

#endif
