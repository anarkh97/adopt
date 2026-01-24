#ifndef _SIMPLIFIED_KD_TREE_
#define _SIMPLIFIED_KD_TREE_

#include <cmath>
#include <queue>
#include <vector>
#include <memory>
#include <limits>
//#include<string>
#include <cassert>
#include <iterator>
//#include<iostream>
#include <algorithm>
#include <type_traits>

//-----------------------------------------------------------------------------
// Helpers
//-----------------------------------------------------------------------------

template <typename CoordinateContainerType>
double EuclideanDistance(int dim, const CoordinateContainerType &p1,
                         const CoordinateContainerType &p2)
{
  double dist;
  for (int i = 0; i < dim; ++i)
    dist += std::pow(p1[i] - p2[i], 2);
  return std::sqrt(dist);
}

template <typename CoordinateContainerType> class DefaultComparator;

// deduce things at compile time
template <typename CoordinateType>
concept IsBidrectional = requires(CoordinateType v) { v.second; };

template <typename CoordinateType>
constexpr decltype(auto) GetCoordinate(const CoordinateType &value)
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

template <typename CoordinateType, typename Compare, typename Iterator>
class SimplifiedKDTree
{

private:
  //! Used for storing candidates based on distance.
  using CandidateDistPair      = std::pair<double, const CoordinateType *>;
  using CandidatePriorityQueue = std::priority_queue<CandidateDistPair>;

  //! Coordinate container deduction at compile time.
  using CoordinateContainerType
    = decltype(GetCoordinate(*std::declval<Iterator>()));
  using CoordinateElementType = typename CoordinateContainerType::value_type;

protected:
  int    num_points;  //! total number of points
  int    num_dims;    //! total number of dimensions
  int    split_dir;   //! current split direction
  double split_value; //! value at which tree was split
  bool   has_child;   //! flag for sub-trees

  //! points at this (sub-)tree.
  Iterator points_b;
  Iterator points_e;

  //! bounds for current (sub-) tree
  std::vector<CoordinateElementType> low;
  std::vector<CoordinateElementType> upp;

  //! sub trees
  std::shared_ptr<SimplifiedKDTree<CoordinateType, Compare, Iterator>> left,
    right;

public:
  SimplifiedKDTree() = default;
  SimplifiedKDTree(int k, const Iterator &input_b, const Iterator &input_e,
                   int depth = 0);
  ~SimplifiedKDTree(){}; //! smart pointers are automatically deleted.

  /*
  //! Give the maximum width of object in this tree.
  //! note this is not the same as the width of the bounding box
  //! of the tree.
  double GetMaxWidth(int dir);
*/

  int GetMaxDepth()
  {
    if (!has_child)
      return 1;
    else
      return std::max(left->GetMaxDepth(), right->GetMaxDepth());
  };
  int GetMaxLeaf()
  {
    if (!has_child)
      return num_points;
    else
      return std::max(left->GetMaxLeaf(), right->GetMaxLeaf());
  };

  //! Finding functions
  int FindCandidates(const CoordinateType        &query,
                     std::vector<CoordinateType> &output,
                     int                          max_candidates) const;
  int FindCandidatesWithin(const CoordinateType        &query,
                           std::vector<CoordinateType> &output,
                           int max_candidates, double distance) const;
  int FindCandidatesInBox(
    const CoordinateType        &qmin, /* query bounding-box min */
    const CoordinateType        &qmax, /* query bounding-box max*/
    std::vector<CoordinateType> &output, int max_candidates) const;

  //void Print(std::ostream& os, int depth=0) const;

private:
  /*
  std::shared_ptr<SimplifiedKDTree> BuildKDTree(int k, 
                                                int num_points, 
                                                double* ipoints,
                                                int depth=0);
*/

  int  GetBestSplit(const Iterator &input_b, const Iterator &input_e, int dir);
  void SetBoundingBox(const Iterator &input_b, const Iterator &input_e);

  // Implementations for find functions.

  int FindCandidatesImpl(const CoordinateType   &query,
                         CandidatePriorityQueue &candidates, int max_candidates,
                         int depth = 0) const;
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

template <typename CoordinateType, typename Compare, typename Iterator>
SimplifiedKDTree<CoordinateType, Compare, Iterator>::SimplifiedKDTree(
  int k, const Iterator &input_b, const Iterator &input_e, int depth /*=0*/)
    : points_b(), points_e(), left(nullptr), right(nullptr)
{

  // empty tree
  if (input_b == input_e)
    return;

  num_dims   = k;
  num_points = std::distance(input_b, input_e);
  split_dir  = depth % k;

  assert(num_points > 0);

  SetBoundingBox(input_b, input_e);

  int                 first_dir = split_dir;
  int                 best_dir  = 0;
  int                 split;
  double              best_quality = 0.0;
  std::vector<double> quality(num_dims, 0.0);

  // Any large max iterations value should be fine.
  int iter = 0;
  for (; iter < 500; ++iter)
  {
    if (num_points <= 4)
    {
      split              = num_points;
      quality[split_dir] = 0.0;
      best_dir           = split_dir;
      break;
    }

    split = GetBestSplit(input_b, input_e, split_dir);

    // calculate the quality of this split
    quality[split_dir] = std::min(split, num_points - split);

    // check if this is the best value seen till now.
    if (split_dir == first_dir || quality[split_dir] > best_quality)
    {
      best_quality = quality[split_dir];
      best_dir     = split_dir;
    }

    // check if we need to improve.
    if (quality[split_dir] > 0.25 * num_points)
      break;
    else
    {
      int new_dir = (split_dir + 1) % num_dims;
      if (new_dir == first_dir)
        break;
      split_dir = new_dir;
    }
  }

  assert(iter < 500); // we should converge, if not, something went wrong.

  // base case --- store and return
  if (quality[best_dir] < 0.25 * num_points)
  {
    // make copies of begin and end iterators.
    has_child = false;
    points_b  = input_b;
    points_e  = input_e;

    return;
  }

  // we have subtrees
  if (split_dir != best_dir)
  {
    split     = GetBestSplit(input_b, input_e, best_dir);
    split_dir = best_dir;
  }

  const Iterator &sp  = std::next(input_b, split);
  const Iterator &spm = std::next(input_b, split - 1);

  split_value
    = (GetCoordinate(*sp)[split_dir] + GetCoordinate(*spm)[split_dir]) / 2.0;

  has_child = true;
  left = std::make_shared<SimplifiedKDTree<CoordinateType, Compare, Iterator>>(
    num_dims, input_b,
    std::next(sp), // non-inclusive
    split_dir + 1);
  right = std::make_shared<SimplifiedKDTree<CoordinateType, Compare, Iterator>>(
    num_dims, std::next(sp),
    input_e, // non-inclusive
    split_dir + 1);
}

//-----------------------------------------------------------------------------

template <typename CoordinateType, typename Compare, typename Iterator>
void SimplifiedKDTree<CoordinateType, Compare, Iterator>::SetBoundingBox(
  const Iterator &input_b, const Iterator &input_e)
{

  for (int k = 0; k < num_dims; ++k)
  {

    // input could be an bidirectional contianer type (e.g., map, set, etc.)
    std::vector<CoordinateContainerType> points;

    for (auto it = input_b; it != input_e; ++it)
      points.push_back(GetCoordinate(*it));

    // sort based on co-ordinate.
    Compare compare(k);
    std::sort(points.begin(), points.end(), compare);

    // store bounds for this dimension.
    low.push_back(points.front()[k]);
    upp.push_back(points.back()[k]);
  }
}

//-----------------------------------------------------------------------------

template <typename CoordinateType, typename Compare, typename Iterator>
int SimplifiedKDTree<CoordinateType, Compare, Iterator>::GetBestSplit(
  const Iterator &input_b, const Iterator &input_e, int dir)
{

  assert(num_dims > 0);
  assert(input_b != input_e);

  int input_size = std::distance(input_b, input_e);

  assert(input_size > 0);

  if (input_size <= 4)
    return input_size;

  // input could be an bidirectional contianer type (e.g., map, set, etc.)
  std::vector<CoordinateContainerType> points;

  for (auto it = input_b; it != input_e; ++it)
    points.push_back(GetCoordinate(*it));

  // sort based on co-ordinate.
  Compare compare(dir);
  std::sort(points.begin(), points.end(), compare);

  int left_split  = input_size / 2;
  int right_split = left_split;

  // Find the best median. Make sure that the split points
  // is such that the object at the median and the one before it
  // are not collocated.
  while (left_split > 0)
  {

    const CoordinateContainerType &lsp  = points[left_split];
    const CoordinateContainerType &lspm = points[left_split - 1];
    if (lsp[dir] == lspm[dir])
      left_split--;
    else
      break;
  }

  // Note that since input_size>4 we are certain that right_split
  // can not be 0.
  while (right_split < input_size)
  {

    const CoordinateContainerType &rsp  = points[right_split];
    const CoordinateContainerType &rspm = points[right_split - 1];
    if (rsp[dir] == rspm[dir])
      right_split++;
    else
      break;
  }

  int split
    = (input_size - left_split < right_split) ? left_split : right_split;
  return split;
}

//-----------------------------------------------------------------------------

template <typename CoordinateType, typename Compare, typename Iterator>
int SimplifiedKDTree<CoordinateType, Compare, Iterator>::FindCandidates(
  const CoordinateType &query, std::vector<CoordinateType> &output,
  int max_candidates) const
{

  assert(max_candidates > 0);

  // internal variable
  CandidatePriorityQueue candidates;

  // get all candidates.
  int num_found = FindCandidatesImpl(query, candidates, max_candidates);

  output.clear();
  int counter = 0;
  while (!candidates.empty() and counter < num_found)
  {
    output.push_back(*candidates.top().second);
    candidates.pop();
  }

  return num_found;
}

//-----------------------------------------------------------------------------

template <typename CoordinateType, typename Compare, typename Iterator>
int SimplifiedKDTree<CoordinateType, Compare, Iterator>::FindCandidatesImpl(
  const CoordinateType &query, CandidatePriorityQueue &candidates,
  int max_candidates, int depth) const
{

  // base case
  if (!has_child)
  {

    int k;
    for (k = 0; k < num_dims; ++k)
    {
      if (GetCoordinate(query)[k] >= low[k]
          and GetCoordinate(query)[k] <= upp[k])
      {
        break;
      }
    }

    if (k != num_dims)
      return 0;

    for (Iterator pit = points_b; pit != points_e; ++pit)
    {

      double d = EuclideanDistance(num_dims, GetCoordinate(query),
                                   GetCoordinate(*pit));

      // store this point if we have space.
      if ((int)candidates.size() < max_candidates)
        candidates.push({d, &(*pit)});

      // check if this point is the new "worst"
      if (d < candidates.top().first)
      {
        candidates.pop();
        candidates.push({d, &(*pit)});
      }
    }

    return (int)candidates.size();
  }

  bool side = GetCoordinate(query)[split_dir] <= split_value;

  // see which tree to traverse first
  using TreeType = SimplifiedKDTree<CoordinateType, Compare, Iterator>;
  std::shared_ptr<const TreeType> first_branch  = (side) ? left : right;
  std::shared_ptr<const TreeType> second_branch = (side) ? right : left;

  first_branch->FindCandidatesImpl(query, candidates, max_candidates,
                                   depth + 1);

  double max_dist = ((int)candidates.size() < max_candidates)
                      ? std::numeric_limits<double>::infinity()
                      : candidates.top().first;

  // go to the second branch if we have space left and
  // there is a possiblity of finding good points.
  if ((int)candidates.size() < max_candidates
      or GetCoordinate(query)[split_dir] <= split_value + max_dist)
    second_branch->FindCandidatesImpl(query, candidates, max_candidates,
                                      depth + 1);

  return (int)candidates.size();
}

//-----------------------------------------------------------------------------

template <typename CoordinateType, typename Compare, typename Iterator>
int SimplifiedKDTree<CoordinateType, Compare, Iterator>::FindCandidatesWithin(
  const CoordinateType &query, std::vector<CoordinateType> &output,
  int max_candidates, double distance) const
{

  return 0;
}

//-----------------------------------------------------------------------------

template <typename CoordinateType, typename Compare, typename Iterator>
int SimplifiedKDTree<CoordinateType, Compare, Iterator>::FindCandidatesInBox(
  const CoordinateType        &qmin, /* query bounding-box min */
  const CoordinateType        &qmax, /* query bounding-box max*/
  std::vector<CoordinateType> &output, int max_candidates) const
{

  return 0;
}

//-----------------------------------------------------------------------------

/*
template<typename CoordinateType, typename Compare, typename Iterator>
void
SimplifiedKDTree<CoordinateType,Compare,Iterator>::
Print(std::ostream &os, int depth) const
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
*/

//-----------------------------------------------------------------------------

template <typename CoordinateContainerType> class DefaultComparator
{
  int dir;

public:
  DefaultComparator(int d) : dir(d)
  {
  }
  bool operator()(const CoordinateContainerType &lhs,
                  const CoordinateContainerType &rhs) const
  {
    return lhs[dir] < rhs[dir];
  }
};

//-----------------------------------------------------------------------------

// CTAD (Class Template Argument Deduction) for compiler refer
// 1. https://en.cppreference.com/w/cpp/language/class_template_argument_deduction
// 2. https://stackoverflow.com/questions/40951697/what-are-template-deduction-guides-and-when-should-we-use-them

// Deduce everything from the iterator type.
template <typename Iterator>
SimplifiedKDTree(int, const Iterator &, const Iterator &, int)
  -> SimplifiedKDTree<
    typename std::iterator_traits<Iterator>::value_type,
    /*DefaultComparator<typename std::iterator_traits<Iterator>::value_type>,*/
    DefaultComparator<
      std::decay_t<decltype(GetCoordinate(*std::declval<Iterator>()))>>,
    Iterator>;

template <typename Iterator>
SimplifiedKDTree(int, const Iterator &, const Iterator &) -> SimplifiedKDTree<
  typename std::iterator_traits<Iterator>::value_type,
  /*DefaultComparator<typename std::iterator_traits<Iterator>::value_type>,*/
  DefaultComparator<
    std::decay_t<decltype(GetCoordinate(*std::declval<Iterator>()))>>,
  Iterator>;

//-----------------------------------------------------------------------------

#endif
