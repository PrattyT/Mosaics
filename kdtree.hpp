/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <algorithm>
#include <utility>

using namespace std;

template <int Dim>
bool KDTree<Dim>::smallerDimVal(const Point<Dim> &first,
                                const Point<Dim> &second, int curDim) const {
  if (first[curDim] == second[curDim])
    return first < second;

  return first[curDim] < second[curDim];
}

template <int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim> &target,
                                const Point<Dim> &currentBest,
                                const Point<Dim> &potential) const {

  double currDistance = distanceSq(target, currentBest);
  double potentialDistance = distanceSq(target, potential);
  if (currDistance == potentialDistance)
    return potential < currentBest;

  return potentialDistance < currDistance;
}

template <int Dim>
double KDTree<Dim>::distanceSq(const Point<Dim> &left,
                               const Point<Dim> &right) const {
  double ret = 0;
  for (int i = 0; i < Dim; i++) {
    ret += pow(left[i] - right[i], 2);
  }

  return ret;
}

template <int Dim> KDTree<Dim>::KDTree(const vector<Point<Dim>> &newPoints) {
  for (Point<Dim> p : newPoints) {
    points.push_back(p);
  }
  root = buildTree(0, 0, points.size() - 1);
}

template <int Dim> void KDTree<Dim>::_copy(const KDTree<Dim> &other) {
  for (Point<Dim> p : other.points) {
    points.push_back(p);
  }
  root = buildTree(0, 0, points.size() - 1);
}

template <int Dim> KDTree<Dim>::KDTree(const KDTree<Dim> &other) {
  _copy(other);
}

template <int Dim>
typename KDTree<Dim>::KDTreeNode *KDTree<Dim>::buildTree(int dim, int left,
                                                         int right) {
  if (left == right)
    return new KDTreeNode(points[left]);
  if (left < right) {
    int middle = (left + right) / 2;
    select(left, right, middle, dim);
    KDTreeNode *node = new KDTreeNode(points[middle]);
    node->left = buildTree((dim + 1) % Dim, left, middle - 1);
    node->right = buildTree((dim + 1) % Dim, middle + 1, right);
    return node;
  }

  return NULL;
}

template <int Dim> int KDTree<Dim>::select(int l, int r, int k, int dim) {
  if (l == r)
    return l;

  int pivot = l + 1;
  pivot = partition(l, r, pivot, dim);
  if (k == pivot)
    return pivot;
  else if (k < pivot)
    return select(l, pivot - 1, k, dim);
  else
    return select(pivot + 1, r, k, dim);
}

template <int Dim>
int KDTree<Dim>::partition(int l, int r, int pivot, int dim) {
  Point<Dim> value = points[pivot];
  swap(pivot, r);
  int storeIndex = l;
  for (int i = l; i < r; i++) {
    if (smallerDimVal(points[i], value, dim)) {
      swap(storeIndex, i);
      storeIndex++;
    }
  }
  swap(r, storeIndex);
  return storeIndex;
}

template <int Dim> void KDTree<Dim>::swap(int l, int r) {
  Point<Dim> temp = points[l];
  points[l] = points[r];
  points[r] = temp;
}

template <int Dim>
const KDTree<Dim> &KDTree<Dim>::operator=(const KDTree<Dim> &rhs) {

  if (this == &rhs)
    return;

  _destroy();
  _copy(rhs);

  return *this;
}

template <int Dim> void KDTree<Dim>::_destroy() { destroy(root); }

template <int Dim> void KDTree<Dim>::destroy(KDTreeNode *node) {
  if (node == NULL)
    return;

  destroy(node->left);
  destroy(node->right);

  delete node;
}

template <int Dim> KDTree<Dim>::~KDTree() { _destroy(); }

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim> &query) const {
  /**
   * @todo Implement this function!
   */
  return findNearestNeighbor(query, 0, root);
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim> &query, int dim,
                                            KDTreeNode *node) const {
  if (node == NULL) {
    return Point<Dim>();
  }
  if (node->left == NULL && node->right == NULL)
    return node->point;

  bool wentLeft = false;
  Point<Dim> nearest;
  if (smallerDimVal(query, node->point, dim)) {
    nearest = findNearestNeighbor(query, (dim + 1) % Dim, node->left);
    wentLeft = true;
  } else
    nearest = findNearestNeighbor(query, (dim + 1) % Dim, node->right);

  if (shouldReplace(query, nearest, node->point))
    nearest = node->point;

  double radius = distanceSq(query, nearest);
  double splitDist = pow(node->point[dim] - query[dim], 2);

  // check if theres a need to search the other subtree.
  Point<Dim> tempNearest;
  if (radius >= splitDist) {

    if (wentLeft)
      tempNearest = findNearestNeighbor(query, (dim + 1) % Dim, node->right);
    else
      tempNearest = findNearestNeighbor(query, (dim + 1) % Dim, node->left);
    if (shouldReplace(query, nearest, tempNearest))
      nearest = tempNearest;
  }

  return nearest;
}
