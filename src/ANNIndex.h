#ifndef ANN_KD_TREE_H
#define ANN_KD_TREE_H

#include <ANN/ANN.h>
#include <Rcpp.h>
#include <IndexInterface.h>

using namespace Rcpp;

class ANNIndex : public IndexInterface
{
private:
  ANNkd_tree *tree;
  ANNidxArray nn_idx; // Allocate near neigh indices
  ANNdistArray dists; // Allocate near neighbor dists
  ANNpointArray dataPts; // Allocate near neighbor dists
  int nPoints;
  int nDimensions;

public:
  ANNIndex(NumericVector x, NumericVector y, int dimensions = 2);
  ~ANNIndex() override;

  IntegerVector searchFixedRadius(const double x, const double y, const double radius) override;
};

#endif