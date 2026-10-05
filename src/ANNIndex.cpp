#include "ANNIndex.h"

ANNIndex::ANNIndex(NumericVector x, NumericVector y, int dimensions)
{
  nDimensions = dimensions;
  nPoints = x.length();
  dataPts = annAllocPts(nPoints, dimensions);

  for (int i = 0; i < nPoints; i++)
  {
    dataPts[i][0] = x[i];
    dataPts[i][1] = y[i];
  }
  tree = new ANNkd_tree(dataPts, nPoints, dimensions);

  nn_idx = new ANNidx[nPoints];
  dists = new ANNdist[nPoints];
}

ANNIndex::~ANNIndex()
{
  delete [] nn_idx;
  delete [] dists;
  delete tree;
  annDeallocPts(dataPts);
  annClose();
}

IntegerVector ANNIndex::searchFixedRadius(const double x, const double y, const double radius)
{
  ANNpoint pt_query = annAllocPt(2);
  
  pt_query[0] = x;
  pt_query[1] = y;
  // ANN expects the squared search radius; callers provide a distance.
  int count = tree->annkFRSearch(pt_query, radius * radius, nPoints, nn_idx, dists);
  
  annDeallocPt(pt_query);
  int *idxx = (int*)nn_idx;
  IntegerVector result(idxx, idxx + count);
  return result;
}
