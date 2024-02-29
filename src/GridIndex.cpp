#include <Rcpp.h>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <inttypes.h>
#include <algorithm>
#include <GridIndex.h>

using namespace Rcpp;

void GridIndex::buildIndex()
{
  mapIndex.clear();
  for (R_xlen_t i = 0; i < X.size(); ++i)
  {
    uint64_t key = convertToIndex(X[i], Y[i]);
    mapIndex[key].push_back(i);
  }
}

Rcpp::XPtr<GridIndex> GridIndex::getPointer()
{
  return Rcpp::XPtr<GridIndex>(this);
}

GridIndex::GridIndex(const NumericVector &xValues, const NumericVector &yValues, const double factor = 1.0)
    : X(xValues), Y(yValues), factor(factor)
{
  NumericVector rngY = Rcpp::range(Y);
  NumericVector rngX = Rcpp::range(X);
  yMin = rngY[0];
  yMax = rngY[1];
  xMin = rngX[0];
  xMax = rngX[1];

  // Calculate bitShift
  bitShift = static_cast<uint8_t>(std::log2((yMax - yMin) / factor) + 1);

  // Build index
  buildIndex();
}

uint64_t GridIndex::convertToIndex(double x, double y)
{
  return static_cast<uint64_t>((x - xMin) / factor) << bitShift | static_cast<uint64_t>((y - yMin) / factor);
}

IntegerVector GridIndex::pointsWithinIndex(uint64_t index)
{
  return mapIndex[index] + 1;
}

IntegerVector GridIndex::searchFixedRadius(const double x, const double y, const double dist)
{
  IntegerVector result;
  double xmin = std::max((x - dist), xMin);
  double xmax = std::min((x + dist + factor / 2), xMax);
  double ymin = std::max((y - dist), yMin);
  double ymax = std::min((y + dist + factor / 2), yMax);
  for (double curY = ymin; curY <= ymax; curY += factor)
  {
    for (double curX = xmin; curX <= xmax; curX += factor)
    {
      uint64_t idx = convertToIndex(curX, curY);
      auto it = mapIndex.find(idx);
      if (it != mapIndex.end())
      {
        for (int index : it->second)
        {
          double distance = std::sqrt((X[index] - x) * (X[index] - x) + (Y[index] - y) * (Y[index] - y));
          if (distance <= dist)
          {
            result.push_back(index);
          }
        }
      }
    }
  }
  return result;
}

RCPP_MODULE(grid_index_module)
{
  class_<GridIndex>("GridIndex")

      .constructor<NumericVector, NumericVector, double>()

      .method("convertToIndex", &GridIndex::convertToIndex)
      .method("pointsWithinIndex", &GridIndex::pointsWithinIndex)
      .method("searchFixedRadius", &GridIndex::searchFixedRadius)
      .method("getPointer", &GridIndex::getPointer);
}
