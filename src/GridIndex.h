#ifndef GRID_INDEX_H
#define GRID_INDEX_H

#include <vector>
#include <unordered_map>
#include <inttypes.h>
#include <Rcpp.h>

using namespace Rcpp;

class GridIndex {
private:
    NumericVector X;
    NumericVector Y;
    double xMin;
    double yMin;
    double xMax;
    double yMax;
    uint8_t bitShift;
    double factor;
    std::unordered_map<uint64_t, IntegerVector> mapIndex;
    void buildIndex();

public:
    GridIndex(const NumericVector& xValues, const NumericVector& yValues, const double factor);

    IntegerVector searchFixedRadius(const double x, const double y, const double dist);
    uint64_t convertToIndex(double x, double y);
    IntegerVector pointsWithinIndex(uint64_t index);
};

#endif // GRID_INDEX_H