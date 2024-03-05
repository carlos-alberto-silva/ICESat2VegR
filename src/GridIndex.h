#ifndef GRID_INDEX_H
#define GRID_INDEX_H

#include <vector>
#include <unordered_map>
#include <Rcpp.h>
#include <IndexInterface.h>

using namespace Rcpp;

class GridIndex : public IndexInterface {
private:
    double grid_size;
    std::unordered_map<int, std::unordered_map<int, IntegerVector>> grid;
    NumericVector x_coords;
    NumericVector y_coords;

public:
    GridIndex(const NumericVector& x, const NumericVector& y, const double grid_size);
    void addToGrid(double x, double y, int index);
    virtual Rcpp::IntegerVector searchFixedRadius(const double x, const double y, const double radius);
};

#endif  // GRID_INDEX_H