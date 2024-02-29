#include <IndexInterface.h>

IndexInterface::~IndexInterface() {}
Rcpp::IntegerVector IndexInterface::searchFixedRadius(const double x, const double y, const double radius) {
      Rcpp::Rcout << "Default implementation in IndexInterface\n";
      return Rcpp::IntegerVector(0);
}