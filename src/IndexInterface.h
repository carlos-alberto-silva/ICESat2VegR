#ifndef INDEX_INTERFACE_H
#define INDEX_INTERFACE_H

#include <Rcpp.h>

class IndexInterface {
public:
    virtual ~IndexInterface() = 0;
    virtual Rcpp::IntegerVector searchFixedRadius(const double x, const double y, const double radius);
};

#endif
