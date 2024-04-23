#ifndef RANDOM_FOREST_TO_STRINGS_H_
#define RANDOM_FOREST_TO_STRINGS_H_

#include <Rcpp.h>

Rcpp::CharacterVector build_forest(Rcpp::List rf);
void visitNode(std::stringstream &, Rcpp::List, int, int, int, int);

#endif // RANDOM_FOREST_TO_STRINGS_H_