#include <Rcpp.h>
#include <random>
#include <ANNIndex.h>

using namespace Rcpp;

IntegerVector sampleMinDistanceRcpp(NumericVector &x, NumericVector &y, const double radius, const int sampleSize, Environment indexPtr, const int seed)
{
  // Set seed for the std random number generator
  srand(seed);
  int nrows = x.length();

  LogicalVector tabooFlag(nrows);
  IntegerVector output(sampleSize);

  // Create a vector of indices from 0 to nrows - 1
  IntegerVector indices = Range(0, nrows - 1);

  SEXP objptr = indexPtr[".pointer"];
  IndexInterface *genericIndex = (IndexInterface *)R_ExternalPtrAddr(objptr);

  // Generate a vector of random numbers
  std::shuffle(indices.begin(), indices.end(), std::mt19937(seed));

  int ii = 0;
  for (int idx : indices)
  {
    // Skip if idx is taboo
    if (tabooFlag[idx] == true)
      continue;

    // Mark as taboo and accept sample
    tabooFlag[idx] = true;
    output[ii++] = idx + 1;

    // Did we reach sampleSize?
    if (ii >= sampleSize)
    {
      break;
    }

    // Add points within radius to tabooList
    IntegerVector pointsWithin = genericIndex->searchFixedRadius(x[idx], y[idx], radius);
    for (int tabooIdx : pointsWithin)
    {
      tabooFlag[tabooIdx] = true;
    }
  }

  IntegerVector final = output[Rcpp::Range(0, ii - 1)];
  return final;
}

RCPP_MODULE(icesat2_module)
{
  function("sampleMinDistanceRcpp", &sampleMinDistanceRcpp);

  class_<ANNIndex>("ANNIndex")

      .constructor<NumericVector, NumericVector>()
      .method("searchFixedRadius", ANNIndex::searchFixedRadius);

}
