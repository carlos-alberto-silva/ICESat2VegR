#include <GridIndex.h>
#include <Rcpp.h>
#include <random>

using namespace Rcpp;

//' Sample within a radius distance
//' @param x NumericVector
//' @param y NumericVector
//' @param radius double
//' @param sampleSize double
//' @export
// [[Rcpp::export]]
IntegerVector findRadius(NumericVector x, NumericVector y, const double radius, const int sampleSize)
{
  // Set seed for the std random number generator
  int seed = (int)R::runif(0, std::numeric_limits<int>::max());
  int nrows = x.length();

  LogicalVector tabooFlag(nrows);
  IntegerVector output(sampleSize);
  // Rcout << "Ok1" << std::endl;

  NumericVector rng = range(x);
  GridIndex gridIndex(x, y, (abs(rng[1] - rng[0]) / 1.0));

  // Create a vector of indices from 0 to nrows - 1
  IntegerVector indices = Range(0, nrows - 1);

  // Generate a vector of random numbers
  std::shuffle(indices.begin(), indices.end(), std::mt19937(seed));

  int ii = 0;
  for (int idx : indices)
  {
    if (tabooFlag[idx] == true)
      continue;

    tabooFlag[idx] = true;
    output[ii++] = idx;

    if (ii >= sampleSize) {
      break;  
    }

    // Add points around to tabooList
    IntegerVector pointsWithin = gridIndex.searchFixedRadius(x[idx], y[idx], radius);
    for (int tabooIdx : pointsWithin)
    {
      if (tabooFlag[tabooIdx] == false)
      {
        tabooFlag[tabooIdx] = true;
      }
    }
  }

  IntegerVector final = output[Rcpp::Range(0, ii - 1)];
  return final;
}
