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
  GridIndex gridIndex(x, y, (abs(rng[1] - rng[0]) / 1000));

  // Create a vector of indices from 0 to nrows - 1
  IntegerVector indices = Range(0, nrows - 1);

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
    if (ii >= sampleSize) {
      break;  
    }

    // Add points within radius to tabooList
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
