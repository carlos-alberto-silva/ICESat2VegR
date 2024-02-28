// [[Rcpp::depends(ICESat2VegR)]]
#include <GridIndex.cpp>
#include <Rcpp.h>


using namespace Rcpp;

// [[Rcpp::export]]
int countGreaterEq(IntegerVector tabooVector, int value, int pos)
{
  // Rprintf("%d %d\nVector:", value, pos);
  // for (int ii = 0; ii < pos; ii++) {
  //   Rprintf(" %d", tabooVector[ii]);
  // }
  // Rprintf("\n");

  int countContainer = 0;
  int initialValue = value;
  int intervalSize = pos > 1 ? pos - 1 : 1;
  LogicalVector intervalContainer(intervalSize);


  int maxTaboo = tabooVector.length();
  for (int ii = 0; ii < pos && ii < maxTaboo; ii++)
  {
    if (tabooVector[ii] <= value)
    {
      value++;
      // As value increments, test if it can be contained by
      // the most incremented possible value after all pos
      // are evaluated
    }
    else if (tabooVector[ii] <= (initialValue + pos))
    {
      int jj = std::min(intervalSize - 1,  tabooVector[ii] - initialValue - 1);
      intervalContainer[jj] = true;
      countContainer++;
    }
  }

  if (countContainer == 0) {
    return value;
  }
  int check = initialValue + 1;
  int ii = 0;
  while (check <= value && ii < (pos - 1) && ii < intervalSize)
  {
    if (intervalContainer[ii] == true)
    {
      value++;
    }
    check++;
    ii++;
  }
  return value;
}



// [[Rcpp::export]]
IntegerVector findRadius3(NumericVector x, NumericVector y, const double radius, const int sampleSize)
{
  int nrows = x.length();

  int tabooSize = 0;
  LogicalVector tabooFlag(x.length());
  IntegerVector tabooList(x.length());
  IntegerVector output(sampleSize);
  // Rcout << "Ok1" << std::endl;

  NumericVector rng = range(x);
  GridIndex gridIndex(x, y, (abs(rng[1] - rng[0]) / 1000.0));

  for (int ii = 0; ii < sampleSize; ii++)
  {
    int sampleIndex = (int)floor(R::runif(0, nrows--));

    int correctedIndex = countGreaterEq(tabooList, sampleIndex, tabooSize);
    // int correctedIndex = sampleIndex;

    output[ii] = correctedIndex;

    if (correctedIndex > nrows)
    {
      // Rprintf("Cannot find anymore samples beyond the provided radius!");
      break;
    }

    if (tabooFlag[correctedIndex] == false)
    {
      tabooFlag[correctedIndex] = true;
      tabooList[tabooSize++] = correctedIndex;
    }

    // Add points around to tabooList
    // Rcout << "Ok" << ii << std::endl;
    IntegerVector pointsWithin = gridIndex.searchFixedRadius(x[correctedIndex], y[correctedIndex], radius);

    for (int jj = 0; jj < pointsWithin.length(); jj++)
    {
      if (tabooFlag[pointsWithin[jj]] == false)
      {
        tabooFlag[pointsWithin[jj]] = true;
        tabooList[tabooSize++] = pointsWithin[jj];
      }
    }

    if (tabooSize == x.length())
    {
      // Rprintf("Cannot find anymore samples beyond the provided radius!");
      break;
    }
  }
  return output;
}


