#ifndef COUNTGREATEREQ_CPP
#define COUNTGREATEREQ_CPP
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
int countGreaterEq(IntegerVector x, int value, int pos)
{
  // Rprintf("%d %d\nVector:", value, pos);
  // for (int ii = 0; ii < pos; ii++) {
    // Rprintf(" %d", x[ii]);
  // }
  // Rprintf("\n");

  int countContainer = 0;
  int initialValue = value;
  LogicalVector intervalContainer(pos > 1 ? pos - 1 : 1);
  for (int ii = 0; ii < pos; ii++)
  {
    if (x[ii] <= value)
    {
      value++;
      // As value increments, test if it can be contained by
      // the most incremented possible value after all pos
      // are evaluated
    }
    else if (x[ii] <= (initialValue + pos))
    {
      intervalContainer[x[ii] - initialValue - 1] = true;
      countContainer++;
    }
  }

  if (countContainer == 0) return value;
  int check = initialValue + 1;
  int ii = 0;
  while (check <= value && ii < (pos - 1))
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

#endif