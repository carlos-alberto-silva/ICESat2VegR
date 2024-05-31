#include <Rcpp.h>
#include <seq_lens_simplify.h>

using namespace Rcpp;

IntegerVector seq_lens_simplify(IntegerVector from, IntegerVector length_out)
{
  IntegerVector output(sum(length_out));
  int pos = 0;
  for (int ii = 0; ii < from.length(); ii++)
  {
    for (int jj = from[ii]; jj < from[ii] + length_out[ii]; jj++)
    {
      output[pos++] = jj;
    }
  }
  return output;
}