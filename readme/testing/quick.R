# source("C:/Users/caiohamamura/src/r/rICESat2Veg/readme/testing/quick.R")

dllPath <- system.file(sprintf("libs/%s", gsub("86-", "", Sys.info()["machine"])), package = "RANN2")
Sys.setenv(PKG_LIBS = sprintf("-L\"%s\" -lRANN2", dllPath))


Rcpp::cppFunction('
NumericVector findRadius(NumericVector x, NumericVector y, const double radius, const int sampleSize){
  Rcpp::RNGScope rngScope;

  const int n_dimensions = 2;
  int nrows = x.length();

  ANNkd_tree	*the_tree;	// Search structure
  ANNpointArray data_pts 	= annAllocPts(nrows, n_dimensions);
  ANNidxArray nn_idx 		= new ANNidx[nrows];		// Allocate near neigh indices
  ANNdistArray dists 		= new ANNdist[nrows];		// Allocate near neighbor dists
  int n_points;

  // now construct the points
  for(int i = 0; i < nrows; i++)
  {
    data_pts[i][0]=x[1];
    data_pts[i][1]=y[i];
  }


  the_tree = new ANNkd_tree( data_pts, nrows, n_dimensions);
  Rcout << "Ok" << std::endl;

  NumericVector result(sampleSize);
  int tabooSize = 0;
  IntegerVector tabooList(x.length());
  Rcout << "Ok1" << std::endl;

  for (int ii = 0; ii < sampleSize; ii++) {
    int sampleIndex = (int)floor(R::runif(0, nrows--));

    int correctedIndex = countGreaterEq(tabooList, sampleIndex, tabooSize);

    result[ii] = correctedIndex;

    if (correctedIndex > nrows) {
      Rprintf("Cannot find anymore samples beond the provided radius!");
      return result;
    }

    tabooList[tabooSize++] = correctedIndex;

    // Add points around to tabooList
    Rcout << "Ok" << ii << std::endl;
    ANNpoint pt_query = annAllocPt(n_dimensions);
    pt_query[0] = x[correctedIndex];
    pt_query[1] = y[correctedIndex];

    n_points = the_tree->annkFRSearch(pt_query, pow(radius, 2), nrows, nn_idx, dists);

    for (int jj = 0; jj < n_points; jj++) {
     tabooList[tabooSize++] = nn_idx[jj];
    }

    if (tabooSize == x.length()) {
      Rprintf("Cannot find anymore samples beyond the provided radius!");
      return result;
    }
  }
  return result;
}',
  depends = c("RANN2"),
  includes = '
  #include <ANN/ANN.h>
  int countGreaterEq(IntegerVector x, int value, int pos) {
  int countContainer = 0;
  int initialValue = value;
  LogicalVector intervalContainer(pos - 1);
  for (int ii = 0; ii < pos; ii++) {
    if (x[ii] <= value) {
      value++;
    // As value increments, test if it can be contained by
    // the most incremented possible value after all pos
    // are evaluated
    } else if (x[ii] <= (initialValue + pos)) {
      Rprintf("x:%d init:%d value:%d\\n", x[ii], initialValue, value);
      intervalContainer[x[ii] - initialValue - 1] = true;
      countContainer++;
    }
  }

  int check = initialValue + 1;
  int ii = 0;
  while (check <= value || ii >= pos) {
    Rprintf("Check:%d value:%d container:%d\\n", check, value, intervalContainer[ii]);
    if (intervalContainer[ii] == true) {
      value++;
    }
    check++;
    ii++;
  }
  return value;
}
  ')
