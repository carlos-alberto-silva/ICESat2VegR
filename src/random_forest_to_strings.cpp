#include <R.h>
#include <Rcpp.h>
#include "random_forest_to_strings.h"

using namespace Rcpp;

CharacterVector build_forest(List rf)
{
  int ntree = rf["ntree"];
  CharacterVector output(ntree);

  for (int ii = 0; ii < ntree; ii++)
  {
    std::stringstream temp;
    temp << "1) root 9999 9999 ";
    visitNode(temp, rf, 1, 1, 1, ii);
    output(ii) = temp.str();
  }

  return output;
}

void visitNode(std::stringstream &temp, List rf, int idx, int previdx = 1, int depth = 1, int tree = 1)
{
  List rfobj = rf["forest"];
  // Accessing columns from DataFrame
  IntegerMatrix leftDaughter = rfobj["leftDaughter"];
  int left_idx = leftDaughter(idx - 1, tree);

  NumericMatrix importance = rf["importance"];
  CharacterVector importanceNames = rownames(importance);
  IntegerMatrix bestvar = rfobj["bestvar"];
  int var_idx = bestvar(idx - 1, tree);
  auto split_var = var_idx > 0 ? Rcpp::as<std::string>(importanceNames[var_idx - 1]) : "";

  NumericMatrix xbestsplit = rfobj["xbestsplit"];
  double split_point = xbestsplit(idx - 1, tree);

  IntegerMatrix rightDaughter = rfobj["rightDaughter"];
  int right_idx = rightDaughter(idx - 1, tree);

  NumericMatrix nodepred = rfobj["nodepred"];
  double node_prediction = nodepred(idx - 1, tree);

  // Check if its a leaf node
  if (left_idx == 0)
  {
    temp << node_prediction << " *";
    return;
  }
  else
  {
    temp << node_prediction;
  }

  // Rule for the current node
  temp << std::endl
       << std::string(depth, ' ')
       << previdx * 2
       << ") "
       << split_var
       << ((left_idx % 2) == 0 ? " <= " : " > ")
       << split_point << " 0 0 ";

  // Recursively visit left_idx and right_idx nodes
  visitNode(temp, rf, left_idx, previdx * 2, depth + 1, tree);
  temp << std::endl
       << std::string(depth, ' ')
       << previdx * 2 + 1
       << ") "
       << split_var
       << ((right_idx % 2) == 0 ? " <= " : " > ")
       << split_point << " 0 0 ";

  visitNode(temp, rf, right_idx, previdx * 2 + 1, depth + 1, tree);
  return;
}
