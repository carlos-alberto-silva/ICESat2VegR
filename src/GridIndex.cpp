#include <GridIndex.h>
#include <unordered_map>
#include <Rcpp.h>

using namespace Rcpp;

GridIndex::GridIndex(const NumericVector& x, const NumericVector& y, double grid_size)
{
  if (x.size() != y.size())
  {
    REprintf("Error: The size of x and y coordinates must be the same.\n");
    return;
  }

  this->grid_size = grid_size;
  this->x_coords = x;
  this->y_coords = y;

  for (size_t i = 0; i < x.size(); ++i)
  {
    addToGrid(x[i], y[i], i);
  }
}

void GridIndex::addToGrid(double x, double y, int index)
{
  int grid_x = static_cast<int>(x / grid_size);
  int grid_y = static_cast<int>(y / grid_size);
  grid[grid_x][grid_y].push_back(index);
}

IntegerVector GridIndex::searchFixedRadius(const double x, const double y, const double radius)
{
  IntegerVector result;
  int grid_x = static_cast<int>(x / grid_size);
  int grid_y = static_cast<int>(y / grid_size);
  int radiusCells = static_cast<int>(ceil(radius / grid_size));

  for (int i = -radiusCells; i <= radiusCells; ++i)
  {
    for (int j = -radiusCells; j <= radiusCells; ++j)
    {
      int neighbor_x = grid_x + i;
      int neighbor_y = grid_y + j;

      if (grid.count(neighbor_x) && grid[neighbor_x].count(neighbor_y))
      {
        for (int index : grid[neighbor_x][neighbor_y])
        {
          double x_diff = x - x_coords[index];
          double y_diff = y - y_coords[index];
          double distance = sqrt(x_diff * x_diff + y_diff * y_diff);

          if (distance <= radius)
          {
            result.push_back(index);
          }
        }
      }
    }
  }

  return result;
}