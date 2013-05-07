#include "initialiserFunctions.h"
#include "math.h"

double boundary_function(int grid_point_row, int grid_point_column, int n)
{
  double grid_value = 4.0;
  return grid_value;
}

double sin_boundary(int grid_point_row, int grid_point_column, int n)
{
  const double pi = 3.14159265;
  double grid_value = sin(pi*grid_point_row/(n-1))*sinh(pi*grid_point_column/(n-1));
  return grid_value;
}

double v_initialiser_function(int grid_point_row, int grid_point_column, int n)
{
  double grid_value = 2.0;
  return grid_value; 
}

double f_initialiser_function(int grid_point_row, int grid_point_column, int n)
{
  double grid_value = 1.0;
  return grid_value; 
}