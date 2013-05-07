#include "initialiserFunctions.h"

double boundary_function(int grid_point_row, int grid_point_column)
{
  double grid_value = 4.0;
  return grid_value;
}

double v_initialiser_function(int grid_point_row, int grid_point_column)
{
  double grid_value = 2.0;
  return grid_value; 
}

double f_initialiser_function(int grid_point_row, int grid_point_column)
{
  double grid_value = 1.0;
  return grid_value; 
}