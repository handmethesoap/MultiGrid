#include "initialiserFunctions.h"
#include "math.h"
#include <iostream>
const double pi = 3.14159265;

double boundary_function_zero(int grid_point_row, int grid_point_column, int n)
{
  double grid_value = 0.0;
  return grid_value;
}

double sin_boundary(int grid_point_row, int grid_point_column, int n)
{
  double grid_value = sin(pi*grid_point_column/(n-1))*sinh(pi*grid_point_row/(n-1));
  return grid_value;
}

double f_initialiser_function_sin(int grid_point_row, int grid_point_column, int n)
{
  double grid_value = sin(pi*grid_point_column/(n-1))*sin(pi*grid_point_row/(n-1));
  return grid_value; 
}

double v_initialiser_function_sin_2(int grid_point_row, int grid_point_column, int n, int k)
{

  double grid_value = sin(k*pi*grid_point_column/(n-1))*sin(k*pi*grid_point_row/(n-1));

  return grid_value; 
}

double v_initialiser_function_zero(int grid_point_row, int grid_point_column, int n)
{
  double grid_value = 0.0;
  return grid_value; 
}

double f_initialiser_function_zero(int grid_point_row, int grid_point_column, int n)
{
  double grid_value = 0.0;
  return grid_value; 
}

double exact_solution( int grid_point_row, int grid_point_column, int n)
{
  double grid_value = sin(pi*grid_point_column/(n-1))*sin(pi*grid_point_row/(n-1))/(2*pi*pi);
  return grid_value;
}