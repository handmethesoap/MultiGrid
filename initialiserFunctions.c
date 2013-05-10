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

double v_initialiser_function_sin_2(int grid_point_row, int grid_point_column, int n)
{
  int k[9] = {1,8,16,24,32,40,48,56,64};
  static int i = 1; 
  double grid_value = sin(k[i]*pi*grid_point_column/(n-1))*sin(k[i]*pi*grid_point_row/(n-1));
  ++i;
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