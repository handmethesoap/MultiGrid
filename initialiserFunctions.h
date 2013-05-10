#ifndef INITIALISER_FUNCTIONS_H
#define INITIALISER_FUNCTIONS_H

double boundary_function_zero(int grid_point_row, int grid_point_column, int n);
double sin_boundary(int grid_point_row, int grid_point_column, int n);
double f_initialiser_function_sin(int grid_point_row, int grid_point_column, int n);
double v_initialiser_function_sin_2(int grid_point_row, int grid_point_column, int n);
double v_initialiser_function_zero(int grid_point_row, int grid_point_column, int n);
double f_initialiser_function_zero(int grid_point_row, int grid_point_column, int n);
double exact_solution( int grid_point_row, int grid_point_column, int n);

#endif