#include "initialiserFunctions.h"
#include "math.h"
#include <iostream>
const double pi = 3.14159265;

double* sin_boundary(int n)
{
  double* init = new double[n*n]();
  
  for( int it = 0; it < n; ++it)
  {
    init[it] = 0;
    init[n*n - n + it] = sin(pi*it/(n-1))*sinh(pi);
    init[it*n] = 0;
    init[it*n + n - 1] = 0;
  }
  
  return init;
}

double* initialiser_function_zero(int n)
{
    double* init = new double[n*n]();
  
  for( int it = 0; it < n*n; ++it)
  {
    init[it] = 0;
  }
  
  return init;
}

double* numbered_initialiser(int n)
{
    double* init = new double[n*n]();
  
  for( int it = 0; it < n*n; ++it)
  {
    init[it] = it;
  }
  
  return init;
}