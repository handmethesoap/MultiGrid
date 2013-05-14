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

double* sin_fxn(int n)
{
  double* init = new double[n*n]();
  
  for( int it_1 = 0; it_1 < n; ++it_1)
  {
      for( int it_2 = 0; it_2 < n; ++it_2)
      {
	init[it_1*n + it_2] = sin(pi*it_2/(n-1))*sinh(pi*it_1/(n-1));
      }
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