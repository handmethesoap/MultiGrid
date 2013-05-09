#include <iostream>
#include <stdlib.h> 
#include <memory>
#include "grid.h"
#include "initialiserFunctions.h"



int main(int argc, char* argv[])
{
 
  int l = atoi(argv[1]);
  
  int n = (1 << l) + 1;
  
  Grid A(n, boundary_function_zero);
  A.print_v();
  A.jacobi_relaxation();
  
  //A.set_initial(v_initialiser_function, f_initialiser_function);
  //A.set_boundary(sin_boundary);
  
  A.print_v();
  A.jacobi_relaxation();
  A.print_v();
  A.jacobi_relaxation();
  A.print_v();
  A.jacobi_relaxation();
  for( int it_col = 0; it_col < 100; ++ it_col )
  {
    A.jacobi_relaxation();
  }
  A.print_v();
  
  
  
  
  
  
  
  
  
  
  
  
  return 0;
}

