#include <iostream>
#include <fstream>
#include <stdlib.h> 
#include <memory>
#include "grid.h"
#include "initialiserFunctions.h"



int main(int argc, char* argv[])
{
 
  int l = atoi(argv[1]);
  
  int n = (1 << l) + 1;
  
  Grid A(n, boundary_function_zero);
  for(int i = 0; i< 50; ++i)
  {
    A.rb_gauss_seidel_relaxation();
  }
  A.print_v();
  

}

