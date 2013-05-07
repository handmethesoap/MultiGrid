#include <iostream>
#include <stdlib.h> 
#include <memory>
#include "grid.h"
#include "initialiserFunctions.h"



int main(int argc, char* argv[])
{
 
  int l = atoi(argv[1]);
  
  int n = (1 << l) + 1;
  
  Grid A(n, sin_boundary);
  
  //A.set_initial(v_initialiser_function, f_initialiser_function);
  //A.set_boundary(sin_boundary);
  
  A.print_v();
  A.print_f();
  
  return 0;
}

