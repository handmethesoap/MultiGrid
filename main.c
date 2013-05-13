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
  
  Grid A(n,  v_initialiser_function_zero, f_initialiser_function_zero);
//   for(int i = 0; i< 20000; ++i)
//   {
//     A.rb_gauss_seidel_relaxation();
//   }

  A.calculate_residual();
  
  
  A.print_v();
  A.print_f();
  A.print_r();
  Grid B((n+1)/2, A.fw_restrict());
  //B.print_v();

  return 0;

}

