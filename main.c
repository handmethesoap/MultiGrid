#include <iostream>
#include <fstream>
#include <stdlib.h> 
#include <memory>
#include "grid.h"
#include "initialiserFunctions.h"



int main(int argc, char* argv[])
{
 
  int l = atoi(argv[1]);
  
  int grid_size = (1 << l) + 1;
  
  Grid **Grids = new Grid*[l];
  
  Grids[1] = new Grid(grid_size,  sin_boundary(grid_size), initialiser_function_zero(grid_size));
  
  //Grid A(n,  sin_boundary(n), initialiser_function_zero(n));
//   for(int i = 0; i< 20000; ++i)
//   {
//     A.rb_gauss_seidel_relaxation();
//   }

  Grids[1]->calculate_residual();
  Grids[1]->rb_gauss_seidel_relaxation();
  
  Grids[1]->print_v();
  Grids[1]->print_f();
  Grids[1]->print_r();
  Grids[2] = new Grid((grid_size+1)/2, Grids[1]->fw_restrict());
  //Grids[2]->print_v();

  return 0;

}

