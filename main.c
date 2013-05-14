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
  
  Grids[0] = new Grid(grid_size,  sin_boundary(grid_size));
  
  //Grid A(n,  sin_boundary(n), initialiser_function_zero(n));
//   for(int i = 0; i< 20000; ++i)
//   {
//     A.rb_gauss_seidel_relaxation();
//   }

  //Grids[1]->calculate_residual();
  //Grids[1]->rb_gauss_seidel_relaxation();
  
  Grids[0]->print_v();
  Grids[0]->print_f();
  Grids[0]->print_r();
  
  Grids[0]->rb_gauss_seidel_relaxation();
  Grids[0]->rb_gauss_seidel_relaxation();
  Grids[0]->calculate_residual();
  
  Grids[0]->print_v();
  Grids[0]->print_f();
  Grids[0]->print_r();
  
  Grids[1] = new Grid((grid_size+1)/2, initialiser_function_zero((grid_size+1)/2), Grids[0]->fw_restrict());
  Grids[1]->rb_gauss_seidel_relaxation();
  Grids[1]->rb_gauss_seidel_relaxation();
  Grids[1]->calculate_residual();
  
  Grids[1]->print_v();
  Grids[1]->print_f();
  Grids[1]->print_r();
  
  Grids[2] = new Grid(((grid_size+1)/2+1)/2, initialiser_function_zero(((grid_size+1)/2+1)/2), Grids[1]->fw_restrict());
  Grids[2]->rb_gauss_seidel_relaxation();
  Grids[2]->rb_gauss_seidel_relaxation();
  
  Grids[1]->add_to_v( Grids[2]->bl_interpolate());
  Grids[1]->rb_gauss_seidel_relaxation();
  
  Grids[0]->add_to_v( Grids[1]->bl_interpolate());
  Grids[0]->rb_gauss_seidel_relaxation();
  
  Grids[0]->print_v();
  Grids[0]->print_f();
  Grids[0]->print_r();
  
  
  
  //Grids[2] = new Grid(grid_size, Grids[1]->bl_interpolate());
  
  
  delete Grids[0];
  delete Grids[1];
  delete Grids[2];
  std::cout << "here" << std::endl;
  delete[] Grids;
  
  return 0;

}

