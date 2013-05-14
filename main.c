#include <iostream>
#include <fstream>
#include <stdlib.h> 
#include <memory>
#include "grid.h"
#include "initialiserFunctions.h"



int main(int argc, char* argv[])
{
  
  
//   Grid A(9, initialiser_function_zero(16));
//   Grid B(17);
//   for(int i = 0; i < 2000; ++i)
//   {
//     B.add_to_v(A.bl_interpolate());
//   }
//   
//   B.print_v();
  
  int l = atoi(argv[1]);
  int n = atoi(argv[2]);
  
  int grid_size = (1 << l) + 1;
  
  double L2_norm;
  
  Grid **Grids = new Grid*[l];
  
  Grids[0] = new Grid(grid_size,  sin_boundary(grid_size));
  
  for(int it = 1; it < l; ++it)
  {
      Grids[it] = new Grid((1 << (l-it)) + 1);
      
  }

  
  for(int vcycles = 0; vcycles < n; ++vcycles)
  {
    for(int it = 0; it < l -1; ++it)
    {
      Grids[it]->rb_gauss_seidel_relaxation();
      Grids[it]->rb_gauss_seidel_relaxation();
      Grids[it]->calculate_residual();
      
      Grids[it+1]->set_f(Grids[0]->fw_restrict());// = new Grid((1 << it) + 1, initialiser_function_zero((1 << it) + 1), Grids[0]->fw_restrict());
    }

    Grids[l-1]->rb_gauss_seidel_relaxation();

    for(int it = l - 1; it > 0; --it)
    {
      Grids[it-1]->add_to_v( Grids[it]->bl_interpolate() );
      Grids[it-1]->rb_gauss_seidel_relaxation();
    }
  }
  
  std::ofstream myfile;
  myfile.open("solution.dat");
  
  myfile << *Grids[0] << std::endl;
  
  myfile.close();
  
  L2_norm = Grids[0]->calculate_L2_norm(sin_fxn(grid_size));
  
  std::cout << L2_norm << std::endl;
  
  //Grids[0]->print_v();
  
  for( int it = 0; it < l; ++it)
  {
    delete Grids[it];
  }
  
  delete[] Grids;
  
  return 0;

}

