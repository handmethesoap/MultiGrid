#include <iostream>
#include <fstream>
#include <stdlib.h> 
#include <memory>
#include <sys/time.h>
#include "grid.h"
#include "initialiserFunctions.h"



int main(int argc, char* argv[])
{
  
  timeval start;
  timeval end;
  gettimeofday(&start, NULL);
  
  int l = atoi(argv[1]);
  int n = atoi(argv[2]);
  
  int grid_size = (1 << l) + 1;
  int operation_time = 0;
  
  double L2_norm = 0;
  double previous_L2_norm = 0;
  double convergence_rate = 0;
  
  Grid **Grids = new Grid*[l];
  
  Grids[0] = new Grid(grid_size,  sin_boundary(grid_size));
  
  previous_L2_norm = Grids[0]->calculate_L2_norm(sin_fxn(grid_size));
  
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
      
      Grids[it+1]->set_f(Grids[0]->fw_restrict());
    }

    Grids[l-1]->rb_gauss_seidel_relaxation();

    for(int it = l - 1; it > 0; --it)
    {
      Grids[it-1]->add_to_v( Grids[it]->bl_interpolate() );
      Grids[it-1]->rb_gauss_seidel_relaxation();
    }
//     L2_norm = Grids[0]->calculate_L2_norm(sin_fxn(grid_size));
//     convergence_rate = L2_norm/previous_L2_norm;
    
//     std::cout << "L2 norm = "<< L2_norm << std::endl;
//     std::cout << "Convergence Rate = "<< convergence_rate << std::endl << std::endl;
//     previous_L2_norm = L2_norm;
  }
  std::cout << "h = 1/" << grid_size-1 << std::endl;
  std::ofstream myfile;
  myfile.open("solution.txt");
  
  myfile << *Grids[0] << std::endl;
  
  myfile.close();
  
   L2_norm = Grids[0]->calculate_L2_norm(sin_fxn(grid_size));
  std::cout << "L2 norm = "<< L2_norm << std::endl;
  
  //Grids[0]->print_v();
  
  for( int it = 0; it < l; ++it)
  {
    delete Grids[it];
  }
  
  delete[] Grids;
  
  gettimeofday(&end, NULL);
  operation_time = (end.tv_sec - start.tv_sec)*1000000 + end.tv_usec - start.tv_usec; 
  std::cout << "Operating time (us): " << operation_time << std::endl;
  
  return 0;

}

