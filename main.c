#include <iostream>
#include <fstream>
#include <stdlib.h> 
#include <memory>
#include <sys/time.h>
#include "grid.h"
#include "initialiserFunctions.h"

void print_exact_solution( int grid_size, double *solution_function);
void print_solution(Grid* solution_grid);
void multigrid_solver(Grid** Grids, int number_vcycles, int levels);

int main(int argc, char* argv[])
{
  //initialise timing variables
  timeval start;
  timeval end;
  gettimeofday(&start, NULL);
  
  //read grid size and number of v cycles
  int l = atoi(argv[1]);
  int n = atoi(argv[2]);
  
  //Initialse variables
  int operation_time = 0;
  
  //Declare array of pointers to grid objects for each level
  Grid** Grids = new Grid*[l];
  
  multigrid_solver(Grids, n, l);
  
  std::cout << "h = 1/" << (1 << l) << std::endl;
  
  //Print computed solution to file and exact solution to file for comparison
  print_solution(Grids[0]);
  print_exact_solution( (1 << l) + 1, sin_fxn((1 << l) + 1));
  
  //delete data structures
  for( int it = 0; it < l; ++it)
  {
    delete Grids[it];
  }
  
  delete[] Grids;
  
  //calculate operation time and print it
  gettimeofday(&end, NULL);
  operation_time = (end.tv_sec - start.tv_sec)*1000000 + end.tv_usec - start.tv_usec; 
  std::cout << "Operating time (us): " << operation_time << std::endl;
  
  return 0;

}


void multigrid_solver(Grid** Grids, int number_vcycles, int levels)
{
  
  //Calculate grid size
  int grid_size = (1 << levels) + 1;
  
  //initialise variables
  double L2_norm = 0;
  double previous_L2_norm = 0;
  double convergence_rate = 0;
  
  //Initialse coarsest grid with known variables 
  Grids[0] = new Grid(grid_size,  sin_boundary(grid_size));
  
  //Calculate the initial error norm with the initial guess
  previous_L2_norm = Grids[0]->calculate_L2_norm(sin_fxn(grid_size));
  
  //Initialse the grids for other levels
  for(int it = 1; it < levels; ++it)
  {
      Grids[it] = new Grid((1 << (levels-it)) + 1);
  }

  //Enter loop executing a single v cycle each pass
  for(int vcycles = 0; vcycles < number_vcycles; ++vcycles)
  {
    //perform two gauss seidel relaxations on each level then calculate the residual and pass it to the next coarset grid
    for(int it = 0; it < levels -1; ++it)
    {
      Grids[it]->rb_gauss_seidel_relaxation();
      Grids[it]->rb_gauss_seidel_relaxation();
      Grids[it]->calculate_residual();
      
      Grids[it+1]->set_f(Grids[0]->fw_restrict());
    }
    
    //Perfom a single gauss seidel relaxation on the coarsest grid in order to find its exact solution
    Grids[levels-1]->rb_gauss_seidel_relaxation();
    
    //Add the calculated error from each grid to the next finest grid and then perform a single gauss seidel relaxation on that grid
    for(int it = levels - 1; it > 0; --it)
    {
      Grids[it-1]->add_to_v( Grids[it]->bl_interpolate() );
      Grids[it-1]->rb_gauss_seidel_relaxation();
    }
    
    //Calculate the L2 norm and the convergence rate for each v cycle and print them
    L2_norm = Grids[0]->calculate_L2_norm(sin_fxn(grid_size));
    convergence_rate = L2_norm/previous_L2_norm;
    
    std::cout << "L2 norm = "<< L2_norm << std::endl;
    std::cout << "Convergence Rate = "<< convergence_rate << std::endl << std::endl;
    previous_L2_norm = L2_norm;
  }
}

void print_solution(Grid* solution_grid)
{
  std::ofstream myfile;
  myfile.open("solution.txt");
  
  myfile << *solution_grid << std::endl;
  
  myfile.close();
}

void print_exact_solution( int grid_size, double *solution_function)
{
  Grid* exact_solution = new Grid(grid_size, solution_function);
  std::ofstream myfile;
  
  myfile.open("exact_solution.txt");
  myfile << *exact_solution << std::endl;
  myfile.close();
  
}
