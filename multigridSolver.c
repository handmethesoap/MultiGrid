#include "multigridSolver.h"
#include <iostream>
#include <fstream>

MultigridSolver:: MultigridSolver(int l)
{
  Grids = new Grid*[l];
  grid_size = (1 << levels) + 1;
  levels = l;
}

MultigridSolver:: MultigridSolver(int l, double *v_initial)
{
  Grids = new Grid*[l];
  grid_size = (1 << l) + 1;
  Grids[0] = new Grid(grid_size,  v_initial);
  levels = l;
}

MultigridSolver:: MultigridSolver(int l, double *v_initial, double *f_initial)
{
  Grids = new Grid*[l];
  grid_size = (1 << l) + 1;
  Grids[0] = new Grid(grid_size,  v_initial, f_initial);
  levels = l;
}

MultigridSolver:: MultigridSolver(int l, double *v_initial, double *f_initial, double *r_initial)
{
  Grids = new Grid*[l];
  grid_size = (1 << l) + 1;
  Grids[0] = new Grid(grid_size,  v_initial, f_initial, r_initial);
  levels = l;
}

MultigridSolver:: ~MultigridSolver()
{
  for( int it = 0; it < levels; ++it)
  {
    delete Grids[it];
  }
  
  delete[] Grids;
}

void MultigridSolver:: solver(int number_vcycles)
{
  
  //initialise variables
  double L2_norm = 0;
  double previous_L2_norm = 1;
  double convergence_rate = 0;
  
  
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
      Grids[it+1]->set_f(Grids[it]->fw_restrict());

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
    L2_norm = Grids[0]->calculate_L2_norm();
    convergence_rate = L2_norm/previous_L2_norm;
    
    std::cout << "L2 norm = "<< L2_norm << std::endl;
    std::cout << "Convergence Rate = "<< convergence_rate << std::endl << std::endl;
    previous_L2_norm = L2_norm;

  }
}

void MultigridSolver:: print_discretisation_error(double *solution_function)
{
  double L2_norm = Grids[0]->calculate_L2_norm(solution_function);    
  std::cout << "discretisation error = "<< L2_norm << std::endl;
}

void MultigridSolver:: print_solution(int level)
{
  std::ofstream myfile;
  myfile.open("solution.txt");
  
  myfile << *Grids[level] << std::endl;
  
  myfile.close();
}

void MultigridSolver:: print_exact_solution( double *solution_function )
{
  Grid* exact_solution = new Grid(grid_size, solution_function);
  std::ofstream myfile;
  
  myfile.open("exact_solution.txt");
  myfile << *exact_solution << std::endl;
  myfile.close();

  delete exact_solution;
  
}