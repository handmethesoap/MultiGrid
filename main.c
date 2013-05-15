#include <stdlib.h> 
#include <sys/time.h>
#include "grid.h"
#include "initialiserFunctions.h" 
#include "multigridSolver.h"

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
  
  //Create object for solving the problem
  MultigridSolver multigridProblem(l, sin_boundary((1 << l) + 1));
  
  //Solve the problem by performing n vcycles
  multigridProblem.solver(n);
  
  //print grid size
  std::cout << "h = 1/" << (1 << l) << std::endl;
  
  //Print computed solution to file and exact solution to file for comparison
  multigridProblem.print_solution();
  multigridProblem.print_exact_solution( sin_fxn((1 << l) + 1));
  
  //Print discretisation error to screen
  multigridProblem.print_discretisation_error( sin_fxn((1 << l) + 1));

  
  //calculate operation time and print it
  gettimeofday(&end, NULL);
  operation_time = (end.tv_sec - start.tv_sec)*1000000 + end.tv_usec - start.tv_usec; 
  std::cout << "Operating time (us): " << operation_time << std::endl;
  
  return 0;

}
