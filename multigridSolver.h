#ifndef GRIDSOLVER_H
#define GRIDSOLVER_H

#include "grid.h"
#include <cmath>

class MultigridSolver
{
  private:
  
    Grid** Grids;
    int grid_size;
    int levels;
  
  public:
    
    //constructor
    MultigridSolver(int l);
    
    MultigridSolver(int l, double *v_initial);
    
    MultigridSolver(int l, double *v_initial, double *f_initial);
    
    MultigridSolver(int l, double *v_initial, double *f_initial, double *r_initial);
    
    ~MultigridSolver();
    
    void solver(int number_vcycles);
    void print_solution(int level);
    void print_exact_solution(double *solution_function);
    void print_discretisation_error(double *solution_function);
};

































#endif