#include <iostream>
#include <fstream>
#include <stdlib.h> 
#include <memory>
#include "grid.h"
#include "initialiserFunctions.h"



int main(int argc, char* argv[])
{
 
  //int l = atoi(argv[1]);
  
  //int n = (1 << l) + 1;
  
  std::ofstream myfile;
  myfile.open("gridSizeVsError.dat");
  
  Grid* A;
  
  for( int n = 7; n <= 102; ++n)
  {
    A = new Grid(n, boundary_function_zero, v_initialiser_function_zero, f_initialiser_function_sin);
    for( int it_col = 0; it_col < 100; ++ it_col )
    {
      A->jacobi_relaxation();
    }
    myfile << n-2 << " " << A->calculate_L_inf_norm(exact_solution)<< std::endl;
    
    delete A;
  }
  
  myfile.close();
  return 0;
}

