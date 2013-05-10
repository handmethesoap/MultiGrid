#include <iostream>
#include <fstream>
#include <stdlib.h> 
#include <memory>
#include "grid.h"
#include "initialiserFunctions.h"



int main(int argc, char* argv[])
{
 
//   int l = atoi(argv[1]);
//   
//   int n = (1 << l) + 1;
//   
//   Grid A(n, boundary_function_zero, v_initialiser_function_sin_2, );
//   for(int i = 0; i< 50; ++i)
//   {
//     A.rb_gauss_seidel_relaxation();
//   }
//   A.print_v();
  //std::ofstream myfile;
  //myfile.open("");
  
  Grid* A;
  double initial_error;
  double current_error;
  
  for( int n = 7; n <= 102; ++n)
  {
    A = new Grid(66, boundary_function_zero, v_initialiser_function_sin_2);
    initial_error = A->calculate_L_inf_norm(v_initialiser_function_zero);
    current_error = initial_error;
    
    while( 100*current_error > initial_error )
    {
      A->jacobi_relaxation();
    }
    //myfile << n-2 << " " << A->calculate_L_inf_norm(exact_solution)<< std::endl;
    
    delete A;
  }
  
  //myfile.close();
  return 0;

}

