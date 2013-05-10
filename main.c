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
  
  int n = 66;
  double initial_error;
  double current_error;
  int iterations;
  int k[9] = {1,8,16,24,32,40,48,56,64};
  
  for(int i = 0; i < 9; i++)
  {

    Grid A(n, boundary_function_zero, v_initialiser_function_sin_2, k[i]);
    initial_error = A.calculate_L_inf_norm(v_initialiser_function_zero);
    current_error = initial_error;
    iterations = 0;
    
    while( (current_error*100 > initial_error) && (iterations < 8000) )
    {
      A.damped_jacobi_relaxation(2/3);
      current_error = A.calculate_L_inf_norm(v_initialiser_function_zero);
      ++iterations;
    }
    std::cout << iterations << ", " << current_error/initial_error << std::endl;
  }
  //myfile.close();
  return 0;

}

