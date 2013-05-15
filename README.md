MultiGrid
=========
To compile mgsolve call:
make

To operate mgsolve call:
./mgsolve l n
Where l is the number of levels for solving the multigrid problems and n is the number of v cycles to perform

To solve a different multigrid problem edit main.  A MultigridSolver object can be created and the constructor can take initialisation 
arrays for the v,f and residual arrays.  Examples of functions generating these arrays can be found in initialiserFunctions.c.  
The solution of the multigrid problem is calculated by calling solver(n) where n is the number of v cycles to perform.
The solution may be printed by a call to print_solution(); and is printed to a file solution.txt in a format compatible with gnuplot.
A similar call to print_exact_solution( double* exact_solution ) prints the exact solution to a file (if it is known) in a format 
compatiable with gnuplot. print_discretisation_error( double* exact solution ); can print the discretisation error of the solution
to the terminal if the exact solution is known.