#ifndef GRID_H
#define GRID_H

#include <iostream>

class Grid
{
  
private:
  
  double* m_v;
  double* m_f;
  double* m_r;
  int m_n;
  
  
  
public:
  
  //default constructor
  Grid(int n)
  {
      m_n = n;
      m_v = new double[n*n];
      m_f = new double[n*n];
      m_r = new double[n*n];
      
      for( int i = 0; i < n*n; i ++)
      {
	m_v[i] = 0;
	m_f[i] = 0;
	m_r[i] = 0;
      }

  }
  
  Grid(int n, double *v_initial)
  {
      m_n = n;
      m_v = v_initial;
      m_f = new double[n*n];
      m_r = new double[n*n];
      
      for( int i = 0; i < n*n; i ++)
      {
	m_f[i] = 0;
	m_r[i] = 0;
      }
  }
  
  Grid(int n, double *v_initial, double *f_initial)
  {
      m_n = n;
      m_v = v_initial;
      m_f = f_initial;
      m_r = new double[n*n];
      
      for( int i = 0; i < n*n; i ++)
      {
	m_r[i] = 0;
      }
  }
  
  // Grid(int n, double(*v_initialiser_function)(int, int, int))
  // {
      // m_n = n;
      // m_v = new double[n*n];
      // m_f = new double[n*n];
      // m_r = new double[n*n];
      
      // for( int i = 0; i < n*n; i ++)
      // {
	// m_v[i] = 0;
	// m_f[i] = 0;
	// m_r[i] = 0;
      // }
      
      // set_initial(v_initialiser_function);
  // }
  
  
  // Grid(int n,  double(*v_initialiser_function)(int, int, int), double(*boundary_function)(int, int, int))
  // {
      // m_n = n;
      // m_v = new double[n*n];
      // m_f = new double[n*n];
      // m_r = new double[n*n];
      
      // for( int i = 0; i < n*n; i ++)
      // {
	// m_f[i] = 0;
	// m_r[i] = 0;
      // }
      
      // set_initial(v_initialiser_function);
      // set_boundary(boundary_function);
  // }
   // Grid(int n, double(*v_initialiser_function)(int, int, int), double(*boundary_function)(int, int, int),  double(*f_initialiser_function)(int, int, int))
  // {
      // m_n = n;
      // m_v = new double[n*n];
      // m_f = new double[n*n];
      // m_r = new double[n*n];
      
      // for( int i = 0; i < n*n; i ++)
      // {
	// m_r[i] = 0;
      // }
      
      // set_initial(v_initialiser_function, f_initialiser_function);
      // set_boundary(boundary_function);
  // }
  
//   Grid(int n, double(*boundary_function)(int, int, int), double(*v_initialiser_function)(int, int, int), double(*f_initialiser_function)(int, int, int))
//   {
//       m_n = n;
//       m_v = new double[n*n];
//       m_f = new double[n*n];
//       
//       set_initial(v_initialiser_function, f_initialiser_function);
//       set_boundary(boundary_function);
//   }
//   
  ~Grid(void)
  {
    delete[] m_v;
    delete[] m_f;
    delete[] m_r;
  }
  
  double* get_v(void)
  {
    return m_v;
  }
  
  double* get_f(void)
  {
    return m_f;
  }
  
  int get_dimension(void)
  {
    return m_n;
  }
  
  
  
  void set_boundary(double(*boundary_function)(int, int, int));
  void set_initial(double(*v_initialiser_function)(int, int, int));
  void set_initial(double(*v_initialiser_function)(int, int, int), double(*f_initialiser_function)(int, int, int));
  
  void rb_gauss_seidel_relaxation(void);
  double jacobi_relaxation(void);
  double damped_jacobi_relaxation(double damping_factor);
  
  //Grid* full_weight_restriction(void);
  double* fw_restrict(void);
  
  double calculate_L_inf_norm(double(*solution_function)(int, int, int));
  void calculate_residual(void);
  
  void print_v(void); // overload ostream operator instead
  void print_f(void);
  void print_r(void);
  
  
};
  
#endif