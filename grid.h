#ifndef GRID_H
#define GRID_H

#include <iostream>

class Grid
{
  
private:
  
  double* m_v;
  double* m_f;
  int m_n;
  
public:
  
  //default constructor
  Grid(int n)
  {
      m_n = n;
      m_v = new double[n*n];
      m_f = new double[n*n];
      
      for( int i = 0; i < n*n; i ++)
      {
	m_v[i] = 0;
	m_f[i] = 0;
      }

  }
  
  Grid(int n, double(*boundary_function)(int, int, int))
  {
      m_n = n;
      m_v = new double[n*n];
      m_f = new double[n*n];
      
      for( int i = 0; i < n*n; i ++)
      {
	m_v[i] = 0;
	m_f[i] = 0;
      }
      
      set_boundary(boundary_function);
  }
  
  Grid(int n, double(*boundary_function)(int, int, int), double(*v_initialiser_function)(int, int, int))
  {
      m_n = n;
      m_v = new double[n*n];
      m_f = new double[n*n];
      
      for( int i = 0; i < n*n; i ++)
      {
	m_f[i] = 0;
      }
      
      set_initial(v_initialiser_function);
      set_boundary(boundary_function);
  }
  
  Grid(int n, double(*boundary_function)(int, int, int), double(*v_initialiser_function)(int, int, int), double(*f_initialiser_function)(int, int, int))
  {
      m_n = n;
      m_v = new double[n*n];
      m_f = new double[n*n];
      
      set_initial(v_initialiser_function, f_initialiser_function);
      set_boundary(boundary_function);
  }
  
  ~Grid(void)
  {
    delete[] m_v;
    delete[] m_f;
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
  void damped_jacobi_relaxation(int damping_factor);
  void full_weight_restriction(void);
  double calculate_L_inf_norm(double(*solution_function)(int, int, int));
  void print_v(void); // overload ostream operator instead
  void print_f(void);
  
  
};
  
#endif