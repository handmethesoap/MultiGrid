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
      m_v = new double[n*n]();
      m_f = new double[n*n]();
      m_r = new double[n*n]();
      
      for( int i = 0; i < n*n; ++i)
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
      m_f = new double[n*n]();
      m_r = new double[n*n]();
      
      for( int i = 0; i < n*n; ++i)
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
      m_r = new double[n*n]();
      
      for( int i = 0; i < n*n; ++i)
      {
	m_r[i] = 0;
      }
  }
  
  Grid(int n, double *v_initial, double *f_initial, double *r_initial)
  {

      m_n = n;
      m_v = v_initial;
      m_f = f_initial;
      m_r = r_initial;
      
  }
  
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
  
  friend std::ostream& operator<< (std::ostream &out, Grid &outputGrid);
  
  double rb_gauss_seidel_relaxation(void);
  double jacobi_relaxation(void);
  double damped_jacobi_relaxation(double damping_factor);
  
  double* fw_restrict(void);
  double* bl_interpolate(void);
  void add_to_v( double * error_correction );
  void set_f( double * new_f ); 
  
  void calculate_residual(void);
  double calculate_L_inf_norm(double(*solution_function)(int, int, int));
  double calculate_L2_norm( double * exact_solution );
  double calculate_L2_norm( void );
  
  void print_v(void); 
  void print_f(void);
  void print_r(void);
  
  
};
  




#endif