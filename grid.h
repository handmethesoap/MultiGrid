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
      std::cout << "test" << std::endl;

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
  
  void set_boundary(double(*boundary_function)(int, int));
  void set_initial(double(*v_initialiser_function)(int, int), double(*f_initialiser_function)(int, int));
  void print_v(void); // overload ostream operator instead
  void print_f(void);
  
  
};
  
#endif