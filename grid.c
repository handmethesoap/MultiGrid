#include "grid.h"

void Grid:: set_boundary(double(*boundary_function)(int, int, int))
{
  for( int it = 0; it < m_n; ++it)
  {
    m_v[it] = boundary_function(0, it, m_n);
    m_v[m_n*m_n - m_n + it] = boundary_function(m_n, it, m_n);
    m_v[it*m_n] = boundary_function(it, 0, m_n);
    m_v[it*m_n + m_n - 1] = boundary_function(it, m_n, m_n);
  }
  
}

void Grid:: set_initial(double(*v_initialiser_function)(int, int, int))
{
   for( int it_row = 0; it_row < (m_n); ++it_row )
   {
     for( int it_col = 0; it_col < (m_n); ++ it_col )
     {
       m_v[it_row * m_n + it_col] = v_initialiser_function(it_row, it_col, m_n);
     }
   }
}

void Grid:: set_initial(double(*v_initialiser_function)(int, int, int), double(*f_initialiser_function)(int, int, int))
{
   for( int it_row = 0; it_row < (m_n); ++it_row )
   {
     for( int it_col = 0; it_col < (m_n); ++ it_col )
     {
       m_v[it_row * m_n + it_col] = v_initialiser_function(it_row, it_col, m_n);
       m_f[it_row * m_n + it_col] = f_initialiser_function(it_row, it_col, m_n);
     }
   }
}

void Grid:: print_v(void)
{
  std::cout << "v=" << std::endl;
  
  for( int i = 0; i < m_n; ++i)
  {
    for( int j = 0; j < m_n; ++j )
    {
      std::cout << m_v[i*m_n + j] << " ";
    }
    std::cout << std::endl;
  }
}

  
void Grid:: print_f(void)
{
  std::cout << "f=" << std::endl;
  
  for( int i = 0; i < m_n; ++i)
  {
    for( int j = 0; j < m_n; ++j )
    {
      std::cout << m_f[i*m_n + j] << " ";
    }
    std::cout << std::endl;
  }
}

