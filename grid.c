#include "grid.h"
#include <cmath>

void Grid:: set_boundary(double(*boundary_function)(int, int, int))
{
  for( int it = 0; it < m_n; ++it)
  {
    m_v[it] = boundary_function(0, it, m_n);
    m_v[m_n*m_n - m_n + it] = boundary_function(m_n, it, m_n);
    m_v[it*m_n] = boundary_function(it, 0, m_n);
    m_v[it*m_n + m_n - 1] = boundary_function(it, m_n - 1, m_n);
    
    m_f[it] = boundary_function(0, it, m_n);
    m_f[m_n*m_n - m_n + it] = boundary_function(m_n, it, m_n);
    m_f[it*m_n] = boundary_function(it, 0, m_n);
    m_f[it*m_n + m_n - 1] = boundary_function(it, m_n - 1, m_n);
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

void Grid:: rb_gauss_seidel_relaxation(void)
{
  int start;
  int r1, r2, r3, r4, r5;
  int h2 = 1.0/(m_n-1)*(m_n-1);
  
  for( int it1 = 1; it1 < ((m_n)-1); ++it1 )
  {
    start = it1%2 +1;
    r1 = it1 * m_n;
    r2 = r1 + 1;
    r3 = r1 - 1;
    r4 = (it1 - 1) * m_n;
    r5 = (it1 +1) * m_n;
    
    for( int it = start; it < ((m_n)-1); it += 2 )
    {
      m_v[ r1 + it] = 0.25*(m_v[ r2 + it ] + m_v[ r3 + it ] + m_v[ r4 + it] + m_v[ r5 + it] + m_f[r1 + it]*h2);
    }
  }
  
  for( int it1 = 1; it1 < ((m_n)-1); ++it1 )
  {
    start = (it1+1)%2 +1;
    r1 = it1 * m_n;
    r2 = r1 + 1;
    r3 = r1 - 1;
    r4 = (it1 - 1) * m_n;
    r5 = (it1 +1) * m_n;
    
    for( int it = start; it < ((m_n)-1); it += 2 )
    {
      m_v[ r1 + it] = 0.25*(m_v[ r2 + it ] + m_v[ r3 + it ] + m_v[ r4 + it] + m_v[ r5 + it] + m_f[r1 + it]*h2);
    }
  }
  
}

double Grid:: jacobi_relaxation(void)
{
  int row_offset, r1, r2, r3, r4;
  double h2 = 1.0/((m_n-1)*(m_n-1));
  double* m_v_temp = new double[m_n*m_n];
  double residual = 0.0;
  
  //copy boundary
  for( int it = 0; it < m_n; ++it)
  {
    m_v_temp[it] = m_v[it];
    m_v_temp[m_n*m_n - m_n + it] = m_v[m_n*m_n - m_n + it];
    m_v_temp[it*m_n] = m_v[it*m_n];
    m_v_temp[it*m_n + m_n - 1] = m_v[it*m_n + m_n - 1];
  }
  
  for( int it_row = 1; it_row < ((m_n)-1); ++it_row )
 {
    row_offset = it_row * m_n;
    r1 = row_offset - m_n;
    r2 = row_offset + m_n;
    r3 = row_offset - 1;
    r4 = row_offset + 1;
    for( int it_col = 1; it_col < ((m_n)-1); ++ it_col )
    {
      m_v_temp[row_offset + it_col] = 0.25*(m_v[r1 + it_col] + m_v[r2 + it_col] + m_v[r4 + it_col] + m_v[r3 + it_col] + m_f[row_offset + it_col]*h2);
      residual += std::abs(m_v_temp[row_offset + it_col] - m_v[row_offset + it_col]);
    }
 }
 delete[] m_v;
 m_v = m_v_temp;
 return residual*h2;
}

double Grid:: damped_jacobi_relaxation(double damping_factor)
{
  int row_offset, r1, r2, r3, r4;
  double h2 = 1.0/((m_n-1)*(m_n-1));
  double* m_v_temp = new double[m_n*m_n];
  double residual = 0.0;
  
  //copy boundary
  for( int it = 0; it < m_n; ++it)
  {
    m_v_temp[it] = m_v[it];
    m_v_temp[m_n*m_n - m_n + it] = m_v[m_n*m_n - m_n + it];
    m_v_temp[it*m_n] = m_v[it*m_n];
    m_v_temp[it*m_n + m_n - 1] = m_v[it*m_n + m_n - 1];
  }
  
  for( int it_row = 1; it_row < ((m_n)-1); ++it_row )
 {
    row_offset = it_row * m_n;
    r1 = row_offset - m_n;
    r2 = row_offset + m_n;
    r3 = row_offset - 1;
    r4 = row_offset + 1;
    for( int it_col = 1; it_col < ((m_n)-1); ++ it_col )
    {
      m_v_temp[row_offset + it_col] = (1.0 - damping_factor)*m_v[row_offset + it_col] + damping_factor*0.25*(m_v[r1 + it_col] + m_v[r2 + it_col] + m_v[r4 + it_col] + m_v[r3 + it_col] + m_f[row_offset + it_col]*h2);
      residual += std::abs(m_v_temp[row_offset + it_col] - m_v[row_offset + it_col]);
      //std::cout << residual << ", " << m_v_temp[row_offset + it_col] << ", " << m_v[row_offset + it_col] << std::endl;
    }
 }
 delete[] m_v;
 m_v = m_v_temp;
 return residual*h2;
}

double Grid:: calculate_L_inf_norm(double(*solution_function)(int, int, int))
{
  double max = 0.0;
  double temp;
  for( int it_row = 0; it_row < (m_n); ++it_row )
  {
    for( int it_col = 0; it_col < (m_n); ++ it_col )
    {
      temp = std::abs(temp = m_v[it_row * m_n + it_col] - solution_function(it_row, it_col, m_n));
      if(temp > max)
      {
	max = temp;
      }
    }
  }
  return max;
}

void Grid:: print_v(void)
{
  std::cout << "v=" << std::endl;
  
  for( int i = 0; i < m_n; ++i)
  {
    for( int j = 0; j < m_n; ++j )
    {
      std::cout << m_v[i*(m_n) + j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

  
void Grid:: print_f(void)
{
  std::cout << "f=" << std::endl;
  
  for( int i = 0; i < m_n; ++i)
  {
    for( int j = 0; j < m_n; ++j )
    {
      std::cout << m_f[i*(m_n) + j] << " ";
    }
    std::cout << std::endl;
  }
}

void Grid:: print_r(void)
{
   std::cout << "r=" << std::endl;
  
  for( int i = 0; i < m_n; ++i)
  {
    for( int j = 0; j < m_n; ++j )
    {
      std::cout << m_r[i*(m_n) + j] << " ";
    }
    std::cout << std::endl;
  }
}

void Grid:: calculate_residual(void)
{
  
  int row_offset, r1, r2, r3, r4;
  double h2 = 1.0/((m_n-1)*(m_n-1));
  
  //set error on boundary to zero
  for( int it = 0; it < m_n; ++it)
  {
    m_r[it] = 0;
    m_r[m_n*m_n - m_n + it] = 0;
    m_r[it*m_n] = 0;
    m_r[it*m_n + m_n - 1] = 0;
  }
  
  for( int it_row = 1; it_row < ((m_n)-1); ++it_row )
 {
    row_offset = it_row * m_n;
    r1 = row_offset - m_n;
    r2 = row_offset + m_n;
    r3 = row_offset - 1;
    r4 = row_offset + 1;
    for( int it_col = 1; it_col < ((m_n)-1); ++ it_col )
    {
      m_r[row_offset + it_col] = h2*(m_v[r1 + it_col] + m_v[r2 + it_col] + m_v[r4 + it_col] + m_v[r3 + it_col] - 4*m_v[row_offset + it_col]) + m_f[row_offset + it_col];
    }
 }

}



double* Grid:: fw_restrict(void)
{
  double *restricted_v = new double[((m_n+1)/2)*((m_n+1)/2)];
  
  //Set boundaries to zero
  for( int it = 0; it < (m_n+1)/2; ++it)
  {
    restricted_v[it] = 0;
    restricted_v[((m_n+1)/2)*((m_n+1)/2) - ((m_n+1)/2) + it] = 0;
    restricted_v[it*((m_n+1)/2)] = 0;
    restricted_v[it*((m_n+1)/2) + ((m_n+1)/2) - 1] = 0;
  }
  
  for( int it_row = 1; it_row < ((m_n-1)/2); ++it_row )
  {
    
    for( int it_col = 1; it_col < ((m_n-1)/2); ++ it_col )
    {
      
      restricted_v[it_row*(m_n+1)/2 + it_col] = 0.0625*( m_r[(2*it_row -1) * m_n + (2*it_col - 1)] + m_r[(2*it_row -1) * m_n + (2*it_col + 1)] + m_r[(2*it_row + 1) * m_n + (2*it_col - 1)] + m_r[(2*it_row + 1) * m_n + (2*it_col + 1)] + 2*( m_r[(2*it_row + 1) * m_n + (2*it_col)] + m_r[(2*it_row - 1) * m_n + (2*it_col)] + m_r[(2*it_row) * m_n + (2*it_col - 1)] + m_r[(2*it_row) * m_n + (2*it_col + 1)]) + 4*m_r[(2*it_row) * m_n + (2*it_col)]);
      //std::cout << "test " << it_row*(m_n+1)/2 + it_col<< std::endl;
    }
    
  }
  
  return restricted_v;
  
}

double* Grid:: bl_interpolate(void)
{
  double* interpolated_v = new double[(m_n*2-1)*(m_n*2-1)];
  
  //Set boundaries to zero
  for( int it = 0; it < (m_n*2-1); ++it)
  {
    interpolated_v[it] = 0;
    
    interpolated_v[(m_n*2-1)*(m_n*2-1) - (m_n*2-1) + it] = 0;
    
    interpolated_v[it*(m_n*2-1)] = 0;
    
    interpolated_v[it*(m_n*2-1) + (m_n*2-1) - 1] = 0;
  }
  
  // set all internal entries to zero temporarily
  for( int it_row = 1; it_row < (m_n*2-1) -1; ++it_row )
  {
    for(int it_col = 1; it_col < (m_n*2-1) -1; ++it_col )
    {
      interpolated_v[it_row*(m_n*2-1) + it_col] = 0;
    }
  }
  
  //copy directly interior values of coarser grid to appropriate locations on finer grid
  for( int it_row = 1; it_row < m_n - 1; ++it_row )
  {
    
    for( int it_col = 1; it_col < m_n - 1; ++ it_col )
    {
      
      interpolated_v[(2*it_row) * (m_n*2-1) + (2*it_col)] = m_v[it_row*m_n + it_col];
    }
    
  }
  
  //interpolate points with only diagonal connections
  for( int it_row = 1; it_row < (m_n*2-1); it_row += 2 )
  {
    
    for( int it_col = 1; it_col < (m_n*2-1) ; it_col += 2 )
    {
      
      interpolated_v[(it_row) * (m_n*2-1) + (it_col)] = 0.25*(interpolated_v[(it_row -1) * (m_n*2-1) + (it_col - 1)] + interpolated_v[(it_row -1) * (m_n*2-1) + (it_col + 1)] + interpolated_v[(it_row + 1) * (m_n*2-1) + (it_col - 1)] + interpolated_v[(it_row +1) * (m_n*2-1) + (it_col + 1)]);
    }
    
  }
  //interpolate points with values above and below
  for( int it_row = 1; it_row < (m_n*2-1); it_row += 2 )
  {
    
    for( int it_col = 2; it_col < (m_n*2-1) ; it_col += 2 )
    {
      
      interpolated_v[(it_row) * (m_n*2-1) + (it_col)] = 0.5*(interpolated_v[(it_row + 1) * (m_n*2-1) + (it_col)] + interpolated_v[(it_row - 1) * (m_n*2-1) + (it_col)]);
      
    }
    
  }
  
  //interpolate points with values left and right
  for( int it_row = 2; it_row < (m_n*2-1); it_row += 2 )
  {
    
    for( int it_col = 1; it_col < (m_n*2-1) ; it_col += 2 )
    {
      
      interpolated_v[(it_row) * (m_n*2-1) + (it_col)] = 0.5*(interpolated_v[(it_row) * (m_n*2-1) + (it_col + 1)] + interpolated_v[(it_row) * (m_n*2-1) + (it_col - 1)]);
      
    }
    
  }
  return interpolated_v;
}


void Grid:: add_to_v( double * error_correction )
{
  for( int i = 0; i < m_n; ++i)
  {
    for( int j = 0; j < m_n; ++j )
    {
      m_v[i*(m_n) + j] += error_correction[i*(m_n) + j];
    }
  }
  delete error_correction;
}

double Grid:: calculate_L2_norm( double * exact_solution )
{
  double L2 = 0.0;
  
  for( int i = 0; i < m_n; ++i)
  {
    for( int j = 0; j < m_n; ++j )
    {
      L2 += sqrt((m_v[i*(m_n) + j] - exact_solution[i*(m_n) + j])*(m_v[i*(m_n) + j] - exact_solution[i*(m_n) + j]));
      std::cout <<  m_v[i*(m_n) + j] << "-" << exact_solution[i*(m_n) + j] << std::endl;
    }

  }
  
  delete[] exact_solution;
  return L2;
}
  
void Grid:: set_f( double * new_f )
{
  delete[] m_f;
  
  m_f = new_f;
}













