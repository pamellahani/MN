#include <stdio.h>
#include <omp.h>
#include "bib3.h"

void init_vector (VecFloat V, float valeur)
{
  int i ;

#pragma parallel for private (i)
  for (i = 0; i < N; i++)
    V [i] = valeur ;
}

void copy (VecFloat V1, VecFloat V2)
{
  int i ;
#pragma omp parallel for private (i)
    for (i = 0 ; i < N ; i+=4)
      {
	V2 [i] = V1 [i] ;
	V2 [i+1] = V1 [i+1] ;
	V2 [i+2] = V1 [i+2] ;
	V2 [i+3] = V1 [i+3] ;
      }
}

float dot (VecFloat V1, VecFloat V2)
{
  int i ;
  float d = 0.0 ;
#pragma parallel omp for private (i) reduction (+:d)  
  for (i = 0 ; i < N ; ++i)
     d += (V2 [i] * V1 [i]) ;
  return d ;
}

void init_matrix (MatFloat m, float v)
{
  int i, j;

  for (i = 0; i < M; ++i)
    for (j = 0 ; j < M; ++j)
      m [i][j] = v ;
}

void mult_mat_vect (MatFloat m, VectMatFloat v, VectMatFloat vr)
{
  int i, j ;
  float d ;

#pragma omp parallel for private (d,i,j)  
  for (i = 0; i < M; ++i)
    {
      d = 0.0 ;
      for (j = 0; j < M; ++j)
	d += m [i][j] * v [j] ;
      vr [i] = d ;
    }
}

void mult_mat_mat (MatFloat m1, MatFloat m2, MatFloat m3)
{
  int i, j, k ;
  float d ;

#pragma omp parallel for private (d,i,j,k)   
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < M; ++j)
	{
	  d = 0.0 ;
	  for (k = 0; k < M; ++k)
	    d += m1 [i][k]* m2 [k][j] ;
	  m3 [i][j] = d ;
	}  
    }
}
 
