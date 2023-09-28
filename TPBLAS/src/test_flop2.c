#include <stdio.h>
#include "flop2.h"

#define N 8192

typedef float VecFloat [N] ;

VecFloat V1, V2;
float d ;


int main (int argc, char **argv)
{
  struct timespec start, end ;
  int i ;
  double duree ;
  
  init_flop_nano () ;

  TOP_NANO (start) ;

     for (i = 0 ; i < N ; i++)
         V2 [i] = V1 [i] ;
     
  TOP_NANO (end) ;

  duree = tdiff_nano (&start, &end) ;
  printf ("duree = %e copie de %ld octets\n", duree,  N*sizeof(float)) ;
  calcul_memop_nano ("Débit mémoire ", N * sizeof (float), duree) ;

  TOP_NANO (start) ;

     for (i = 0 ; i < N ; i++)
       d = (V2 [i] * V1 [i]) + 10.0 ;
     
  TOP_NANO (end) ;

  duree = tdiff_nano (&start, &end) ;
  printf ("duree = %e dot de %ld octets\n", duree,  (long int) 2 * N) ;
  calcul_flop_nano ("Produit scalaire ", 2 * N , duree) ;
}
