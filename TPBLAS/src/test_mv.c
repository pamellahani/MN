#include <stdio.h>
#include "flop3.h"
#include "bib3.h"


MatFloat M1;
VectMatFloat V3, V4 ;

int main (int argc, char **argv)
{
  struct timespec start, end ;
  int i ;
  double duree ;
  float d ;
  init_flop_nano () ;
  
  for (i = 0; i < 10; i++)
    {
      init_matrix (M1, 2.0) ;

      TOP_NANO (start) ;

          mult_mat_vect (M1, V3, V4) ;
    
      TOP_NANO (end) ;

      duree = tdiff_nano (&start, &end) ;
      printf ("duree = %e dot de %ld opérations flottantes \n", duree,  (long int) 2 * M * M) ;
      calcul_flop_nano ("Produit mat vect ", 2 * M * M , duree) ;
    }
}
