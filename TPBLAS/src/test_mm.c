#include <stdio.h>
#include "flop3.h"
#include "bib3.h"


MatFloat M1, M2, M3 ;

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
      init_matrix (M2, 3.0) ;
      init_matrix (M3, 0.0) ;

      TOP_NANO (start) ;

          mult_mat_mat (M1, M2, M3) ;
    
      TOP_NANO (end) ;

      duree = tdiff_nano (&start, &end) ;
      printf ("duree = %e dot de %ld opérations flottantes \n", duree,  (long int) 2 * M * M) ;
      calcul_flop_nano ("Produit mat mat ",  2 * M * M * M , duree) ;  
    }

}
