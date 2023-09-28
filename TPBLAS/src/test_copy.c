#include <stdio.h>
#include "flop3.h"
#include "bib3.h"

VecFloat V1, V2 ;

int main (int argc, char **argv)
{
  struct timespec start, end ;
  int i ;
  double duree ;
  float d ;
  init_flop_nano () ;

  for (i = 0; i < 10; i++)
    {
      init_vector (V1, 1.0) ;
      init_vector (V2, 2.0) ;

      TOP_NANO (start) ;

         copy (V1, V2) ;
  
	 TOP_NANO (end) ;

      duree = tdiff_nano (&start, &end) ;
      printf ("duree = %e copie de %ld octets\n", duree,  N*sizeof(float)) ;
      calcul_memop_nano ("Débit mémoire ", N * sizeof (float), duree) ;
    }

}
