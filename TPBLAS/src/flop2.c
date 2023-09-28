#include <stdio.h>
#include "flop2.h"

// 2.6 Ghz sur la machine corte
// duree du cycle processeur en nano seconde (10e9)

static const double duree_cycle = (double) 1 / (double) 2.6 ;

static double residu_nano ;

double tdiff_nano (struct timespec *start,
	     struct timespec *end)
 {
   double t ;

   t = (end->tv_sec - start->tv_sec)
         +
        1e-9*(end->tv_nsec - start->tv_nsec) ;
   return (t - residu_nano) ;
 }

void init_flop_nano ()
{
  struct timespec start, end ;

  TOP_NANO(start) ;

  TOP_NANO(end) ;

  residu_nano = tdiff_nano (&start, &end) ; 
  printf ("residu_nano = %e \n", residu_nano) ;
}

void calcul_flop_nano (char *message, long int nb_operations_flottantes, double duree)
{
  printf ("%s %ld operations %.7f secondes Performance %5.3f GFLOP/s\n", message, nb_operations_flottantes, duree, ((((double) nb_operations_flottantes) / duree))/1e9) ;
  return ;
}

void calcul_memop_nano (char *message, long int nb_octets, double duree)
{
  printf ("%s %ld octets %.7f secondes Performance %5.3f GB/s\n", message, nb_octets, duree, ((((double) nb_octets) / duree))/1e9) ;
  return ;
}


