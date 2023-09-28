#include <stdio.h>
#include <stdlib.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define    NB_FOIS        512

int main (int argc, char **argv)
{
 complexe_float_t c1 = {1.0, 2.0} ;
 complexe_float_t c2 = {3.0, 6.0} ;

 complexe_double_t cd1 ;
 complexe_double_t cd2 ;

 complexe_double_t cm1 = (complexe_double_t) {1.0, 17.0} ;
 complexe_double_t cm2 = (complexe_double_t) {10.0, 4.0} ;

 complexe_float_t cmf1 = {2.0, 4.0}; 
 complexe_float_t cmf2 = {6.0, 6.0}; 

 complexe_double_t cdiv1 = (complexe_double_t) {5.0, 3.0} ;
 complexe_double_t cdiv2 = (complexe_double_t) {12.0, 9.0} ;

 complexe_float_t cdivf1 = {9.0, 8.0};
 complexe_float_t cdivf2 = {7.0, 6.0};

 struct timeval start, end ;
 
 int i ;

 init_flop_micro () ;
 
 c1 = add_complexe_float (c1, c2) ;

 printf ("c1.r %f c&.i %f\n", c1.real, c1.imaginary) ;



 cd1 = (complexe_double_t) {10.0, 7.0} ;
 cd2 = (complexe_double_t) {25.0, 32.0} ;


 TOP_MICRO(start) ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     cd1 = add_complexe_double (cd1, cd2) ;
   }

 TOP_MICRO(end) ;

 printf ("apres boucle cd1.real %f cd1.imaginary %f duree %f \n", cd1.real, cd1.imaginary, tdiff_micro (&start, &end)) ;

 calcul_flop_micro ("calcul complexe ", NB_FOIS*2, tdiff_micro(&start, &end)) ;
 exit (0) ;



//mult_complex_float
TOP_MICRO(start) ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     cmf1 = mult_complexe_float (cmf1, cmf2) ;
   }

 TOP_MICRO(end) ;

 printf ("apres boucle cd1.real %f cd1.imaginary %f duree %f \n", cmf1.real, cmf1.imaginary, tdiff_micro (&start, &end)) ;

 calcul_flop_micro ("calcul complexe ", NB_FOIS*2, tdiff_micro(&start, &end)) ;
 exit (0) ;


 //mult_complex_double
 TOP_MICRO(start) ;

  for (i = 0 ; i < NB_FOIS; i++)
   {
     cm1 = mult_complexe_double (cm1, cm2) ;
   }

 TOP_MICRO(end) ;

 printf ("apres boucle cd1.real %f cd1.imaginary %f duree %f \n", cm1.real, cm1.imaginary, tdiff_micro (&start, &end)) ;

 calcul_flop_micro ("calcul complexe ", NB_FOIS*2, tdiff_micro(&start, &end)) ;
 exit (0) ;


 // div complex float 
TOP_MICRO(start) ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     cdivf1 = mult_complexe_float (cdivf1, cdivf2) ;
   }

 TOP_MICRO(end) ;

 printf ("apres boucle cd1.real %f cd1.imaginary %f duree %f \n", cdivf1.real, cdivf1.imaginary, tdiff_micro (&start, &end)) ;

 calcul_flop_micro ("calcul complexe ", NB_FOIS*2, tdiff_micro(&start, &end)) ;
 exit (0) ;



 //div complex double 
TOP_MICRO(start) ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     cdiv1 = mult_complexe_double (cdiv1, cdiv2) ;
   }

 TOP_MICRO(end) ;

 printf ("apres boucle cd1.real %f cd1.imaginary %f duree %f \n", cdiv1.real, cdiv1.imaginary, tdiff_micro (&start, &end)) ;

 calcul_flop_micro ("calcul complexe ", NB_FOIS*2, tdiff_micro(&start, &end)) ;
 exit (0) ;

}
