#include <stdio.h>
#include <stdlib.h>

#include "poly.h"


int main (int argc, char **argv)
{
  p_polyf_t p1, p2; 
  p_polyf_t p3,p4, p5,p6;  
  
  if (argc != 3)
    {
      fprintf (stderr, "deux paramètres (polynomes,fichiers) sont à passer \n") ;
      exit (-1) ;
    }
      
  p1 = lire_polynome_float (argv [1]) ;
  p2 = lire_polynome_float (argv [2]) ;

  printf ("p1 ="); 
  ecrire_polynome_float (p1) ;
  printf ("p2 ="); 
  ecrire_polynome_float (p2) ;

  printf("multiplication entre 2 polynomes p1 et p1\n") ;
  p3 = multiplication_polynomes (p1, p1); 
  ecrire_polynome_float (p3) ;

  printf("multiplication polynome p1 et scalaire=(2)\n") ;
  p4 =  multiplication_polynome_scalaire (p1, 2.0); 
  ecrire_polynome_float (p4) ;

  printf("puissance du polynome p1 a la puissance 2:\n"); 
  p5 = puissance_polynome(p1,2); 
  ecrire_polynome_float (p5) ;

  printf ("composition du polynome p1 par p2\n") ;
  p6 = composition_polynome (p1, p2);
  ecrire_polynome_float (p6) ;
  /*
    ajouter du code pour tester les fonctions
    sur les polynomes
  */
}
