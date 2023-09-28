#ifndef _POLY_CREUX_H_
#define _POLY_CREUX_H_

typedef struct {
  int taille; 
  int *degres;          
  float *coeffs;  
} polyf_creux_t ,*p_polyf_creux_t ;


p_polyf_creux_t creer_polyf_creux(int degre); 
void init_polyf_creux(p_polyf_creux_t p, float x); 
void detruire_polyf_creux(p_polyf_creux_t p); 

void ecrire_polynome_float_creux(p_polyf_creux_t p); 
p_polyf_creux_t lire_polynome_float_creux(char *name); 

int egalite_polynome_float_creux(p_polyf_creux_t p1, p_polyf_creux_t p2);   
p_polyf_creux_t addition_polyf_creux(p_polyf_creux_t p1, p_polyf_creux_t p2); 
p_polyf_creux_t mult_polyf_creux_scalaire(p_polyf_creux_t p, float n); 
p_polyf_creux_t mult_polyf(p_polyf_creux_t p1, p_polyf_creux_t p2); 
p_polyf_creux_t puiss_poly_creux(p_polyf_creux_t p, int n); 


#endif