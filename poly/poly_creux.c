#include <stdio.h>
#include <stdlib.h>

#include "polycreux.h"
#include <math.h>

#include <x86intrin.h>

p_polyf_creux_t creer_polyf_creux(int taille){
    p_polyf_creux_t p ; 
    p = (p_polyf_creux_t)malloc(sizeof(polyf_creux_t));
    p->coeffs = (float *)malloc(sizeof(float) * (taille));
    p->degres = (int *)malloc(sizeof(int) * taille);
    p->taille = taille; 

    return p; 
}

void init_polyf_creux(p_polyf_creux_t p, float x){
    register unsigned int i ;
    for (i = 0 ; i < p->taille; i++){
        p->coeffs [i] = x ;
        p->degres [i] = i ;
    }
}

void detruire_polyf_creux(p_polyf_creux_t p){
    free(p->coeffs);
    free(p->degres);
    free(p);
}

void ecrire_polynome_float_creux(p_polyf_creux_t p){
    if (p->degres[0]==0){
        printf("%f", p->coeffs[0]);
    }
    for (int i = 0; i< p->taille; i++){
         printf("+ %f X^%d", p->coeffs[i], i);
    }
    printf("\n"); 
}

p_polyf_creux_t lire_polyf_creux(char *name){
    FILE *f;
    p_polyf_creux_t p;
    int nb_degre;
    int i;
    int cr;
  
    f = fopen (name, "r") ;
    if (f == NULL){
        fprintf (stderr, "erreur ouverture %s \n", name);
        exit (-1);
        }
    
    cr = fscanf (f, "%d", &nb_degre);
    if (cr != 1){
        fprintf (stderr, "erreur lecture du nombre de degre\n");
        exit (-1);
        }

    p = creer_polynome_creux (nb_degre);
    
    for (i = 0 ; i < nb_degre; i++) {
            cr = fscanf (f, "%f", &p->degres[i]);
            if (cr != 1) {
                fprintf (stderr, "erreur lecture coefficient ou degree %d\n", i);
                exit (-1);
            }
            cr = fscanf (f, "%f", &p->coeffs[i]);
        }

    fclose (f);
    return p;
}

int egalite_polynome_float_creux(p_polyf_creux_t p1, p_polyf_creux_t p2){
    if (p1->taille!= p2->taille)return 0; 
    for(int i = 0 ; i<p1->taille;i++){
        if (p1->degres[i]!= p2->degres[i])return 0;
        if (p1->coeffs[i]!= p2->coeffs[i])return 0;
    }
    return  1;
}

p_polyf_creux_t addition_polyf_creux(p_polyf_creux_t p1, p_polyf_creux_t p2){
    //recherche du nombre de degres necessaires au resultat
	int taille_nouveau = p1->taille;
	int present_dans_deux = 0;
	for(int i=0; i<p2->taille; i++){
		present_dans_deux = 0;
		for(int y=0; y<p1->taille; y++){
			if(p2->degres[i] == p1->degres[y]) present_dans_deux=1;
		}
		if(!present_dans_deux) taille_nouveau++;
	}
	p_polyf_creux_t res = creer_polynome_creux(taille_nouveau);
	//rechreche du degre par lequel commencer
	int degre_min = p1->degres[0];
	if(p2->degres[0] < degre_min) degre_min=p2->degres[0];

	for(int i=0; i<res->taille; i++){
		res->degres[i]=degre_min;
		res->coeffs[i]=0;
		//ajout ce qu'il y a avec ce degre dans p1
		for(int x=0; x<p1->taille; x++){
			if(p1->degres[x]==degre_min) res->coeffs[i]+=p1->coeffs[x];
		}
		//ajout ce qu'il y a avec ce degre dans p2
		for(int x=0; x<p2->taille; x++){
			if(p2->degres[x]==degre_min) res->coeffs[i]+=p2->coeffs[x];
		}

		//recherche du prochain degre
		int pro_degre_p1 = degre_min;
		for(int x=0; x<p1->taille-1; x++){
			if(p1->degres[x]<=degre_min && p1->degres[x+1]>=degre_min) pro_degre_p1=p1->degres[x+1];
		}
		int pro_degre_p2 = degre_min;
		for(int x=0; x<p2->taille-1; x++){
			if(p2->degres[x]<=degre_min && p2->degres[x+1]>=degre_min) pro_degre_p2=p2->degres[x+1];
		}
		if(pro_degre_p2<pro_degre_p1)degre_min=pro_degre_p2; else degre_min=pro_degre_p1;
	}
}

p_polyf_creux_t mult_polyf_creux_scalaire(p_polyf_creux_t p, float n){
    p_polyf_creux_t res = creer_polyf_creux(p->taille); 
    res->degres = p->degres; 
    for (int i = 0; i < p->taille; i++){
        res->coeffs[i] = p->coeffs[i] * n;
    }
    return res;
} 


p_polyf_creux_t mult_polyf(p_polyf_creux_t p1, p_polyf_creux_t p2){
    int new_taille = p1->taille + p2->taille; 
    p_polyf_creux_t res = creer_polyf_creux(new_taille); 
    for (int i =0; i<p1->taille; i++){
        for (int j = 0; j<p2->taille; j++){
            res->degres[i+j] = p1->degres[i]+p2->coeffs[j];
            res->coeffs[i+j] = p1->coeffs[i]+p2->coeffs[j]; 
        }
    }
    return res;
}

p_polyf_creux_t puiss_poly_creux(p_polyf_creux_t p, int n){
    p_polyf_creux_t res = p; 
    for (int i =0; i<p->taille; i++){
        res = mult_polyf(res,p); 
    }
    return res; 
}