
#define N 8192

typedef float VecFloat [N] ;

void init_vector (VecFloat V, float valeur) ;

void copy (VecFloat V1, VecFloat V2) ;

float dot (VecFloat V1, VecFloat V2) ;

#define M 512

typedef float MatFloat [M][M] ;
typedef float VectMatFloat [M] ;

void init_matrix (MatFloat m, float v) ;

void mult_mat_vect (MatFloat M1, VectMatFloat V1, VectMatFloat VR) ; 

void mult_mat_mat (MatFloat M1, MatFloat M2, MatFloat M3) ;




