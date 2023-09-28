#include <time.h>
#include <sys/time.h>


#define TOP_NANO(x)  clock_gettime (CLOCK_MONOTONIC, &x) ;


void init_flop_nano () ;

double tdiff_nano (struct timespec *start, struct timespec *end) ;

void calcul_flop_nano (char *message, long int nb_operations_flottantes, double duree) ;

void calcul_memop_nano (char *message, long int nb_operations_flottantes, double duree) ;

