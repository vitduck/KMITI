#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

#define SEED 1234
#define min(x,y) (((x) < (y)) ? (x) : (y))

#if defined DOUBLE
    typedef double PREC;
#else
    typedef float PREC;
#endif

PREC random_number(); 

void  random_matrix(PREC*, int, int); 
void  zero_matrix(PREC*, int, int) ; 
void  print_matrix(PREC*, int , int, const char*); 
void  mat_mul(PREC*, PREC*, PREC*, int, int, int); 

int main(int argc, char **argv) { 
    int  m, n, p; 
    double t1, t2, elapsed_time; 
    
    // matrix size 
    if (argc != 4) {
        m = 4; n = 4; p = 4; 
    } else {  
        m = atoi(argv[1]); n = atoi(argv[2]); p = atoi(argv[3]); 
    } 

    // allocation
    PREC* A = (PREC*) malloc(sizeof(PREC) * m * p); 
    PREC* B = (PREC*) malloc(sizeof(PREC) * p * n); 
    PREC* C = (PREC*) malloc(sizeof(PREC) * m * n); 

    //initialize A, B 
    srand(SEED); 
    random_matrix(A, m, p); 
    random_matrix(B, p, n); 

    //initialize C
    zero_matrix(C, m, n); 

    // start timing 
    
    // C = A * B 
    t1 = omp_get_wtime(); 
    mat_mul(A, B, C, m, n, p); 
    t2 = omp_get_wtime(); 

    // walltime 
    elapsed_time = t2 - t1; 
    printf("Timing: %10.3f (s)\n", elapsed_time); 

    // gflops 
    PREC gflops = (2.0*m*n*p - 1.0*m*p)*1E-9; 
    printf("Performance: %10.3f (GFlops)\n", gflops/elapsed_time);
   
    // debug
    print_matrix(A, m, p, "A ="); 
    print_matrix(B, p, n, "B ="); 
    print_matrix(C, m, n, "C ="); 
    
    //deallocate 
    free(A); 
    free(B); 
    free(C); 

    return 0;
}

void random_matrix(PREC *matrix, int m, int n) { 
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            matrix[i*n + j] = random_number();
} 

void zero_matrix(PREC *matrix, int m, int n ) { 
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            matrix[i*n + j] = 0.0; 
} 

void mat_mul(PREC *A, PREC *B, PREC *C, int m, int n, int p) { 
    int i, j, k; 
    #pragma omp parallel for shared(A, B, C) private(i, j, k)
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)  
            for (k = 0; k < p; k++) 
                C[i*n+j] += A[i*p+k] * B[k*n+j]; 
} 

void print_matrix(PREC *matrix, int m , int n, const char *name ) { 
    printf("%s\n", name); 
    
    for (int i=0; i<min(m,4); i++) {
        for (int j=0; j<min(n,4); j++) {
            printf ("%12.5f", matrix[i*n+j]);
        }
        printf ("\n");
    }
} 

PREC random_number() { 
    return ((PREC)rand() / (PREC)RAND_MAX);   
} 
