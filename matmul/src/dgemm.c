#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <mkl.h>

#define SEED 1234
#define min(x,y) (((x) < (y)) ? (x) : (y))

double random_number(); 

void  random_matrix(double*, int, int); 
void  zero_matrix(double*, int, int) ; 
void  print_matrix(double*, int , int, const char*); 

int main(int argc, char **argv) { 
    int    m, n, p; 
    struct timeval t1, t2; 
    double elapsed_time; 

    // matrix size 
    if (argc != 4) {
        m = 4; n = 4; p = 4; 
    } else {  
        m = atoi(argv[1]); n = atoi(argv[2]); p = atoi(argv[3]); 
    } 

    // scaling factors
    double alpha = 1.0, beta = 0.0; 
    
    // allocation 
    double* A = (double*) malloc(sizeof(double) * m * p);
    double* B = (double*) malloc(sizeof(double) * p * n);
    double* C = (double*) malloc(sizeof(double) * m * n);

    // initialize C 
    zero_matrix(C, m, n); 

    // initialize A, B 
    // incorrect result ! 
    srand(SEED); 
    random_matrix(A, m, p); 
    random_matrix(B, p, n); 
    
    // start timing 
    gettimeofday(&t1, NULL);

    // blas3 
    /* dgemm("N", "N", &m, &n, &p, &alpha, A, &p, B, &n, &beta, C, &n); */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, p, alpha, A, p, B, n, beta, C, n);

    // end timing 
    gettimeofday(&t2, NULL);

    // walltime 
    elapsed_time = (t2.tv_usec - t1.tv_usec)*1e-6 + (t2.tv_sec - t1.tv_sec);
    printf("Timing: %10.3f (s)\n", elapsed_time); 

    // gflops
    double gflops = (2.0*m*n*p - 1.0*m*p)*1E-9;
    printf("Performance: %10.3f (GFlops)\n", gflops/elapsed_time);
    
    // debug
    print_matrix(A, m, p, "A ="); 
    print_matrix(B, p, n, "B ="); 
    print_matrix(C, m, n, "C ="); 
    
    free(A); 
    free(B); 
    free(C); 

    return 0;
}

void random_matrix(double *matrix, int m, int n) { 
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            matrix[i*n + j] = random_number();
} 

void zero_matrix(double *matrix, int m, int n ) { 
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            matrix[i*n + j] = 0.0; 
} 

void print_matrix(double *matrix, int m , int n, const char *name ) { 
    printf("%s\n", name); 
    for (int i=0; i<min(m,4); i++) {
        for (int j=0; j<min(n,4); j++) {
            printf ("%12.5f", matrix[i*n+j]);
        }
        printf ("\n");
    }
} 

double random_number() { 
    return ((double)rand() / (double)RAND_MAX);   
} 
