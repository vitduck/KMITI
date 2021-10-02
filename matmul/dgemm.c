#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>

#include <mkl.h>

#define SEED 1337
#define SIZE 4000

double random_number(); 
void   random_matrix(double *matrix, int row, int col); 
void   zero_matrix(double *matrix, int row, int col); 
void   print_matrix(double *matrix, int row , int col, const char *name); 
void   mat_mul(double *A, double *B, double *C, int row_A, int col_A, int row_B, int col_B); 

int main() { 
    struct timeval t1, t2; 
    double elapsed_time; 

    double alpha = 1.0, beta = 1.0; 
    char   transA = 'N', transB = 'N'; 

    int row_A = SIZE, col_A = SIZE; 
    int row_B = SIZE, col_B = SIZE; 

    assert(col_A == row_B); 

    double* A = (double*) malloc(sizeof(double) * row_A * col_A); 
    double* B = (double*) malloc(sizeof(double) * row_B * col_B); 
    double* C = (double*) malloc(sizeof(double) * row_A * col_B); 

    // initialize A, B 
    srand(SEED); 
    random_matrix(A, row_A, col_A); 
    random_matrix(B, row_B, col_B); 

    // initialize C 
    zero_matrix(C, row_A, col_B); 

    gettimeofday(&t1, NULL);

    dgemm(
        &transA, &transB,
        &row_A, &col_B, &col_A, 
        &alpha, A, &row_A, B, &row_B, &beta, C, &row_A);

    gettimeofday(&t2, NULL);

    // walltime 
    elapsed_time = (t2.tv_usec - t1.tv_usec)*1e-6 + (t2.tv_sec - t1.tv_sec);
    printf("Timing: %10.3f (s)\n", elapsed_time); 

    // gflops 
    double gflops = (2.0*row_A*col_B*col_A + 1.0*row_A*col_B)*1E-9; 
    printf("Performance: %10.3f (GFlops)\n", gflops/elapsed_time);

    // debug
    /* print_matrix(A, row_A, col_A, "A =");  */
    /* print_matrix(B, row_B, col_B, "B =");  */
    /* print_matrix(C, row_A, col_B, "C =");  */

    free(A); 
    free(B); 
    free(C); 

    return 0;
}

void random_matrix(double *matrix, int row, int col) { 
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            matrix[i * col + j] = random_number();
} 

void zero_matrix(double *matrix, int row, int col ) { 
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            matrix[i * col + j] = 0.0; 
} 

void print_matrix(double *matrix, int row , int col, const char *name ) { 
    printf("%s\n", name); 
    
    for (int i = 0; i < row*col; i++) { 
        if (i > 0 && i % col == 0) { 
            printf("\n"); 
        }
        printf("%12.5f", matrix[i]); 
    }

    printf("\n"); 
} 

void mat_mul(double *A, double *B, double *C, int row_A, int col_A, int row_B, int col_B) { 
    for (int i = 0; i < row_A; i++)
        for (int j = 0; j < col_B; j++)  
            for (int k = 0; k < col_A; k++) 
                C[i * col_B +j] += A[i * col_A + k] * B[k * col_B + j]; 
} 

double random_number() { 
    return ((double)rand() / (double)RAND_MAX);   
} 
