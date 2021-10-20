#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define SEED 1234

#define N 20000
#define G 6.674e-11
#define SOFT 1e-20
#define DISTANCE 1e9 

#define T 10
#define DT 1 

int main() { 
    // mass
    double* m  = (double*) malloc(N*sizeof(double)); 
    
    // position
    double* x  = (double*) malloc(N*sizeof(double)); 
    double* y  = (double*) malloc(N*sizeof(double)); 
    double* z  = (double*) malloc(N*sizeof(double)); 

    // velocity
    double* vx = (double*) malloc(N*sizeof(double)); 
    double* vy = (double*) malloc(N*sizeof(double)); 
    double* vz = (double*) malloc(N*sizeof(double)); 
    
    // force 
    double* fx = (double*) malloc(N*sizeof(double)); 
    double* fy = (double*) malloc(N*sizeof(double)); 
    double* fz = (double*) malloc(N*sizeof(double)); 

    // verlet update 
    double* dpx = (double*) malloc(N*sizeof(double)); 
    double* dpy = (double*) malloc(N*sizeof(double)); 
    double* dpz = (double*) malloc(N*sizeof(double)); 

    // initialization 
    srand(SEED); 

    for (int i = 0; i < N; i++) {
        // mass of earth
        m[i]  = 5.97e24; 
        
        x[i]  = DISTANCE * (double)rand() / (double)RAND_MAX; 
        y[i]  = DISTANCE * (double)rand() / (double)RAND_MAX; 
        z[i]  = DISTANCE * (double)rand() / (double)RAND_MAX; 

        vx[i] = 0.0; 
        vy[i] = 0.0; 
        vz[i] = 0.0;  
    } 
   
    int t, i, j; 
    for (t = 1; t <= T; t += DT) {
        double dx, dy, dz; 
        double dvx, dvy, dvz; 
        double ax, ay, az; 
        double d2, d32;  
        double F; 

        for (i = 0; i < N; i++) { 
            // force initialization 
            fx[i] = 0.0;  
            fy[i] = 0.0;  
            fz[i] = 0.0;  

            // force calculation 
            for (j = 0; j < N; j++) { 
                // (G * m[i] * m[j]/d2) * (d_j - d_j)/|d|
                // x component: (G *m[i] * m[j]) * dx * 1/d^3
                // y component: (G *m[i] * m[j]) * dy * 1/d^3
                // z component: (G *m[i] * m[j]) * dz * 1/d^3
                F  = G*m[i] * m[j]; 
                dx = x[j] - x[i]; 
                dy = y[j] - y[i]; 
                dz = z[j] - z[i]; 

                d2 = dx*dx + dy*dy + dz*dz + SOFT;
                d32 = 1/pow(d2, 3.0/2.0);

                fx[i] += F * dx * d32; 
                fy[i] += F * dy * d32; 
                fz[i] += F * dz * d32; 
            } 
        
            // acceleration of ith particle
            ax = fx[i] / m[i]; 
            ay = fy[i] / m[i]; 
            az = fz[i] / m[i]; 

            // Verlet
            dvx = vx[i] + ax*DT/2; 
            dvy = vy[i] + ay*DT/2; 
            dvz = vz[i] + az*DT/2; 

            dpx[i] = dvx * DT;  
            dpy[i] = dvy * DT;  
            dpz[i] = dvz * DT;  

            vx[i] = dvx; 
            vy[i] = dvy; 
            vz[i] = dvz; 
        }
    
        // update position 
        for (int i = 0; i < N; i++) { 
            x[i] += dpx[i]; 
            y[i] += dpy[i]; 
            z[i] += dpz[i]; 
        }

        printf("%-4d %12.5f %12.5f %12.5f\n", t, x[0], y[0], z[0]); 
    }

    free(m); 
    free(x); 
    free(y); 
    free(z); 
    free(vx); 
    free(vy); 
    free(vz); 
    free(fx); 
    free(fy); 
    free(fz); 
    free(dpx); 
    free(dpy); 
    free(dpz); 

    return 0; 
} 
