#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define N 10000
#define Max 500001
#define Max_Loop 1

double Uniform(void) {
    return ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
}

void model(double **f2, double **f1, double **x, double V, double dt, double a) {
    int i, j, m;

    for (i = 0; i < N; i++) {
        f2[i][0] = (x[i][1] + a * f1[i][1]) * dt;
        f2[i][1] = (-V - 3 * pow(x[i][0] + a * f1[i][0] - 1, 2) * (x[i][0] + a * f1[i][0] - 2) +
                    (-1 + (2 / (1 + 20 * pow(x[i][0] + a * f1[i][0] - 1, 2)))) * (x[i][1] + a * f1[i][1])) * dt - 
                    2 * 0.1 * (x[i][1] + a * f1[i][1]) * dt - 2 * (x[i][0] + a * f1[i][0]) * dt;

        for (j = -1; j <= 1; j++) {
            if (j != 0) {
                m = i + j;
                if (m < 0) {
                    m = N + m;
                } else if (m >= N) {
                    m = m - N;
                }
                f2[i][1] += 0.1 * (x[m][1] + a * f1[m][1]) * dt + (x[m][0] + a * f1[m][0]) * dt;
            }
        }
    }
}

int main(void) {
    FILE *data;
    char data_file[256];
    double V, h = 0.01;
    int i, j, n0, n;
    double **x, **I, **f1, **f2, **f3, **f4;
    int loop;

    x = (double **)malloc(sizeof(double *) * N);
    I = (double **)malloc(sizeof(double *) * N);
    f1 = (double **)malloc(sizeof(double *) * N);
    f2 = (double **)malloc(sizeof(double *) * N);
    f3 = (double **)malloc(sizeof(double *) * N);
    f4 = (double **)malloc(sizeof(double *) * N);

    for (i = 0; i < N; i++) {
        x[i] = (double *)malloc(sizeof(double) * 10);
        I[i] = (double *)malloc(sizeof(double) * 10);
        f1[i] = (double *)malloc(sizeof(double) * 10);
        f2[i] = (double *)malloc(sizeof(double) * 10);
        f3[i] = (double *)malloc(sizeof(double) * 10);
        f4[i] = (double *)malloc(sizeof(double) * 10);
    }

    sprintf(data_file, "output.dat");
    printf("%s", data_file);
    data = fopen(data_file, "a");

    for (loop = 0; loop < Max_Loop; loop++) {
        srand(loop);

        // Initial condition
        for (i = 0; i < N; i++) {
            I[i][0] = 0;
            I[i][1] = 0;
        }

        V = 0.309;

        // Initialize
        for (i = 0; i < N; i++) {
            for (j = 0; j < 2; j++) {
                x[i][j] = I[i][j];
                f1[i][j] = I[i][j];
                f2[i][j] = I[i][j];
                f3[i][j] = I[i][j];
                f4[i][j] = I[i][j];
            }
        }

        // Initial dynamics
        for (n0 = 0; n0 < Max; n0++) {
            // 4th Runge-Kutta method
            model(f1, x, x, V, h, 0);
            model(f2, f1, x, V, h, 0.5);
            model(f3, f2, x, V, h, 0.5);
            model(f4, f3, x, V, h, 1);
            for (i = 0; i < N; i++) {
                x[i][0] = x[i][0] + (f1[i][0] + 2 * f2[i][0] + 2 * f3[i][0] + f4[i][0]) / 6;
                x[i][1] = x[i][1] + (f1[i][1] + 2 * f2[i][1] + 2 * f3[i][1] + f4[i][1]) / 6;
                // Noise
                if (0.0001 > Uniform()) {
                    x[i][0] = 0;
                }
            }
        }

        // Dynamics for output
        for (n = 0; n < Max; n++) {
            // 4th Runge-Kutta method
            model(f1, x, x, V, h, 0);
            model(f2, f1, x, V, h, 0.5);
            model(f3, f2, x, V, h, 0.5);
            model(f4, f3, x, V, h, 1);

            for (i = 0; i < N; i++) {
                x[i][0] = x[i][0] + (f1[i][0] + 2 * f2[i][0] + 2 * f3[i][0] + f4[i][0]) / 6;
                x[i][1] = x[i][1] + (f1[i][1] + 2 * f2[i][1] + 2 * f3[i][1] + f4[i][1]) / 6;
                // Noise
                if (0.0001 > Uniform()) {
                    x[i][0] = 0;
                }
            }

            if(n%10==0){
                for(i=0;i<N;i++){
                    if (x[i][0]<1){
                        fprintf(data,"0\t");
                    }else{
                        fprintf(data,"255\t");
                    }
                }
                fprintf(data,"\n");
            }
        }
    }
    fclose(data);

    // Free allocated memory
    for (i = 0; i < N; i++) {
        free(x[i]);
        free(I[i]);
        free(f1[i]);
        free(f2[i]);
        free(f3[i]);
        free(f4[i]);
    }
    free(x);
    free(I);
    free(f1);
    free(f2);
    free(f3);
    free(f4);

    return 0;
}

