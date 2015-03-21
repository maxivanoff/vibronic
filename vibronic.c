#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <string.h>
#include "qbasic.h"

double VDiag(int v[], double w[], double E[], int nmodes){
    double value = E[v[0]];
    for (int i=0; i<nmodes;i++){
        value += w[2*i]*(v[1 + 2*i] + 0.5) + w[2*i+1]*(v[2*i + 2] + 0.5);
    }
    return value;
}


int Hamiltonian(int d, double M[][d]) {

    const size_t size_q = 3;
    int q[] = {2,2,2};
    int dim = get_prod(q, size_q);
    int q_base[size_q];
    get_q_base(q_base, q, size_q);
    
    print_1d(q, size_q);
    printf("%i\n", dim);
    print_1d(q_base, size_q);
    
    int nmodes = 1;
    double w[2] = {10, 10};
    double E[2] = {0,0};
    double Vab = 99;
    double b[2] = {1,1};
    
    int *I, *J;
    double *VALUES;

    I = (int *)malloc(sizeof(int)*dim*4);
    J = (int *)malloc(sizeof(int)*dim*4);
    VALUES = (double *)malloc(sizeof(double)*dim*4);
    
    for (int i=0;i<dim*2;i++) VALUES[i] = 1.1;



  

    int v[size_q] = {0};
    //print_1d(v, size_q);
    //printf("%i\n",pack_to_index(v, q_base, size_q));
    int Vsize = 0;
    int v2[size_q];
    int Ij;
    for (int i=0; i < dim; i++) {
        
        Ij = pack_to_index(v, q_base, size_q);
        // diagonal
        VALUES[Vsize] = VDiag(v, w, E, nmodes);
        I[Vsize] = Ij;
        J[Vsize] = I[Vsize];
        Vsize++;

        // coupling
        memcpy(v2, v, sizeof(v));
        v2[0] = 1 - v[0];
        VALUES[Vsize] = Vab;
        I[Vsize] = Ij;
        J[Vsize] = pack_to_index(v2, q_base, size_q);
        Vsize++;

        // n' = n - 1
        memcpy(v2, v, sizeof(v));
        int t = 1 + v2[0];
        if (v2[t]-- > 0) {

            VALUES[Vsize] = b[v2[0]] * w[v2[0]] * sqrt(v[t]);
            I[Vsize] = Ij;
            J[Vsize] = pack_to_index(v2, q_base, size_q);
            Vsize++;

            // n' = n + 1 is implented using symmetry
            VALUES[Vsize] = VALUES[Vsize-1];
            I[Vsize] = J[Vsize-1];
            J[Vsize] = I[Vsize-1];
            Vsize++;
        }



        increase(v, q, size_q);

    }

    // Expand to matrix
    //double **M = (double **)malloc(sizeof(double *)*dim);
    //for (int i=0; i < dim; i++) M[i] = (double *)malloc(sizeof(double)*dim);
    for(int i=0; i < dim; i++) {
          for(int j=0; j < dim; j++) {
                   M[i][j] = 0;
                     }
    }
    for (int i=0; i < Vsize; i++) M[I[i]][J[i]] = VALUES[i];

    // print matrix
    //for(int i=0; i < dim; i++) {
    //      for(int j=0; j < dim; j++) printf("%10f ",M[i][j]);
    //      printf("\n");
    //}

    // for (int i=0;i<Vsize;i++) printf("%d %d %f\n", I[i], J[i], VALUES[i]);

    free(VALUES);
    free(I);
    free(J);
    //free(M);
    return 0;
}
