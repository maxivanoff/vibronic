#include "multi.h"
#include <math.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>

int main(){
   
    int nmodes[] = {2,2};
    int size_q = 5;
    double Vab = 99;
    int q[5] = {2, 2,1, 2,1};// quantum numbers
    double E[2] = {0,0};
    // intra modes
    double w[] = {10, 20, 10, 20};//frequencies
    double b[] = {1, 0.5, 1, 0.5};//shifts
    // inter modes
    int lmodes = 1;
    double w_sym[] = {2};
    double w_asym[] = {2};
    double b_sym[] = {1};
    double b_sym[] = {1};

    double *w1, *w2, *b1, *b2;
    w1 = (double *)malloc(sizeof(int)*lmodes);
    w2 = (double *)malloc(sizeof(int)*lmodes);
    b1 = (double *)malloc(sizeof(int)*lmodes);
    b2 = (double *)malloc(sizeof(int)*lmodes);

    for (int i=0; i<lmodes; i++){
        w1[i] = sqrt((pow(w_sym, 2) + pow(w_asym, 2))/2);
        w2[i] = sqrt((pow(w_sym, 2) - pow(w_asym, 2))/2);
        b1[i] = (b_sym[i] + b_asym[i])/2;
        b2[i] = (b_sym[i] - b_asym[i])/2;
    }
   
    // Memory for sparse Hamiltonian matrix 
    int *I, *J;
    double *VALUES;
    int numStates = get_prod(q, size_q); 
    printf("Number of states: %d\n", numStates);
    int numElems = numStates*4;
    I = (int *)malloc(sizeof(int)*numElems);
    J = (int *)malloc(sizeof(int)*numElems);
    VALUES = (double *)malloc(sizeof(double)*numElems);
    
    // Compute Hamiltonian
    int elems;
    SparseHamiltonian(nmodes, q, size_q, w, b, E, Vab, I, J, VALUES, numStates, &elems, w1, w2, b1, b2, lmodes);
    
    // Expand to matrix
    double **M = (double **)malloc(sizeof(double *)*numStates);
    for (int i=0; i < numStates; i++) M[i] = (double *)malloc(sizeof(double)*numStates);
    for(int i=0; i < numStates; i++) {
          for(int j=0; j < numStates; j++) {
                   M[i][j] = 0;
                     }
    }
    
    for (int i=0; i < elems; i++) M[I[i]][J[i]] = VALUES[i];
    
    // print matrix
    for(int i=0; i < numStates; i++) {
          for(int j=0; j < numStates; j++) printf("%3f ",M[i][j]);
          printf("\n");
    }

    free(I);
    free(J);
    free(VALUES);
    free(M);
 
}
