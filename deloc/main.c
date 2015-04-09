#include <time.h>
#include "vibronic.h"
#include <math.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(){
   
    // inter modes
    int nInter = 1;
    double wInter[] = {150, 250}; // mode 1 sym asym; mode 2 sym asym ...
    double bInter[] = {1.0, 1.0};

    double E[] = {0,0};
    double Vab = 300;
    int size_q = 2;
    int q[] = {2, 2};// quantum numbers

    // Memory for sparse Hamiltonian matrix 
    int *I, *J;
    double *VALUES;
    int numStates = get_prod(q, size_q); 
    printf("Number of states: %d\n", numStates);
    int numElems = 2*numStates + nInter*numStates*4;
    I = (int *)malloc(sizeof(int)*numElems);
    J = (int *)malloc(sizeof(int)*numElems);
    VALUES = (double *)malloc(sizeof(double)*numElems);

    // Compute Hamiltonian elements
    int elems=0;
    clock_t start = clock(), diff;
    int sec;
    SparseHamiltonian(nInter, wInter, bInter, q, size_q, E, Vab, I, J, VALUES, numStates, &elems);
    diff = clock() - start;
    sec = diff / CLOCKS_PER_SEC;
    printf("Sparse Hamiltonian computation time: %d seconds %d milliseconds\n", sec, sec/1000);
    
    printf("filling in zeros\n");
    // Expand to matrix
    double **M = (double **)malloc(sizeof(double *)*numStates);
    for (int i=0; i < numStates; i++) M[i] = (double *)malloc(sizeof(double)*numStates);
    for(int i=0; i < numStates; i++) {
          for(int j=0; j < numStates; j++) {
                   M[i][j] = 0;
                     }
    }
    
    for (int i=0; i < elems; i++) M[I[i]][J[i]] = VALUES[i];
    printf("printing\n");    
    // print matrix
    for(int i=0; i < numStates; i++) {
          for(int j=0; j < numStates; j++) printf("%3f ",M[i][j]);
          printf("\n");
    }

    free(I);
    free(J);
    free(VALUES);
    free(M);
    return 0; 
}
