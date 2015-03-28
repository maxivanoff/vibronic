#include "multi.h"
#include <math.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>

int main(){
   
    int nIntra[] = {0,0};
    // intra modes
    double wIntra[] = {0, 0};//frequencies
    double bIntra[] = {0, 0};//shifts
    // inter modes
    int nInter = 1;
    double wInter[] = {100,100}; // mode 1 sym asym; mode 2 sym asym ...
    double bInter[] = {0.8,0.8};

    double E[] = {0,0};
    double Vab = 300;
    int size_q = 2;
    int q[] = {2, 5};// quantum numbers

    // Memory for sparse Hamiltonian matrix 
    int *I, *J;
    double *VALUES;
    int numStates = get_prod(q, size_q); 
    printf("Number of states: %d\n", numStates);
    int numElems = 2*numStates + nIntra[0]*2*numStates + nInter*numStates*4;
    I = (int *)malloc(sizeof(int)*numElems);
    J = (int *)malloc(sizeof(int)*numElems);
    VALUES = (double *)malloc(sizeof(double)*numElems);
    
    // Compute Hamiltonian
    int elems;
    SparseHamiltonian(nIntra, nInter, wIntra, bIntra, wInter, bInter, q, size_q, E, Vab, I, J, VALUES, numStates, &elems);
    
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
