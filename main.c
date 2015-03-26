#include "multi.h"
#include <math.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>

int main(){
   
    struct data params;
    params.nmodesA = 2;
    params.nmodesB = 2;
    params.size_q = 5;
    params.Vab = 99;
    int q[5] = {2, 2,1, 2,1};// quantum numbers
    double wA[2] = {10, 20};//frequencies
    double wB[2] = {10, 20};//frequencies
    double bA[2] = {1, 0.5};//shifts
    double bB[2] = {1, 0.5};//shifts
    double E[2] = {0,0};
    memcpy(params.q, q, sizeof params.q);
    memcpy(params.wA, wA, sizeof params.wA);
    memcpy(params.wB, wB, sizeof params.wB);
    memcpy(params.bA, bA, sizeof params.bA);
    memcpy(params.bB, bB, sizeof params.bB);
    memcpy(params.E, E, sizeof params.E);
    
   
    // Memory for sparse Hamiltonian matrix 
    int *I, *J;
    double *VALUES;
    int numStates = get_prod(params.q, params.size_q); 
    printf("Number of states: %d\n", numStates);
    int numElems = numStates*4;
    I = (int *)malloc(sizeof(int)*numElems);
    J = (int *)malloc(sizeof(int)*numElems);
    VALUES = (double *)malloc(sizeof(double)*numElems);
    
    // Compute Hamiltonian
    int elems;
    elems = SparseHamiltonian(&params, I, J, VALUES, numStates);
    //elems = SparseHamiltonian(nmodes, q, size_q, w, b, E, Vab, I, J, VALUES, numStates);
    
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
