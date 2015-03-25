#include "multi.h"
#include <math.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>


int main(){
    
    int nmodes = 2;
    int size_q = 5;
    int q[] = {2, 2,2, 2,2};// quantum numbers
    double w[] = {10, 10, 20,20};//frequencies
    double b[] = {1, 1, 0.5, 0.5};//shifts
    double E[] = {0,0};
    double Vab = 99;
    
   
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
    elems = SparseHamiltonian(nmodes, q, size_q, w, b, E, Vab, I, J, VALUES, numStates);
    
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
