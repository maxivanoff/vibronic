#include <math.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>

int get_q_base(int q_out[], int q[], int size_q) {
    // calculates factors needed to calculate an index from array and reverse
    // q_out will contain the result
    q_out[size_q-1] = 1;
    for (int i = size_q-1; i > 0; i--) q_out[i-1] = q_out[i] * q[i];
    return 0;
}


int increase(int v[], int q[], int size_q) {
    // finds the next element of v
    // q: number of possible values per each position
    // v and q are assumed to be of the same length
    for (int i=size_q-1; i >= 0; i--) {
        if (q[i] == ++v[i]) v[i] = 0; else return 0;
    }
    return 1;
}

int print_1d(int v[], int size_v) {
    // prints a one-dimensional array v
    for (int i=0; i < size_v; i++) {
        printf("%i ",v[i]);
    }
    printf("\n");
    return 0;
}
int get_prod(int q[], int size_q) {
    // returns a product of all elements in q
    int qprod = 1;
    for (int i=0; i < size_q; i++) qprod *= q[i];
    return qprod;
}
int pack_to_index(int v[], int q_base[], int size_q) {
    // gets an index for v
    // q_base: factors needed to calculate an index from array and reverse
    // v and q_base are assumed to be of the same length
    int result = 0;
    for (int i = size_q-1; i >= 0; i--) result += v[i]*q_base[i];
    return result;
}

double EDiag(int state[], double w[], double E[], int nmodes){
    double En = E[state[0]];
    for (int i=0; i<nmodes;i++){
        En += w[2*i]*(state[1 + 2*i] + 0.5) + w[2*i+1]*(state[2*i + 2] + 0.5);
    }
    return En;
}

int main(){
    
    int nmodes = 2;
    int size_q = 5;
    int q[5] = {2, 2,2, 2,2};// quantum numbers
    double w[4] = {10, 10, 20, 20};//frequencies
    double b[4] = {1, 1, 0.5, 0.5};//shifts
    double E[2] = {0,0};
    double Vab = 99;
    
    int *q_base;
    q_base = (int *)malloc(sizeof(int)*size_q);
    get_q_base(q_base, q, size_q);
    
    // Memory for sparse Hamiltonian matrix 
    int *I, *J;
    double *VALUES;
    int numStates = get_prod(q, size_q); 
    printf("Number of states: %d\n", numStates);
    int numElems = numStates*4;
    I = (int *)malloc(sizeof(int)*numElems);
    J = (int *)malloc(sizeof(int)*numElems);
    VALUES = (double *)malloc(sizeof(double)*numElems);
    
    // Initialize first state
    int *state, *stateCopy;
    state = (int *)malloc(sizeof(int)*size_q);
    stateCopy = (int *)malloc(sizeof(int)*size_q);
    for(int i=0; i<size_q; i++) state[i]=0;

    int Ij;
    int pos;
    int elems=0; // at the loop end contains elemsber of nonzero elements
    for (int istate=0; istate<numStates; istate++){
        Ij = pack_to_index(state, q_base, size_q);
        
        printf("i=%d ", Ij);
        print_1d(state, size_q);
        
        // diagonal
        VALUES[elems] = EDiag(state, w, E, nmodes);
        I[elems] = Ij;
        J[elems] = Ij;
        elems++;
        
        // coupling
        for (int j=0; j<size_q; j++) stateCopy[j] = state[j];
        stateCopy[0] = 1 - state[0];
        VALUES[elems] = Vab;
        I[elems] = Ij;
        J[elems] = pack_to_index(stateCopy, q_base, size_q);
        elems++;
        
        // n' = n - 1
        for (int imode=0; imode<nmodes; imode++){
            pos = 1 + state[0] + 2*imode;
            for (int j=0; j<size_q; j++) stateCopy[j] = state[j];
            if (stateCopy[pos]-- > 0) {

                VALUES[elems] = b[pos-1] * w[pos-1] * sqrt(state[pos]);
                I[elems] = Ij;
                J[elems] = pack_to_index(stateCopy, q_base, size_q);
                elems++;

                // n' = n + 1 is implemented using symmetry
                VALUES[elems] = VALUES[elems-1];
                I[elems] = J[elems-1];
                J[elems] = I[elems-1];
                elems++;
            }
        }

        increase(state, q, size_q);

}

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
          for(int j=0; j < numStates; j++) printf("%10f ",M[i][j]);
          printf("\n");
    }

   
    free(q_base);
    free(stateCopy);
    free(state);
    free(I);
    free(J);
    free(VALUES);
    free(M);
} 
