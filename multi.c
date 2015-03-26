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

double Energy(int state[], double w[], int start, int nmodes){
    double E=0;
    for (int i=0; i<nmodes; i++){
        E += w[i]*(state[start + i] + 0.5);
    }
    return E;
}
int ladder(int state[], double w[], double b[], int start, int nmodes, int I[], int J[], double VALUES[], int size_q, int Ij, int q_base[], int *elems){
    int pos;
    int *stateCopy;
    stateCopy = (int *)malloc(sizeof(int)*size_q);
    for (int imode=0; imode<nmodes; imode++){
        pos = start + imode;
        for (int j=0; j<size_q; j++) stateCopy[j] = state[j];
        if (stateCopy[pos]-- > 0) {

            // n' = n - 1
            VALUES[*elems] = b[imode] * w[imode] * sqrt(state[pos]);
            I[*elems] = Ij;
            J[*elems] = pack_to_index(stateCopy, q_base, size_q);
            *elems++;

            // n' = n + 1 is implemented using symmetry
            VALUES[*elems] = VALUES[*elems-1];
            I[*elems] = J[*elems-1];
            J[*elems] = I[*elems-1];
            *elems++;
        }
    }
    return 0;
}
//int SparseHamiltonian(int nmodes, int q[], int size_q, double w[], double b[], double E[], double Vab, int I[], int J[], double VALUES[], int numStates){
int SparseHamiltonian( int *nmodes, int *q, int size_q, double *w, double *b, double *E, double Vab, int *I, int *J, double *VALUES, int numStates, int *elems){
    
    int *q_base;
    q_base = (int *)malloc(sizeof(int)*size_q);
    get_q_base(q_base, q, size_q);
    
    // Initialize first state
    int *state;
    state = (int *)malloc(sizeof(int)*size_q);
    for(int i=0; i<size_q; i++) state[i]=0;
    int *stateCopy;
    stateCopy = (int *)malloc(sizeof(int)*size_q);

    int Ij;
    int start;
    for (int istate=0; istate<numStates; istate++){
        Ij = pack_to_index(state, q_base, size_q);
        // diagonal
        VALUES[*elems] = E[*state] + Energy(state, w, 1, nmodes[0]+nmodes[1]);
        I[*elems] = Ij;
        J[*elems] = Ij;
        *elems++;
        
        // coupling
        for (int j=0; j<size_q; j++) stateCopy[j] = state[j];
        stateCopy[0] = 1 - *state;
        VALUES[*elems] = Vab;
        I[*elems] = Ij;
        J[*elems] = pack_to_index(stateCopy, q_base, size_q);
        *elems++;
        
        // subdiagonal
        start = 1;
        for (int i=0; i < *state; i++) start += nmodes[i];
        ladder(state, w, b, start, nmodes[*state], I, J, VALUES, size_q, Ij, q_base, elems);

        increase(state, q, size_q);

    }
    free(q_base);
    free(stateCopy);
    free(state);
    return 0;
} 
