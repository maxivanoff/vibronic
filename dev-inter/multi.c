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

double EIntra(int state[], double w[], double E[], int nIntra){
    double En = E[state[0]];
    int i;
    for (i=0; i<nIntra;i++){
        En += w[2*i]*(state[1 + 2*i] + 0.5) + w[2*i+1]*(state[2*i + 2] + 0.5);
    }
    return En;
}

double EInter(int state[], double w1[], double E[],  int nInter){
    double En = E[state[0]];
    int i;
    for (i=0; i<nInter; i++){
        En += w[i]*(state[1 + 2*nIntra + i] + 0.5); 
    }
    return En;
}


int SparseHamiltonian(int nIntra, int q[], int size_q, double w[], double b[], double E[], double Vab, int I[], int J[], double VALUES[], int numStates){
    
    int *q_base;
    q_base = (int *)malloc(sizeof(int)*size_q);
    get_q_base(q_base, q, size_q);
    
    // Initialize first state
    int *state, *stateCopy;
    state = (int *)malloc(sizeof(int)*size_q);
    stateCopy = (int *)malloc(sizeof(int)*size_q);
    for(int i=0; i<size_q; i++) state[i]=0;

    int Ij;
    int pos;
    int elems=0; // at the loop end contains number of nonzero elements

    //inter modes variables
    int nInter = 2;
    double l1[2] = [1, 1];
    double l2[2] = [1, 1];
    double w1[2] = [1, 1];
    double w2[2] = [1, 1];

    for (int istate=0; istate<numStates; istate++){
        Ij = pack_to_index(state, q_base, size_q);
        
        // diagonal block - diagonal
        VALUES[elems] = EIntra(state, w, E, nIntra) + EInter(state, w1, E, nInter);
        I[elems] = Ij;
        J[elems] = Ij;
        elems++;
        
        // diagonal block - subdiagonal INTRA
        for (int imode=0; imode<nIntra; imode++){
            pos = 1 + state[0] + 2*imode;
            for (int j=0; j<size_q; j++) stateCopy[j] = state[j];
            if (stateCopy[pos]-- > 0) { // n' = n - 1
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
        
        // diagonal block - subdiagonal INTER
        for (int imode=0; imode<nInter; imode++){
            pos = 1 + 2*nIntra + imode;
            for (int j=0; j<size_q; j++) stateCopy[j] = state[j];
            if (stateCopy[pos]-- > 0) { // n' = n - 1
                VALUES[elems] = l1[imode] * w1[imode] * sqrt(state[pos]);
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
         
        // coupling block - diagonal
        VALUES[elems] = Vab;
        for (int i=0; i<nInter; i++){
            VALUES[elems] += w2[i]*(state[1 + 2*nIntra + i] + 0.5)
        }
        for (int j=0; j<size_q; j++) stateCopy[j] = state[j];
        stateCopy[0] = 1 - state[0];
        I[elems] = Ij;
        J[elems] = pack_to_index(stateCopy, q_base, size_q);
        elems++;
        // coupling block - subdiagonal INTER
        for (int imode=0; imode<nInter; imode++){
            for (int j=0; j<size_q; j++) stateCopy[j] = state[j];
            pos = 1 + 2*nIntra + imode;
            if (stateCopy[pos]-- > 0){ // n' = n - 1
                VALUES[elems] = l2[imode] * w2[imode] * sqrt(state[pos]);
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

   
    free(q_base);
    free(stateCopy);
    free(state);
    return elems;
} 
