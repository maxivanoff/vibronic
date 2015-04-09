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

double Energy(int state[], double w[], int start, int nmodes){
    double E=0;
    for (int i=0; i<nmodes; i++){
        E += w[i]*(state[start + i] + 0.5);
    }
    return E;
}
int ladder(int state[], double bw[], int start, int nmodes, int I[], int J[], double VALUES[], int size_q, int Ij, int q_base[], int *elems, int q[]){
    int pos;
    int *stateCopy;
    stateCopy = (int *)malloc(sizeof(int)*size_q);
    for (int imode=0; imode<nmodes; imode++){
        pos = start + imode;
        // n' = n - 1
        for (int j=0; j<size_q; j++) stateCopy[j] = state[j];
        if (stateCopy[pos]-- > 0) {
            VALUES[*elems] = bw[imode] * sqrt(state[pos]);
            I[*elems] = Ij;
            J[*elems] = pack_to_index(stateCopy, q_base, size_q);
            ++*elems;
        }
        // n' = n + 1
        for (int j=0; j<size_q; j++) stateCopy[j] = state[j];
        if (stateCopy[pos]++ < q[pos]-1) {
            VALUES[*elems] = bw[imode] * sqrt(state[pos]+1);
            I[*elems] = Ij;
            J[*elems] = pack_to_index(stateCopy, q_base, size_q);
            ++*elems;
        }
    }
    return 0;
}
int SparseHamiltonian( int nIntra, int nInter, double wIntra[], double bIntra[], double wInter[], double bInter[], int q[], int size_q, double E[], double Vab, int I[], int J[], double VALUES[], int numStates, int *elems){
    int nUnits=2;
    
    int *q_base;
    q_base = (int *)malloc(sizeof(int)*size_q);
    get_q_base(q_base, q, size_q);
    
    int *Istate, *Jstate;
    Istate = (int *)malloc(sizeof(int)*size_q);
    Jstate = (int *)malloc(sizeof(int)*size_q);
    // Initialize first state
    for(int i=0; i<size_q; i++) Istate[i]=0;

    double *bwIntra, *w;
    bwIntra = (double *)malloc(sizeof(double) * nIntra);
    w = (double *)malloc(sizeof(double) * nIntra);
    
    // Convert inter
    double *wPlus, *wMinus, *bwPlus, *bwMinus;
    wPlus = (double *)malloc(sizeof(double) * nInter);
    wMinus = (double *)malloc(sizeof(double) * nInter);
    bwPlus = (double *)malloc(sizeof(double) * nInter);
    bwMinus = (double *)malloc(sizeof(double) * nInter);
    for (int i=0; i<nInter; i++){
        wPlus[i] = (wInter[2 * i] + wInter[2 * i + 1]) * 0.5;
        wMinus[i] = (wInter[2 * i] - wInter[2 * i + 1]) * 0.5;
        bwPlus[i] = (bInter[2 * i] * wInter[2 * i] + bInter[2 * i + 1] * wInter[2 * i + 1]) * 0.5;
        bwMinus[i] = (bInter[2 * i] * wInter[2 * i] - bInter[2 * i + 1] * wInter[2 * i + 1]) * 0.5;
    }

    int Iindex;
    int start;
    for (int i=0; i<numStates; i++){
        Iindex = pack_to_index(Istate, q_base, size_q);
        for (int u=0; u<nUnits; u++){
            if (u == *Istate){ // diagonal block
                // diagonal
                // E0 + inter vibrations
                VALUES[*elems] = E[u] + Energy(Istate, wPlus, 1 + nUnits * nIntra, nInter);
                // intra from all units
                for (int k=0; k<nUnits; k++){
                    for(int kk=0; kk<nIntra; kk++) w[kk] = wIntra[kk * nUnits + k];
                    VALUES[*elems] += Energy(Istate, w, 1 + k * nIntra, nIntra);
                }
                I[*elems] = Iindex;
                J[*elems] = Iindex;
                ++*elems;
                
                // subdiagonal
                // INTRA
                start = 1 + u * nIntra;
                for(int k=0; k<nIntra; k++) bwIntra[k] = wIntra[k * nUnits + u] * bIntra[k * nUnits + u];
                ladder(Istate, bwIntra, start, nIntra, I, J, VALUES, size_q, Iindex, q_base, elems, q);
                // INTER
                ladder(Istate, bwPlus, 1 + nUnits * nIntra, nInter, I, J, VALUES, size_q, Iindex, q_base, elems, q);
            }
            else{ // coupling block
                for (int j=0; j<size_q; j++) Jstate[j] = Istate[j];
                *Jstate = u;

                // Vab constant + INTER 
                VALUES[*elems] = Vab + Energy(Istate, wMinus, 1 + nUnits * nIntra, nInter);
                I[*elems] = Iindex;
                J[*elems] = pack_to_index(Jstate, q_base, size_q);
                ++*elems;
                
                // INTER only
                ladder(Jstate, bwMinus, 1 + nUnits * nIntra, nInter, I, J, VALUES, size_q, Iindex, q_base, elems, q);
            }
        }

        increase(Istate, q, size_q);

    }
    free(q_base);
    free(Istate);
    free(Jstate);
    free(bwIntra);
    free(w);
    free(wPlus);
    free(wMinus);
    free(bwPlus);
    free(bwMinus);
    return 0;
} 
