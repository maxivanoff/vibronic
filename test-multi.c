#include <stdlib.h> 
#include <stdio.h>
#include <string.h>

int print_1d(int v[], int size_v) {
    // prints a one-dimensional array v
    for (int i=0; i < size_v; i++) {
        printf("%i ",v[i]);
    }
    printf("\n");
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
int get_prod(int q[], int size_q) {
    // returns a product of all elements in q
    int qprod = 1;
    for (int i=0; i < size_q; i++) qprod *= q[i];
    return qprod;
}
int zeros(int v[], int s){
    for (int i=0; i<s; i++) v[i] = 0;
    return 0;
}

int set_qRem(int qRem[], int q[], int imode, int nmodes){
    int j=0;
    for (int i=1; i<nmodes*2+1; i++){
        if ((i != 1+imode*2) && (i != imode*2+2)) {
            qRem[j] = q[i];
            j++;
        }
    }
    return 0;
}
int pack_to_index(int v[], int q_base[], int size_q) {
    // gets an index for v
    // q_base: factors needed to calculate an index from array and reverse
    // v and q_base are assumed to be of the same length
    int result = 0;
    for (int i = size_q-1; i >= 0; i--) result += v[i]*q_base[i];
    return result;
}


double VDiag(int mode[], double w[], int imode){
    return w[2*imode]*(mode[0] + 0.5) + w[2*imode+1]*(mode[1] + 0.5);
}


int main(){
    
    int nmodes = 2;
    int size_q = 5;
    int q[5] = {1, 3,3, 2,2};// quantum numbers
    double w[4] = {10, 10, 20, 20};
    int nRem = (nmodes - 1)*2;
    int iv, imode;
    int qMode[2];
    
    int *RemProd, *MODE, *FullProd, *qRem;
    MODE = (int *)malloc(sizeof(int)*2);
    RemProd = (int *)malloc(sizeof(int)*nRem);
    FullProd = (int *)malloc(sizeof(int)*size_q);
    qRem = (int *)malloc(sizeof(int)*nRem);
    
    int numStates = get_prod(q, nmodes*2+1);
    int *I, *J;
    double *VALUES;

    I = (int *)malloc(sizeof(int)*numStates);
    J = (int *)malloc(sizeof(int)*numStates);
    VALUES = (double *)malloc(sizeof(double)*numStates);
    for (int i=0; i<numStates; i++) VALUES[i] = 0;

    int statesPerMode;
    int stateInd;
    FullProd[0] = 0; // First monomer
    for (imode=0; imode<nmodes; imode++){
        zeros(MODE, 2);
        qMode[0] = q[1+imode*2];
        qMode[1] = q[2+imode*2];
        for (iv=0; iv<q[1+imode*2]*q[1+imode*2]; iv++){
            printf("current mode:\n");
            print_1d(MODE, 2);
            
            zeros(RemProd, nRem);
            zeros(FullProd, nmodes*2);
            set_qRem(qRem, q, imode, nmodes);
            print_1d(qRem, nRem);
            statesPerMode = get_prod(qRem, nRem);
            printf("States per mode: %d\n", statesPerMode);
            for(int j=0; j<statesPerMode; j++){
                // Place current MODE at proper position                
                // Place remaining modes
                int irem=0;
                for (int ifull=0; ifull<nmodes; ifull++){
                    if (ifull == imode){
                        FullProd[1+2*imode] = MODE[0];
                        FullProd[2+2*imode] = MODE[1];}
                    else{
                        FullProd[1+2*ifull] = RemProd[2*irem];
                        FullProd[2+2*ifull] = RemProd[2*irem+1];
                        irem++;}
                }

                // diagonal
                stateInd = pack_to_index(FullProd, q, size_q);
                VALUES[stateInd] += VDiag(MODE, w, imode);      
                I[stateInd] = stateInd; 
                J[stateInd] = stateInd; 

                print_1d(FullProd, size_q);
                increase(RemProd, qRem, nRem);
            }

            increase(MODE, qMode, 2);

        }
    }

    // Expand to matrix
    double **H = (double **)malloc(sizeof(double *)*numStates);
    for (int i=0; i < numStates; i++) H[i] = (double *)malloc(sizeof(double)*numStates);
    for(int i=0; i < numStates; i++) {
          for(int j=0; j < numStates; j++) {
                   H[i][j] = 0;
                     }
    }
    for (int i=0; i < numStates; i++) {
        H[I[i]][J[i]] = VALUES[i];
        printf("%d %f\n",i, VALUES[i]); }

    // print matrix
    for(int i=0; i < numStates; i++) {
          for(int j=0; j < numStates; j++) printf("%10f ",H[i][j]);
          printf("\n");
    }


    free(RemProd);
    free(MODE);
    free(FullProd);
    free(I);
    free(J);
    free(VALUES);
    free(H);

    return 0;
}


