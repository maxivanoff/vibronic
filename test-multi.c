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
int modeUp(int v[], int q) {
    // finds the next element of v
    // q: number of excitations
    // single mode is defined by two quantum numbers, thus
    // v is assumed to be of length = 2
    for (int i=1; i >= 0; i--) {
        if (q == ++v[i]) v[i] = 0; else return 0;
    }
    return 1;
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
}

int main(){
    
    int nmodes = 3;
    int q[7] = {2, 3,3, 2,2, 1,1};// quantum numbers
    int nRem = (nmodes - 1)*2;
    int iv, imode;
    
    int *RemProd, *MODE, *FullProd, *qRem;
    MODE = (int *)malloc(sizeof(int)*2);
    RemProd = (int *)malloc(sizeof(int)*nRem);
    FullProd = (int *)malloc(sizeof(int)*nmodes*2);
    qRem = (int *)malloc(sizeof(int)*nRem);

    for (imode=0; imode<nmodes; imode++){
        zeros(MODE, 2);
        for (iv=0; iv<q[1+imode*2]*q[1+imode*2]; iv++){
            printf("current mode:\n");
            print_1d(MODE, 2);
            
            zeros(RemProd, nRem);
            zeros(FullProd, nmodes*2);
            set_qRem(qRem, q, imode, nmodes);
            print_1d(qRem, nRem);
            for(int j=0; j<get_prod(qRem, nRem); j++){
                // Place current MODE at proper position                
                FullProd[2*imode] = MODE[0];
                FullProd[2*imode+1] = MODE[1];

                // Place remaining modes
                
                //print_1d(FullProd, nmodes*2);
                increase(RemProd, qRem, nRem);
            }

            modeUp(MODE, q[1+imode]);

        }
    }

    free(RemProd);
    free(MODE);
    free(FullProd);

    return 0;
}


