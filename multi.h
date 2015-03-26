#include <stdio.h>

struct data{
    int nmodesA;
    int nmodesB;
    int size_q;
    int q[5];// quantum numbers
    double wA[2];//frequencies
    double wB[2];//frequencies
    double bA[2];//shifts
    double bB[2];//shifts
    double E[2];
    double Vab;
};


extern int get_prod(int q[], int size_q);

extern int SparseHamiltonian( struct data *params, int I[], int J[], double VALUES[], int numStates);

