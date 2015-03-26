#include <stdio.h>

struct data{
    int nmodes[2];
    int size_q;
    int q[5];// quantum numbers
    double w[4];//frequencies
    double b[4];//shifts
    double E[2];
    double Vab;
};


extern int get_prod(int q[], int size_q);

extern int SparseHamiltonian( int *nmodes, int *q, int size_q, double *w, double *b, double *E, double Vab, int *I, int *J, double *VALUES, int numStates);

