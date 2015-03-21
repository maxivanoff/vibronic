#include <stdio.h>
 

int get_prod(int q[], size_t size_q) {
    // returns a product of all elements in q
    int qprod = 1;
    for (int i=0; i < size_q; i++) qprod *= q[i];
    return qprod;
}


int get_q_base(int q_out[], int q[], size_t size_q) {
    // calculates factors needed to calculate an index from array and reverse
    // q_out will contain the result
    q_out[size_q-1] = 1;
    for (int i = size_q-1; i > 0; i--) q_out[i-1] = q_out[i] * q[i];
    return 0;
}


int increase(int v[], int q[], size_t size_q) {
    // finds the next element of v
    // q: number of possible values per each position
    // v and q are assumed to be of the same length
    int x;
    for (int i=size_q-1; i >= 0; i--) {
        // printf("%i\n",i);
        if (q[i] == ++v[i]) v[i] = 0; else return 0;
    }
    return 1;
}



int print_1d(int v[], size_t size_v) {
    // prints a one-dimensional array v
    for (int i=0; i < size_v; i++) {
        printf("%i ",v[i]);
    }
    printf("\n");
    return 0;
}



int pack_to_index(int v[], int q_base[], size_t size_q) {
    // gets an index for v
    // q_base: factors needed to calculate an index from array and reverse
    // v and q_base are assumed to be of the same length
    int result = 0;
    for (int i = size_q-1; i >= 0; i--) result += v[i]*q_base[i];
    return result;
}



int unpack_index(int k, int out_v[], int q_base[], size_t size_q) {
    // reverses convert_to_index
    // k: number to parse
    // out_v: output array
    // q: number of possible values per each position
    // out_v and q are assumed to be of the same length

    for (int i = 0; i < size_q; i++) {
        out_v[i] = k / q_base[i];
        k %= q_base[i];
    }

    return 0;
}
