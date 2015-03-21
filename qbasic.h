#include <stdio.h>
 

extern int get_prod(int q[], size_t size_q);

extern int increase(int v[], int q[], size_t size_q);

extern int get_q_base(int q_out[], int q[], size_t size_q);

extern int pack_to_index(int v[], int q_base[], size_t size_q);

extern int unpack_index(int k, int out_v[], int q_base[], size_t size_q);

extern int print_1d(int v[], size_t size_v);

