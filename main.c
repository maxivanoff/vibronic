#include "vibronic.h"


int main(){
    
    const int dim=8;
    double M[dim][dim];
    Hamiltonian(dim, M);
    
    // print matrix
    for(int i=0; i < dim; i++) {
          for(int j=0; j < dim; j++) printf("%10f ",M[i][j]);
          printf("\n");
    }

  
}
