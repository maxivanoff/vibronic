#include "multi.h"
#include <slepceps.h>
#include <time.h>
#include <unistd.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  SlepcInitialize(&argc,&argv,(char*)0,NULL);
    int nmodes[] = {3,3};
    int size_q = 7;
    double Vab = 20;
    int q[] = {2, 2,2,2, 2,2,2};// quantum numbers
    double w[] = {10, 20, 30, 10, 20, 30};//frequencies
    double b[] = {1, 0.5, 0.3, 1, 0.5, 0.3};//shifts
    double E[] = {0,0};
  double Msym = 1, Masym = 1;

  // Memory for sparse Hamiltonian matrix 
  int *I, *J;
  double *VALUES;
  int numStates = get_prod(q, size_q); 
  printf("Number of states: %d\n", numStates);
  int numElems = 2*numStates + nmodes[0]*2*numStates;
  I = (int *)malloc(sizeof(int)*numElems);
  J = (int *)malloc(sizeof(int)*numElems);
  VALUES = (double *)malloc(sizeof(double)*numElems);
  
  // Compute Hamiltonian elements
  int elems;
  clock_t start = clock(), diff;
  int sec;
  elems = SparseHamiltonian(nmodes, q, size_q, w, b, E, Vab, I, J, VALUES, numStates);
  diff = clock() - start;
  sec = diff / CLOCKS_PER_SEC;
  printf("Sparse Hamiltonian computation time: %d seconds %d milliseconds\n", sec, sec/1000);

  // Define matrix
  start = clock();
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  PetscInt       Istart,Iend,nev;
  PetscErrorCode ierr;

  PetscInt N = (PetscInt)numStates;
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  
  // Move data to the matrix
  for (int i=0; i<elems; i++){
  ierr = MatSetValue(A, I[i], J[i], VALUES[i], INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  diff = clock() - start;
  sec = diff / CLOCKS_PER_SEC;
  printf("Time to set up matrix: %d seconds %d milliseconds\n", sec, sec/1000);

  // Create eigensolver context
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

  // Set operators. In this case, it is a standard eigenvalue problem
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);

  // Set solver parameters at runtime
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
  // Smallest to largest eigenvalues 
  ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);CHKERRQ(ierr);

  // Solve the eigensystem
  start = clock();
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  diff = clock() - start;
  sec = diff / CLOCKS_PER_SEC;
  printf("Diagonalization time: %d seconds %d milliseconds\n", sec, sec/1000);
  
  // Print energy levels and intensities
  PetscReal levels;
  Vec Vr;
  PetscScalar *xx;
  double Isym, Iasym;

  VecCreateSeq(PETSC_COMM_SELF, numStates, &Vr);

  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  printf("E    Isym    Iasym\n");
  for (int i=0; i<nev; i++){
    ierr = EPSGetEigenvalue(eps, i, (PetscScalar*)&levels, PETSC_NULL);CHKERRQ(ierr); 
    ierr = EPSGetEigenvector(eps, i, Vr, PETSC_NULL);CHKERRQ(ierr);
    VecGetArray(Vr, &xx);
    Isym = Msym*pow(xx[0] + xx[numStates/2], 2);
    Iasym = Masym*pow(xx[0] - xx[numStates/2], 2);
    printf("%f %f %f\n", levels, Isym, Iasym);
  }

  // clean up
  free(I);
  free(J);
  free(VALUES);
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = VecDestroy(&Vr);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return 0;
}

