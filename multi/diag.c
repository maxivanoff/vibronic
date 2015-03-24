/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2014, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.

   SLEPc is free software: you can redistribute it and/or modify it under  the
   terms of version 3 of the GNU Lesser General Public License as published by
   the Free Software Foundation.

   SLEPc  is  distributed in the hope that it will be useful, but WITHOUT  ANY
   WARRANTY;  without even the implied warranty of MERCHANTABILITY or  FITNESS
   FOR  A  PARTICULAR PURPOSE. See the GNU Lesser General Public  License  for
   more details.

   You  should have received a copy of the GNU Lesser General  Public  License
   along with SLEPc. If not, see <http://www.gnu.org/licenses/>.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Standard symmetric eigenproblem corresponding to the Laplacian operator in 2 dimensions.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x dimension.\n"
  "  -m <m>, where <m> = number of grid subdivisions in y dimension.\n\n";

#include "multi.h"
#include <slepceps.h>
#include <time.h>
#include <unistd.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  SlepcInitialize(&argc,&argv,(char*)0,help);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute Hamiltonian matrix 
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   
  int nmodes = 3;
  int size_q = 7;
    
  int q[] = {2, 3,3, 3,3, 2,2};// quantum numbers
  double w[] = {10, 10, 20, 20, 30, 30};//frequencies
  double b[] = {1, 1, 0.5, 0.5, 0.3, 0.3};//shifts
  double E[] = {0,0};
  double Vab = 20;

  // Memory for sparse Hamiltonian matrix 
  int *I, *J;
  double *VALUES;
  int numStates = get_prod(q, size_q); 
  printf("Number of states: %d\n", numStates);
  int numElems = numStates*50;
  I = (int *)malloc(sizeof(int)*numElems);
  J = (int *)malloc(sizeof(int)*numElems);
  VALUES = (double *)malloc(sizeof(double)*numElems);
  
  // Compute Hamiltonian elemts
  int elems;
  clock_t start = clock(), diff;
  double sec;
  elems = SparseHamiltonian(nmodes, q, size_q, w, b, E, Vab, I, J, VALUES, numStates);
  diff = clock() - start;
  sec = diff / CLOCKS_PER_SEC;
  printf("Sparse Hamiltonian computation time: %f seconds %f milliseconds\n", sec, sec/1000);


  // Expand to matrix
  double **M = (double **)malloc(sizeof(double *)*numStates);
  for (int i=0; i < numStates; i++) M[i] = (double *)malloc(sizeof(double)*numStates);
  for(int i=0; i < numStates; i++) {
      for(int j=0; j < numStates; j++) {
          M[i][j] = 0;
      }
  }
  for (int i=0; i < elems; i++) M[I[i]][J[i]] = VALUES[i];
  
  //  Define EPS matrix
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  EPSType        type;
  PetscInt       Istart,Iend,nev,ii,jj;
  PetscErrorCode ierr;

  PetscInt N = (PetscInt)numStates;
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  
  for (ii=0; ii<numStates; ii++){
    for (jj=0; jj<numStates; jj++){
        ierr = MatSetValue(A, ii, jj, M[ii][jj], INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);

  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);CHKERRQ(ierr);
  start = clock();
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  diff = clock() - start;
  sec = diff / CLOCKS_PER_SEC;
  printf("Diagonalization time: %f seconds %f milliseconds\n", sec, sec/1000);

  /*
     Optional: Get some information from the solver and display it
  */
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);
  PetscReal levels;
  Vec Vr;
  PetscScalar *xx;
  double Isym, Iasym;
  int iB;

  VecCreateSeq(PETSC_COMM_SELF, numStates, &Vr);

  printf("E    Isym    Iasym\n");
  for (ii=0; ii<nev; ii++){
    ierr = EPSGetEigenvalue(eps, ii, (PetscScalar*)&levels, PETSC_NULL);CHKERRQ(ierr); 
    ierr = EPSGetEigenvector(eps, ii, Vr, PETSC_NULL);CHKERRQ(ierr);
    VecGetArray(Vr, &xx);
    iB = numStates/2;
    Isym = pow(xx[0] + xx[iB], 2);
    Iasym = pow(xx[0] - xx[iB], 2);
    printf("%f %f %f\n", levels, Isym, Iasym);
  }


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  //ierr = EPSPrintSolution(eps,NULL);CHKERRQ(ierr);
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = VecDestroy(&Vr);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return 0;
}

