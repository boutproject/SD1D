/**************************************************************************
 * Serial solver used for calculating sheath boundary conditions in non-local electron closures
 * solves A.x=b for an n by n matrix A using the PETSc library
 *
 **************************************************************************
 * Copyright 2014 J.T.Omotani
 *
 * Contact: John Omotani, john.omotani@ccfe.ac.uk
 * 
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

class SerialMatrixSolve;

#ifndef __SERIAL_MATRIX_SOLVE_H__
#define __SERIAL_MATRIX_SOLVE_H__

#ifndef BOUT_HAS_PETSC

#include <boutexception.hxx>

class SerialMatrixSolve {
public:
  SerialMatrixSolve() { throw BoutException("PETSc not available"); }
};

#else

#include <globals.hxx>
#include <output.hxx>
#include <petscksp.h>
#include <options.hxx>
#include <bout/petsclib.hxx>
#include <boutexception.hxx>

class SerialMatrixSolve {
public:
  SerialMatrixSolve(const int &pass_n, void* ctx, PetscErrorCode mult(Mat,Vec,Vec), Options *opt = NULL);
  ~SerialMatrixSolve() {
    KSPDestroy(&ksp);
    MatDestroy(&A);
    VecDestroy(&b);
    VecDestroy(&x);
  }
  
  PetscErrorCode solve(BoutReal* bdata, BoutReal* xdata);
  
private:
  int n; // The size of the matrix system to be solved
  PetscLib lib;
  KSP ksp;
  Mat A; // The matrix
  Vec b; // The RHS vector
  Vec x; // The solution vector
  
  int* idx; // Array of indices. Since the matrix here is dense, this is just 0..(n-1)
  KSPConvergedReason ksp_converged_reason;
};

#endif //BOUT_HAS_PETSC

#endif //__SERIAL_MATRIX_SOLVE_H__
