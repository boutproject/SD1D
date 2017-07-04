/*!
 * \file non-local_parallel_serial_matrixsolve.cxx
 *
 * \brief Serial solver used for calculating sheath boundary conditions in non-local electron closures
 * solves A.x=b for an n by n matrix A using the PETSc library
 *
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
 */

#include "non-local_parallel_serial_matrixsolve.hxx"

SerialMatrixSolve::SerialMatrixSolve(const int &pass_n, void* ctx, PetscErrorCode mult(Mat,Vec,Vec), Options *opt) {
  
  PetscErrorCode ierr;
  
  n=pass_n;
  BoutReal rtol, atol, dtol;
  int maxits;
  BoutReal richardson_damping_factor, chebyshev_max,chebyshev_min;
  int gmres_max_steps;
  string ksptype;
  string pctype;

  // Get Options in Laplace Section
  Options *opts;
  if (!opt) opts = Options::getRoot()->getSection("serialmatrixsolve");
  else opts=opt->getSection("serialmatrixsolve");
  // Get Tolerances for KSP solver
  OPTION(opts,rtol,1.e-10);
  OPTION(opts,atol,1.e-50);
  OPTION(opts,dtol,1.e50);
  OPTION(opts,maxits,10000);
  // Get KSP Solver Type
  opts->get("ksptype", ksptype, "gmres");
  
  // Get preconditioner type
  // WARNING: only a few of these options actually make sense: see the PETSc documentation to work out which they are (possibly pbjacobi, sor might be useful choices?)
  opts->get("pctype", pctype, "none");
  
  // Get Options specific to particular solver types
  opts->get("richardson_damping_factor",richardson_damping_factor,1.0,true);
  opts->get("chebyshev_max",chebyshev_max,100,true);
  opts->get("chebyshev_min",chebyshev_min,0.01,true);
  opts->get("gmres_max_steps",gmres_max_steps,30,true);
 
  ierr = KSPCreate(PETSC_COMM_SELF, &ksp);CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = KSPSetInitialGuessNonzero(ksp, PetscBool(true));CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = KSPSetTolerances(ksp, rtol, atol, dtol, maxits);CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = KSPSetType(ksp, ksptype.c_str());CHKERRABORT(MPI_COMM_WORLD,ierr);
  if( ksptype == KSPRICHARDSON ) KSPRichardsonSetScale( ksp, richardson_damping_factor );
  #ifdef KSPCHEBYSHEV
    else if( ksptype == KSPCHEBYSHEV ) KSPChebyshevSetEigenvalues( ksp, chebyshev_max, chebyshev_min );
  #endif
  else if( ksptype == KSPGMRES ) KSPGMRESSetRestart( ksp, gmres_max_steps );
  
  PC pc;
  ierr = KSPGetPC(ksp, &pc);CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = PCSetType(pc, pctype.c_str());CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRABORT(MPI_COMM_WORLD,ierr);
  
  ierr = MatCreateShell(PETSC_COMM_SELF, n, n, n, n, ctx, &A);CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = MatShellSetOperation(A, MATOP_MULT, (void(*)(void))mult );CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, n, &b);CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecDuplicate(b, &x);CHKERRABORT(MPI_COMM_WORLD,ierr);
  
  idx = new int[n];
  for (int i=0; i<n; i++)
    idx[i]=i;
  
}

PetscErrorCode SerialMatrixSolve::solve(BoutReal* bdata, BoutReal* xdata) {
  //NB 'matrix' information is stored in the memory pointed to by ctx given to MatCreateShell
  
  PetscErrorCode ierr;
  
  ierr = VecSetValues(b, n, idx, bdata, INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecSetValues(x, n, idx, xdata, INSERT_VALUES);CHKERRQ(ierr);
  
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x);CHKERRQ(ierr);
  
#if PETSC_VERSION_GE(3,5,0)
  KSPSetOperators( ksp, A, A ); CHKERRQ(ierr);
#else
  KSPSetOperators( ksp, A, A,SAME_NONZERO_PATTERN ); CHKERRQ(ierr);
#endif
  
  ierr = KSPSolve(ksp, b, x);CHKERRQ(ierr);
  
  ierr = KSPGetConvergedReason(ksp, &ksp_converged_reason);CHKERRQ(ierr);
  if (ksp_converged_reason<0) {
    output<<"Failed to converge with ksp_converged_reason="<<ksp_converged_reason<<endl;
    throw BoutException("SerialMatrixSolve failed to converge");
  }
  
  ierr = KSPGetSolution(ksp, &x);CHKERRQ(ierr);
  
  ierr = VecGetValues(x, n, idx, xdata);CHKERRQ(ierr);
  
}
