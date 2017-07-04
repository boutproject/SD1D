/**************************************************************************
 * Solver to impose sheath boundary conditions in non-local electron closures
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

class SheathSolver;

#ifndef __SHEATH_SOLVER_H__
#define __SHEATH_SOLVER_H__

#ifndef BOUT_HAS_PETSC

#include <boutexception.hxx>

class SheathSolver {
public:
  SheathSolver() { throw BoutException("PETSc not available"); }
};

#else

#include <globals.hxx>
#include <output.hxx>
#include <petscksp.h>
#include <options.hxx>
#include <bout/petsclib.hxx>
#include <boutexception.hxx>
#include "non-local_parallel_serial_matrixsolve.hxx"

class Ematrix;

class SheathSolver {
public:
  SheathSolver(const int pass_n, const int pass_n_perp, BoutReal* pass_eigenvalues, const int moments_number_Hermite, const int moments_number_Laguerre,
	       BoutReal* const pass_deltaz, BoutReal* const pass_n_lower, BoutReal* const pass_n_upper, BoutReal* const pass_n_plus, BoutReal* const pass_n_minus,
	       BoutReal* const pass_n_plus_integral, BoutReal* const pass_n_minus_integral, BoutReal* const pass_n_zero_lower, BoutReal* const pass_n_zero_upper,
	       BoutReal* const pass_nS_lower, BoutReal* const pass_nS_upper, BoutReal* const pass_ddsslower_nS_lower, BoutReal* const pass_ddssupper_nS_lower,
	       BoutReal* const pass_ddsslower_nS_upper, BoutReal* const pass_ddssupper_nS_upper, BoutReal* const pass_ss_lower, BoutReal* const pass_ss_upper, Options* opt = NULL);
  ~SheathSolver();
  PetscErrorCode mult(Vec vec_in, Vec vec_out);
  PetscErrorCode solve();
  
private:
  SerialMatrixSolve* matrixsolver;
  bool approximate_ss_derivs;
  int n; // The size of the matrix system to be solved
  int n_perp; // The number of perpendicular (x-z) points to solve for on this processor
  BoutReal * eigenvalues;
  BoutReal * exp_deltaz_over_eigenvalues;
  
  BoutReal * deltaz;
  BoutReal * n_lower;
  BoutReal * n_upper;
  BoutReal * n_plus;
  BoutReal * n_minus;
  BoutReal * n_plus_integral;
  BoutReal * n_minus_integral;
  BoutReal * n_zero_lower;
  BoutReal * n_zero_upper;
  BoutReal * nS_lower;
  BoutReal * nS_upper;
  BoutReal * ddsslower_nS_lower;
  BoutReal * ddssupper_nS_lower;
  BoutReal * ddssupper_nS_upper;
  BoutReal * ddsslower_nS_upper;
  BoutReal * ss_lower;
  BoutReal * ss_upper;
  
  Ematrix * ematrix;
  
  BoutReal * nbuffer; // length-n buffer used for temporary storage during calculation of rhs and mult
  BoutReal * ddsslower_nminus; // arrays used for temporary storage in computation of ss_lower and ss_upper derivatives of nS_lower and nS_upper
  BoutReal * ddssupper_nminus;
  BoutReal * ddssupper_nplus;
  BoutReal * ddsslower_nplus;
  BoutReal * n_minus_previous;
  BoutReal * rhs;
  PetscErrorCode solve_one(const BoutReal &single_deltaz, const BoutReal single_n_lower, const BoutReal single_n_upper, BoutReal* single_n_plus, BoutReal* single_n_minus,
			   BoutReal* single_n_plus_integral, BoutReal* single_n_minus_integral, BoutReal* single_n_minus_previous, BoutReal &single_n_zero_lower, BoutReal &single_n_zero_upper,
			   BoutReal &single_nS_lower, BoutReal &single_nS_upper, BoutReal &single_ddsslower_nS_lower, BoutReal &single_ddssupper_nS_lower,
			   BoutReal &single_ddsslower_nS_upper, BoutReal &single_ddssupper_nS_upper);
};

class Ematrix {
public:
  Ematrix(const int n, const int moments_number_Hermite, const int moments_number_Laguerre);
  ~Ematrix();
  BoutReal * E_lower;
  BoutReal * E_upper;
  BoutReal * E0_lower;
  BoutReal * E0_upper;
  BoutReal * zero_eigenvalue_E_lower;
  BoutReal * zero_eigenvalue_E_upper;
  BoutReal zero_eigenvalue_E0_lower;
  BoutReal zero_eigenvalue_E0_upper;
  BoutReal * nS_E_lower;
  BoutReal * nS_E_upper;
  BoutReal nS_E0_lower;
  BoutReal nS_E0_upper;
  BoutReal * E_lower_prime;
  BoutReal * E_upper_prime;
  BoutReal * E0_lower_prime;
  BoutReal * E0_upper_prime;
  BoutReal * nS_E_lower_prime;
  BoutReal * nS_E_upper_prime;
  BoutReal nS_E0_lower_prime;
  BoutReal nS_E0_upper_prime;
  void calculate_E(const BoutReal &ss_lower, const BoutReal& ss_upper);
private:
  int n;
  int nfiles;
  int interpolation_index_lower;
  int interpolation_index_upper;
  BoutReal * ss_min_array;
  BoutReal * ss_max_array;
  BoutReal * E_array;
  BoutReal * Eprime_array;
  BoutReal * E0_array;
  BoutReal * E0prime_array;
  BoutReal * zero_eigenvalue_E_array;
  BoutReal * zero_eigenvalue_Eprime_array;
  BoutReal * zero_eigenvalue_E0_array;
  BoutReal * zero_eigenvalue_E0prime_array;
  BoutReal * nS_E_array;
  BoutReal * nS_Eprime_array;
  BoutReal * nS_E0_array;
  BoutReal * nS_E0prime_array;
};

#endif //BOUT_HAS_PETSC

#endif //__SHEATH_SOLVER_H__
