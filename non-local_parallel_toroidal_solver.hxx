/**************************************************************************
 * Inversion for non-local parallel closures on toroidal flux surfaces. 
 *                           Using MUMPS Solver
 *
 **************************************************************************
 * Copyright 2013 Copyright 2013 J. Omotani
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
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
class NonLocalParallelToroidalSolver;

#ifndef __MUMPS_LAPLACE_H__
#define __MUMPS_LAPLACE_H__

// #ifndef BOUT_HAS_MUMPS
// 
// #include <boutexception.hxx>
// 
// class NonLocalParallelToroidalSolver {
// public:
//   NonLocalParallelToroidalSolver() { throw BoutException("Mumps library not available"); }
// };
// 
// #else

#include <globals.hxx>
#include <output.hxx>
// #include <bout/constants.hxx>
// #include <options.hxx>
#include <boutexception.hxx>
#include "dmumps_c.h"

#define MUMPS_JOB_INIT -1
#define MUMPS_JOB_END -2
#define MUMPS_JOB_ANALYSIS 1
#define MUMPS_JOB_FACTORIZATION 2
#define MUMPS_JOB_SOLUTION 3
#define MUMPS_JOB_ANALYSIS_AND_FACTORIZATION 4 // combines job 1 and job 2
#define MUMPS_JOB_BOTH 5 // combines job 2 and job 3
#define MUMPS_JOB_ALL 6 // combines job 1, job 2 and job 3

class NonLocalParallelToroidalSolver {
public:
  NonLocalParallelToroidalSolver(bool pass_reverse);
  ~NonLocalParallelToroidalSolver() {
    if (is_initialized) {
      delete [] is_root;
      for (int i=0; i<nx; i++) {
	mumps_struc[i].job = -2;
	dmumps_c(&mumps_struc[i]);
	delete [] mumps_struc[i].irn_loc;
	delete [] mumps_struc[i].jcn_loc;
	delete [] mumps_struc[i].a_loc;
	delete [] mumps_struc[i].isol_loc;
	delete [] mumps_struc[i].rhs;
	if (reverse) MPI_Comm_free(&comms[i]);
      }
    }
  }
  
//   const FieldPerp solve(const FieldPerp &deltaI, const FieldPerp &decayfactor);
  void solve(BoutReal* deltaI, BoutReal* decayfactor, int jx, BoutReal shiftangle);


private:
  
//   FieldPerp sol;              // solution Field
  
  int n; // Order of the matrix
  int nx; // Number of points in x = number of solvers (mumps_struc's)
  
  // Istart is the first row of MatA owned by the process, Iend is 1 greater than the last row.
  int Istart, Iend; 

  bool reverse; // True if data starts on last y-processor in the mesh rather than first
  MPI_Comm* comms;

//   Options *opts;              // Options Object
  
  bool* is_root; // True on the root processor for this solver
  bool is_initialized;
  DMUMPS_STRUC_C* mumps_struc;
};

// #endif //BOUT_HAS_MUMPS

#endif //MUMPS_LAPLACE_H
