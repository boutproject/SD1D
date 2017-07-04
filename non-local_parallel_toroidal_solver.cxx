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
#include "non-local_parallel_toroidal_solver.hxx"

// #include "mpi.h"
#include <bout/sys/timer.hxx>
#include <boutcomm.hxx>
#include <msg_stack.hxx>
#include <cmath>

NonLocalParallelToroidalSolver::NonLocalParallelToroidalSolver(bool pass_reverse) {
  
  reverse = pass_reverse;
  
//   // Get Options in Laplace Section
//   if (!opt) opts = Options::getRoot()->getSection("laplace");
//   else opts=opt;
  
  // Get communicator for group of processors in Y.
  comms = new MPI_Comm[mesh->xend-mesh->xstart+1];
  is_root = new bool[mesh->xend-mesh->xstart+1];
  for (int i=0; i<(mesh->xend-mesh->xstart+1); i++) {
    int rank;
    MPI_Comm_rank(mesh->getYcomm(i+mesh->xstart),&rank);
    int commsize;
    MPI_Comm_size(mesh->getYcomm(i+mesh->xstart),&commsize);
    if (reverse) {
      MPI_Group group;
      MPI_Comm_group(mesh->getYcomm(i+mesh->xstart),&group);
      MPI_Group newgroup;
      int* ranks = new int[commsize];
      for (int j=0; j<commsize; j++)
	ranks[j] = commsize - j - 1;
      MPI_Group_incl(group,commsize,ranks,&newgroup);
      MPI_Comm_create(mesh->getYcomm(i+mesh->xstart),newgroup,&comms[i]);
      MPI_Group_free(&newgroup);
      MPI_Group_free(&group);
      if (rank == commsize-1) is_root[i] = true;
      else is_root[i] = false;
    }
    else {
      comms[i] = mesh->getYcomm(i+mesh->xstart);
      if (rank == 0) is_root[i] = true;
      else is_root[i] = false;
    }
  }

#ifdef CHECK
  if (mesh->LocalNz != mesh->GlobalNz) throw BoutException("NonLocalParallelToroidalSolve cannot handle grid split in z-direction at present");
#endif
  
  // Get implementation specific options
//   opts->get("fourth_order", fourth_order, false);

//   sol.allocate();
//   // The value of sol's y-index should be irrelevant, but needs to be set to something so checkData() does not complain
//   if(reverse) sol.setIndex(mesh->ystart);
//   else sol.setIndex(mesh->yend);
  
  // n is the order of the matrix to be inverted
  n = mesh->LocalNz;
  nx = mesh->xend-mesh->xstart+1;
  mumps_struc = new DMUMPS_STRUC_C[nx];
  for (int i=0; i<nx; i++) {
    mumps_struc[i].comm_fortran = (MUMPS_INT) MPI_Comm_c2f(comms[i]); // MPI communicator for MUMPS, in fortran format
    mumps_struc[i].sym = 0; // Solve using unsymmetric matrix
    mumps_struc[i].par = 1; // Use the host processor (rank 0) to do work for the solution
    
    // Initialize the mumps struc
    mumps_struc[i].job = MUMPS_JOB_INIT;
    dmumps_c(&mumps_struc[i]);
    
    mumps_struc[i].icntl[2] = 0; // Suppress output of global information
  //   mumps_struc[i].icntl[3] = 0; // gives no information for debugging
    mumps_struc[i].icntl[3] = 1; // error messages only
  //   mumps_struc[i].icntl[3] = 4; // gives lots of information for debugging
  //   opts->get("mumps_increase_working_space", mumps_struc[i].icntl[13], 20);
    mumps_struc[i].icntl[13] = 2000; // memory inrease (in %) needed for working space, default is 20
    mumps_struc[i].icntl[17] = 0; // provide matrix only on host
  // Don't use this: gives solution in some sructure convenient for the solver, not as desired by BOUT++       mumps_struc[i].icntl[20] = 1; // leave the solution distributed across the processors
  //   mumps_struc[i].icntl[22] = ; // Apparently should provide a value here significantly larger than infog[15] when running in parallel
  //   mumps_struc[i].icntl[27] = 2; // Force parallel analysis
  //   mumps_struc[i].icntl[28] = 0; // set to 1 to force pt-scotch and to 2 to force parmetis to be used for the ordering
  //   mumps_struc[i].cntl[0] = 0.01; // relative threshold for numerical pivoting. 0.0 is appropriate if matrix is diagonally dominant (there are no zero pivots). High values increase fill-in but lead to more accurate factorization. Default is 0.01
  // for (int i=0; i<40; i++) output<<i<<" "<<mumps_struc[i].icntl[i]<<endl;
  // exit(1);
    
    // Set diagonal indices for decayfactor part, cannot set indices for twist-shift part here as they may change with x-position
    if ( is_root[i] ) {
      mumps_struc[i].n = n;
      // nz is the total number of non-zero elements in the matrix
      mumps_struc[i].nz = n*4+n; // 4 values for interpolation to rotated position, 1 for diagonal decayfactor part
      mumps_struc[i].irn = new MUMPS_INT[mumps_struc[i].nz]; // list of row indices of matrix entries
      mumps_struc[i].jcn = new MUMPS_INT[mumps_struc[i].nz]; // list of column indices of matrix entries
      mumps_struc[i].a = new BoutReal[mumps_struc[i].nz]; // the matrix entries
      mumps_struc[i].nrhs = 1; // number of right hand side vectors
      mumps_struc[i].lrhs = mumps_struc[i].n; // leading dimension of rhs (i.e. length of vector)
      for (int j=0; j<n; j++) {
	mumps_struc[i].irn[j] = j+1; // Indices for fortran arrays that start at 1
	mumps_struc[i].jcn[j] = j+1;
      }
    }
      
    mumps_struc[i].job = MUMPS_JOB_ALL;
  }
}


// const FieldPerp NonLocalParallelToroidalSolver::solve(const FieldPerp &deltaI, const FieldPerp &decayfactor) {
//   BoutReal* deltaI_data = *deltaI.getData();
//   BoutReal* decayfactor_data = *decayfactor.getData();
//   deltaI_data += mesh->xstart*mesh->ngz;
//   decayfactor_data += mesh->xstart*mesh->ngz;
//   for (int jx=mesh->xstart; jx<=mesh->xend; jx++) {
//     BoutReal shiftangle;
//     if (mesh->periodicY(jx,shiftangle)) {
//       solve(deltaI_data, decayfactor_data, jx, shiftangle);
//       for (int jz=0; jz<n; jz++)
// 	sol[jx][jz]=mumps_struc[jx-mesh->xstart].rhs[jz];
//     }
//     else {
//       for (int jz=0; jz<n; jz++)
// 	sol[jx][jz]=0.;
//     }
//     deltaI_data += mesh->ngz;
//     decayfactor_data += mesh->ngz;
//   }
//   
//   return sol;
// }

void NonLocalParallelToroidalSolver::solve(BoutReal* deltaI, BoutReal* decayfactor, int jx, BoutReal shiftangle) {
  
  int i = jx-mesh->xstart;
  
{ Timer timer("mumpssetup");

  if (is_root[i]) {
    
// { Timer timer("mumpssetup");
    int j = 0; // Step through the MUMPS arrays of row/column indices and matrix element values
    
    for (int jz=0; jz<n; jz++) {
      mumps_struc[i].a[j] = decayfactor[jz];
      j++;
    }
    
    // mesh->periodicY(jx,shiftangle) gives the total shift from the first y-point to the last y-point on the grid
    if (!reverse) shiftangle *= -1.; // When interpolating from lastY to firstY we need to shift in the opposite direction
    {
      BoutReal totalangle = n*mesh->coordinates()->dz;
      int ishift=floor(shiftangle/totalangle);
      shiftangle -= ishift*totalangle;
    }
    int indexshift = int(floor(shiftangle/mesh->coordinates()->dz));
    BoutReal deltazind = shiftangle/mesh->coordinates()->dz-indexshift; // position to interpolate to after newjz, in index space (i.e. 0<=deltazind<=1)
    // Use Lagrange interpolation polynomials
    BoutReal Lagrange_coefficient0 = deltazind*(deltazind-1.)*(deltazind-2.)/(-1.)/(-2.)/(-3.);
    BoutReal Lagrange_coefficient1 = (deltazind+1.)*(deltazind-1.)*(deltazind-2.)/(1.)/(-1.)/(-2.);
    BoutReal Lagrange_coefficient2 = (deltazind+1.)*deltazind*(deltazind-2.)/(2.)/(1.)/(-1.);
    BoutReal Lagrange_coefficient3 = (deltazind+1.)*deltazind*(deltazind-1.)/(3.)/(2.)/(1.);
    for (int jz=0; jz<n; jz++) {
      int newjz = jz+indexshift;
      if (newjz<0) newjz+=n;
      else if (newjz>=n) newjz-=n;
      int newjzm = jz+indexshift-1;
      if (newjzm<0) newjzm+=n;
      else if (newjzm>=n) newjzm-=n;
      int newjzp = jz+indexshift+1;
      if (newjzp<0) newjzp+=n;
      else if (newjzp>=n) newjzp-=n;
      int newjzpp = jz+indexshift+2;
      if (newjzpp<0) newjzpp+=n;
      else if (newjzpp>=n) newjzpp-=n;
      // Use Lagrange interpolation polynomials
      mumps_struc[i].irn[j]=jz+1; // N.B. irn and icn store fortran array indices (start at 1)
      mumps_struc[i].jcn[j]=newjzm+1;
      mumps_struc[i].a[j]=Lagrange_coefficient0;
      j++;
      mumps_struc[i].irn[j]=jz+1;
      mumps_struc[i].jcn[j]=newjz+1;
      mumps_struc[i].a[j]=Lagrange_coefficient1;
      j++;
      mumps_struc[i].irn[j]=jz+1;
      mumps_struc[i].jcn[j]=newjzp+1;
      mumps_struc[i].a[j]=Lagrange_coefficient2;
      j++;
      mumps_struc[i].irn[j]=jz+1;
      mumps_struc[i].jcn[j]=newjzpp+1;
      mumps_struc[i].a[j]=Lagrange_coefficient3;
      j++;
    }
    
    #ifdef CHECK
      if (j!=mumps_struc[i].nz) throw BoutException("NonLocalParallelToroidalSolver: wrong number of matrix elements entered");
    #endif
      
    mumps_struc[i].rhs = deltaI;
  }
  
}
{ Timer timer("mumpssolve");
// output<<"irn	jcn	a"<<endl;
// for (int k=0; k<mumps_struc[i].nz;k++) output<<mumps_struc[i].irn[k]<<"	"<<mumps_struc[i].jcn[k]<<"	"<<mumps_struc[i].a[k]<<endl;
// output<<endl<<"rhs"<<endl;
// for (int k=0; k<mumps_struc[i].n;k++)output<<mumps_struc[i].rhs[k]<<endl;
// Solve the system
  dmumps_c( &mumps_struc[i] );
}
  if (mumps_struc[i].job==MUMPS_JOB_ALL) mumps_struc[i].job=MUMPS_JOB_BOTH;
//   }

}
