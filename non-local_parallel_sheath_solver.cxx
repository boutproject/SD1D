/*!
 * \file non-local_parallel_sheath_solver.cxx
 *
 * \brief Solver to impose sheath boundary conditions in non-local electron closures
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

#include "non-local_parallel_sheath_solver.hxx"

PetscErrorCode sheathsolver_mult_apply(Mat A, Vec vec_in, Vec vec_out) {
  SheathSolver* ctx;
  PetscErrorCode ierr = MatShellGetContext(A, &ctx);CHKERRQ(ierr);
  PetscFunctionReturn(ctx->mult(vec_in, vec_out));
}

SheathSolver::SheathSolver(const int pass_n, const int pass_n_perp, BoutReal* pass_eigenvalues,
			   const int moments_number_Hermite, const int moments_number_Laguerre,
			   BoutReal* const pass_deltaz, BoutReal* const pass_n_lower, BoutReal* const pass_n_upper,
			   BoutReal* const pass_n_plus, BoutReal* const pass_n_minus,
			   BoutReal* const pass_n_plus_integral, BoutReal* const pass_n_minus_integral,
			   BoutReal* const pass_n_zero_lower, BoutReal* const pass_n_zero_upper,
			   BoutReal* const pass_nS_lower, BoutReal* const pass_nS_upper,
			   BoutReal* const pass_ddsslower_nS_lower, BoutReal* const pass_ddssupper_nS_lower,
			   BoutReal* const pass_ddsslower_nS_upper, BoutReal* const pass_ddssupper_nS_upper,
			   BoutReal* const pass_ss_lower, BoutReal* const pass_ss_upper,
			   Options* opt) {
  
  n = pass_n;
  n_perp = pass_n_perp;
  eigenvalues = pass_eigenvalues;
  deltaz = pass_deltaz;
  n_lower = pass_n_lower;
  n_upper = pass_n_upper;
  n_plus = pass_n_plus;
  n_minus = pass_n_minus;
  n_plus_integral = pass_n_plus_integral;
  n_minus_integral = pass_n_minus_integral;
  n_zero_lower = pass_n_zero_lower;
  n_zero_upper = pass_n_zero_upper;
  nS_lower = pass_nS_lower;
  nS_upper = pass_nS_upper;
  ddsslower_nS_lower = pass_ddsslower_nS_lower;
  ddssupper_nS_lower = pass_ddssupper_nS_lower;
  ddsslower_nS_upper = pass_ddsslower_nS_upper;
  ddssupper_nS_upper = pass_ddssupper_nS_upper;
  ss_lower = pass_ss_lower;
  ss_upper = pass_ss_upper;
  
  OPTION(opt, approximate_ss_derivs, true); // If true drop terms in the ddss_lower and ddss_upper calculations that require another inversion to calculate
  
  nbuffer = new BoutReal[n];
  ddsslower_nminus = new BoutReal[n];
  ddssupper_nminus = new BoutReal[n];
  ddssupper_nplus = new BoutReal[n];
  ddsslower_nplus = new BoutReal[n];
  n_minus_previous = new BoutReal[n_perp*n];
  for (int i=0; i<n_perp*n; i++)
    n_minus_previous[i] = 0.;
  rhs = new BoutReal[n];
  
  ematrix = new Ematrix(n, moments_number_Hermite, moments_number_Laguerre);
  
  exp_deltaz_over_eigenvalues = new BoutReal[n];
  
  matrixsolver = new SerialMatrixSolve(n, this, &sheathsolver_mult_apply, opt); // The context pointer is the pointer to this SheathSolver, the static member staticmult(A, vec_in, vec_out) calls ctx->mult(A, vec_in, vec_out)
  
}

SheathSolver::~SheathSolver() {
  delete [] nbuffer;
  delete [] ddsslower_nminus;
  delete [] ddssupper_nminus;
  delete [] ddssupper_nplus;
  delete [] ddsslower_nplus;
  delete [] n_minus_previous;
  delete [] n_minus_previous;
  delete [] rhs;
  delete ematrix;
  delete [] exp_deltaz_over_eigenvalues;
  delete matrixsolver;
}

PetscErrorCode SheathSolver::solve() {
  
// output<<"n-prev ";for (int i=0; i<n; i++)output<<n_minus_previous[0*n+i]<<" ";output<<endl;
  PetscErrorCode ierr;
  
  for (int i=0; i<n_perp; i++) {
    ematrix->calculate_E(ss_lower[i], ss_upper[i]);
    ierr = solve_one(deltaz[i], n_lower[i], n_upper[i], n_plus+i*n, n_minus+i*n, n_plus_integral+i*n, n_minus_integral+i*n, n_minus_previous+i*n, n_zero_lower[i], n_zero_upper[i],
		     nS_lower[i], nS_upper[i], ddsslower_nS_lower[i], ddssupper_nS_lower[i], ddsslower_nS_upper[i], ddssupper_nS_upper[i]);CHKERRQ(ierr);
  }
// output<<"in solve() "<<nS_lower[0]<<endl;
  
  PetscFunctionReturn(0);
}

PetscErrorCode SheathSolver::solve_one(const BoutReal &single_deltaz, const BoutReal single_n_lower, const BoutReal single_n_upper, BoutReal* single_n_plus, BoutReal* single_n_minus,
				       BoutReal* single_n_plus_integral, BoutReal* single_n_minus_integral, BoutReal* single_n_minus_previous,
				       BoutReal &single_n_zero_lower, BoutReal &single_n_zero_upper, BoutReal &single_nS_lower, BoutReal &single_nS_upper,
				       BoutReal &single_ddsslower_nS_lower, BoutReal &single_ddssupper_nS_lower, BoutReal &single_ddsslower_nS_upper, BoutReal &single_ddssupper_nS_upper) {
  
  for (int i=0; i<n; i++)
    exp_deltaz_over_eigenvalues[i] = exp(single_deltaz/eigenvalues[i]);
  
  for (int i=0; i<n; i++) {
    nbuffer[i] = ematrix->E0_upper[i]*single_n_upper;
    for (int j=0; j<n; j++)
      nbuffer[i] += ematrix->E_upper[i*n+j]*single_n_minus_integral[j];
    nbuffer[i] *= exp_deltaz_over_eigenvalues[i];
    nbuffer[i] += single_n_plus_integral[i];
  }
  
  for (int i=0; i<n; i++) {
    rhs[i] = ematrix->E0_lower[i]*single_n_lower;
    for (int j=0; j<n; j++)
      rhs[i] += ematrix->E_lower[i*n+j]*nbuffer[j];
  }
  
// output<<"solve 1"<<endl;
  PetscErrorCode ierr = matrixsolver->solve(rhs, single_n_minus_previous);CHKERRQ(ierr);
  
  // Copy the solution into n_minus
  for (int i=0; i<n; i++) {
    single_n_minus[i] = single_n_minus_previous[i];
  }
  
  // Compute n_plus using deltan_minus stored in n_minus before overwriting n_minus
  for (int i=0; i<n; i++) {
    single_n_plus[i] = ematrix->E0_upper[i]*single_n_upper;
    for (int j=0; j<n; j++)
      single_n_plus[i] += ematrix->E_upper[i*n+j]*(single_n_minus_integral[j]+exp_deltaz_over_eigenvalues[j]*single_n_minus[j]);
  }
  
  // Compute the zero eigenvalue and 1,0 moment (n*S_parallel) boundary values
  single_n_zero_lower = -ematrix->zero_eigenvalue_E0_lower*single_n_lower; // on the lower boundary the component of n_(b) is -n_zero
  single_n_zero_upper = ematrix->zero_eigenvalue_E0_upper*single_n_upper;
  single_nS_lower = -ematrix->nS_E0_lower*single_n_lower; // on the lower boundary the component of n_(b) is -nS
  single_nS_upper = ematrix->nS_E0_upper*single_n_upper;
  for (int i=0; i<n; i++) {
    single_n_zero_lower -= ematrix->zero_eigenvalue_E_lower[i]*( single_n_plus_integral[i] + exp_deltaz_over_eigenvalues[i]*single_n_plus[i] ); // on the lower boundary the component of n_(b) is -n_zero
    single_n_zero_upper += ematrix->zero_eigenvalue_E_upper[i]*( single_n_minus_integral[i] + exp_deltaz_over_eigenvalues[i]*single_n_minus[i]/*[=n_minus|lower-boundary]*/ );
    single_nS_lower -= ematrix->nS_E_lower[i]*( single_n_plus_integral[i] + exp_deltaz_over_eigenvalues[i]*single_n_plus[i] ); // on the lower boundary the component of n_(b) is -nS
    single_nS_upper += ematrix->nS_E_upper[i]*( single_n_minus_integral[i] + exp_deltaz_over_eigenvalues[i]*single_n_minus[i]/*[=n_minus|lower-boundary]*/ );
  }
  
  // Compute the derivatives with ss boundary values of the 1,0 moment, n*S_parallel
  if (approximate_ss_derivs) {
    // This calculates approximate derivatives of n_plus and n_minus as their derivatives on the rhs are neglected
    for (int i=0; i<n; i++) {
      rhs[i] = ematrix->E0_upper_prime[i]*single_n_upper;
      nbuffer[i] = ematrix->E0_lower_prime[i]*single_n_lower;
      for (int j=0; j<n; j++) {
	rhs[i] += ematrix->E_upper_prime[i*n+j]*(single_n_minus_integral[j]+exp_deltaz_over_eigenvalues[j]*single_n_minus[j]);
	nbuffer[i] += ematrix->E_lower_prime[i*n+j]*(single_n_plus_integral[j]+exp_deltaz_over_eigenvalues[j]*single_n_plus[j]);
      }
      rhs[i] *= exp_deltaz_over_eigenvalues[i];
      nbuffer[i] *= exp_deltaz_over_eigenvalues[i];
    }
    for (int i=0; i<n; i++){
      ddssupper_nminus[i] = ematrix->E_lower[i*n+0]*rhs[0];
      ddsslower_nplus[i] = ematrix->E_upper[i*n+0]*nbuffer[0];
      for (int j=1; j<n; j++) {
	ddssupper_nminus[i] += ematrix->E_lower[i*n+j]*rhs[j];
	ddsslower_nplus[i] += ematrix->E_upper[i*n+j]*nbuffer[j];
      }
    }
    // Now use these derivatives to calculate the other pair (which would be exact if the previous ones were not approximate)
    for (int i=0; i<n; i++) {
      ddsslower_nminus[i] = ematrix->E0_lower_prime[i]*single_n_lower;
      ddssupper_nplus[i] = ematrix->E0_upper_prime[i]*single_n_upper;
      for (int j=0; j<n; j++) {
	ddsslower_nminus[i] += ematrix->E_lower_prime[i*n+j]*(single_n_plus_integral[j]+exp_deltaz_over_eigenvalues[j]*single_n_plus[j]) + ematrix->E_lower[i*n+j]*exp_deltaz_over_eigenvalues[j]*ddsslower_nplus[j];
	ddssupper_nplus[i] += ematrix->E_upper_prime[i*n+j]*(single_n_minus_integral[j]+exp_deltaz_over_eigenvalues[j]*single_n_minus[j]) + ematrix->E_upper[i*n+j]*exp_deltaz_over_eigenvalues[j]*ddssupper_nminus[j];
      }
    }
  }
  else {
    for (int i=0; i<n; i++) {
      rhs[i] = ematrix->E0_lower_prime[i]*single_n_lower;
      nbuffer[i] = ematrix->E0_upper[i]*single_n_upper; // nbuffer is the first part of the calculation of rhs for the second solve (to find ddssupper_nminus)
      for (int j=0; j<n; j++) {
	rhs[i] += ematrix->E_lower_prime[i*n+j]*(single_n_plus_integral[i]+exp_deltaz_over_eigenvalues[j]*single_n_plus[j]);
	nbuffer[i] += ematrix->E_upper_prime[i*n+j]*(single_n_minus_integral[j]+exp_deltaz_over_eigenvalues[j]*single_n_minus[j]);
      }
      nbuffer[i] *= exp_deltaz_over_eigenvalues[i];
    } 
    
// output<<"solve 2"<<endl;
    ierr = matrixsolver->solve(rhs, ddsslower_nminus);CHKERRQ(ierr);
    
    for (int i=0; i<n; i++) {
      rhs[i] = ematrix->E_lower[i*n+0]*nbuffer[0];
// output<<0<<"; ";output<<rhs[i]<<" ";output<<endl; 
      for (int j=1; j<n; j++) {
	rhs[i] += ematrix->E_lower[i*n+j]*nbuffer[j];
// output<<j<<"; ";output<<rhs[i]<<" ";output<<endl; 
      }
    }
    
// output<<"solve 3"<<endl;
// for (int i=0; i<n; i++) output<<ematrix->E_lower[0*n+i]<<" ";output<<endl;
    ierr = matrixsolver->solve(rhs, ddssupper_nminus);CHKERRQ(ierr);
    
    // Now use these derivatives to calculate the other pair
    for (int i=0; i<n; i++) {
      ddsslower_nplus[i] = 0.;
      ddssupper_nplus[i] = ematrix->E0_upper_prime[i]*single_n_upper;
      for (int j=0; j<n; j++) {
	ddsslower_nplus[i] += ematrix->E_upper[i*n+j]*exp_deltaz_over_eigenvalues[j]*ddsslower_nminus[j];
	ddssupper_nplus[i] += ematrix->E_upper_prime[i*n+j]*(single_n_minus_integral[j]+exp_deltaz_over_eigenvalues[j]*single_n_minus[j])
			      + ematrix->E_upper[i*n+j]*exp_deltaz_over_eigenvalues[j]*ddssupper_nminus[j];
      }
    }
  }
  single_ddsslower_nS_lower = -ematrix->nS_E0_lower_prime*single_n_lower; // on the lower boundary the component of n_(b) is -nS
  single_ddssupper_nS_lower = 0.;
  single_ddssupper_nS_upper = ematrix->nS_E0_upper_prime*single_n_upper;
  single_ddsslower_nS_upper = 0.;
// output<<"negative here? "<<ematrix->nS_E0_upper_prime<<" "<<single_ddssupper_nS_upper<<endl;
  for (int i=0; i<n; i++) {
    single_ddsslower_nS_lower -= ematrix->nS_E_lower_prime[i]*(single_n_plus_integral[i]+exp_deltaz_over_eigenvalues[i]*single_n_plus[i])
				  + ematrix->nS_E_lower[i]*exp_deltaz_over_eigenvalues[i]*ddsslower_nplus[i]; // on the lower boundary the component of n_(b) is -nS
    single_ddssupper_nS_lower -= ematrix->nS_E_lower[i]*exp_deltaz_over_eigenvalues[i]*ddssupper_nplus[i]; // on the lower boundary the component of n_(b) is -nS
    single_ddssupper_nS_upper += ematrix->nS_E_upper_prime[i]*(single_n_minus_integral[i]+exp_deltaz_over_eigenvalues[i]*single_n_minus[i])
				  + ematrix->nS_E_upper[i]*exp_deltaz_over_eigenvalues[i]*ddssupper_nminus[i];
    single_ddsslower_nS_upper += ematrix->nS_E_upper[i]*exp_deltaz_over_eigenvalues[i]*ddsslower_nminus[i];
  }
// output<<"negative here? "<<single_ddssupper_nS_upper<<endl;
  
  PetscFunctionReturn(0);
}

PetscErrorCode SheathSolver::mult(Vec vec_in, Vec vec_out) {
  
  PetscErrorCode ierr;
  BoutReal* data_in;
  BoutReal* data_out;
  ierr = VecGetArray(vec_in, &data_in);CHKERRQ(ierr);
  ierr = VecGetArray(vec_out, &data_out);CHKERRQ(ierr);
  
// output<<endl<<"a"<<endl;for (int i=0; i<n; i++)output<<data_in[i]<<",";
  for (int i=0; i<n; i++)
    nbuffer[i] = exp_deltaz_over_eigenvalues[i]*data_in[i];
  
// output<<endl<<"b"<<endl;for (int i=0; i<n; i++)output<<nbuffer[i]<<",";
  for (int i=0; i<n; i++) {
    data_out[i] = ematrix->E_upper[i*n+0]*nbuffer[0];
    for (int j=1; j<n; j++)
      data_out[i] += ematrix->E_upper[i*n+j]*nbuffer[j];
  }
  
// output<<endl<<"c"<<endl;for (int i=0; i<n; i++)output<<data_out[i]<<",";
  for (int i=0; i<n; i++)
    nbuffer[i] = exp_deltaz_over_eigenvalues[i]*data_out[i];
  
// output<<endl<<"d"<<endl;for (int i=0; i<n; i++)output<<nbuffer[i]<<",";
  for (int i=0; i<n; i++) {
    data_out[i] = data_in[i];
    for (int j=0; j<n; j++)
      data_out[i] -= ematrix->E_lower[i*n+j]*nbuffer[j];
  }
// output<<endl<<"e"<<endl;for (int i=0; i<n; i++)output<<data_out[i]<<",";
  
  ierr = VecRestoreArray(vec_in, &data_in);CHKERRQ(ierr);
  ierr = VecRestoreArray(vec_out, &data_out);CHKERRQ(ierr);
  
  
  PetscFunctionReturn(0);
}

Ematrix::Ematrix(const int pass_n, const int moments_number_Hermite, const int moments_number_Laguerre) {
  
  n = pass_n;
  
  interpolation_index_lower = 0;
  interpolation_index_upper = 0;
  
  std::stringstream infilename;
  infilename<<"nonlocal_coefficients/Ematrix"<<moments_number_Hermite<<","<<moments_number_Laguerre<<".0";
  std::ifstream infile( infilename.str().c_str() );
  if (!infile.is_open())
    throw BoutException("Could not open Ematrix file");
  infile >> nfiles;
  
  E_lower = new BoutReal[n*n];
  E_upper = new BoutReal[n*n];
  E0_lower = new BoutReal[n];
  E0_upper = new BoutReal[n];
  zero_eigenvalue_E_lower = new BoutReal[n];
  zero_eigenvalue_E_upper = new BoutReal[n];
  nS_E_lower = new BoutReal[n];
  nS_E_upper = new BoutReal[n];
  ss_min_array = new BoutReal[nfiles];
  ss_max_array = new BoutReal[nfiles];
  E_array = new BoutReal[nfiles*n*n];
  Eprime_array = new BoutReal[nfiles*n*n];
  E0_array = new BoutReal[nfiles*n];
  E0prime_array = new BoutReal[nfiles*n];
  zero_eigenvalue_E_array = new BoutReal[nfiles*n];
  zero_eigenvalue_Eprime_array = new BoutReal[nfiles*n];
  zero_eigenvalue_E0_array = new BoutReal[nfiles];
  zero_eigenvalue_E0prime_array = new BoutReal[nfiles];
  nS_E_array = new BoutReal[nfiles*n];
  nS_Eprime_array = new BoutReal[nfiles*n];
  nS_E0_array = new BoutReal[nfiles];
  nS_E0prime_array = new BoutReal[nfiles];
  
  for (int i=0; i<nfiles; i++) {
    if (i>0) {
      infilename.str(""); // Need to clear the old stuff
      infilename.clear(); // from infilename
      infilename<<"nonlocal_coefficients/Ematrix"<<moments_number_Hermite<<","<<moments_number_Laguerre<<"."<<i;
      infile.open( infilename.str().c_str() );
      if (!infile.is_open())
	throw BoutException("Could not open Ematrix file");
    }
    
    if (infile.eof())
      throw BoutException("reached end of Ematrix file unexpectedly");
    infile>>ss_min_array[i];
    if (infile.eof())
      throw BoutException("reached end of Ematrix file unexpectedly");
    infile>>ss_max_array[i];
    for (int j=0; j<n*n; j++) {
      if (infile.eof())
	throw BoutException("reached end of Ematrix file unexpectedly");
      infile>>E_array[i*n*n+j];
    }
    for (int j=0; j<n*n; j++) {
      if (infile.eof())
	throw BoutException("reached end of Ematrix file unexpectedly");
      infile>>Eprime_array[i*n*n+j];
    }
    for (int j=0; j<n; j++) {
      if (infile.eof())
	throw BoutException("reached end of Ematrix file unexpectedly");
      infile>>E0_array[i*n+j];
    }
    for (int j=0; j<n; j++) {
      if (infile.eof())
	throw BoutException("reached end of Ematrix file unexpectedly");
      infile>>E0prime_array[i*n+j];
    }
    for (int j=0; j<n; j++) {
      if (infile.eof())
	throw BoutException("reached end of Ematrix file unexpectedly");
      infile>>zero_eigenvalue_E_array[i*n+j];
    }
    for (int j=0; j<n; j++) {
      if (infile.eof())
	throw BoutException("reached end of Ematrix file unexpectedly");
      infile>>zero_eigenvalue_Eprime_array[i*n+j];
    }
    if (infile.eof())
      throw BoutException("reached end of Ematrix file unexpectedly");
    infile>>zero_eigenvalue_E0_array[i];
    if (infile.eof())
      throw BoutException("reached end of Ematrix file unexpectedly");
    infile>>zero_eigenvalue_E0prime_array[i];
    for (int j=0; j<n; j++) {
      if (infile.eof())
	throw BoutException("reached end of Ematrix file unexpectedly");
      infile>>nS_E_array[i*n+j];
    }
    for (int j=0; j<n; j++) {
      if (infile.eof())
	throw BoutException("reached end of Ematrix file unexpectedly");
      infile>>nS_Eprime_array[i*n+j];
    }
    if (infile.eof())
      throw BoutException("reached end of Ematrix file unexpectedly");
    infile>>nS_E0_array[i];
    if (infile.eof())
      throw BoutException("reached end of Ematrix file unexpectedly");
    infile>>nS_E0prime_array[i];
    infile.close();
  }
}

Ematrix::~Ematrix() {
  delete [] E_lower;
  delete [] E_upper;
  delete [] E0_lower;
  delete [] E0_upper;
  delete [] zero_eigenvalue_E_lower;
  delete [] zero_eigenvalue_E_upper;
  delete [] nS_E_lower;
  delete [] nS_E_upper;
  delete [] E_array;
  delete [] Eprime_array;
  delete [] E0_array;
  delete [] E0prime_array;
  delete [] zero_eigenvalue_E_array;
  delete [] zero_eigenvalue_Eprime_array;
  delete [] zero_eigenvalue_E0_array;
  delete [] zero_eigenvalue_E0prime_array;
  delete [] nS_E_array;
  delete [] nS_Eprime_array;
  delete [] nS_E0_array;
  delete [] nS_E0prime_array;
}

void Ematrix::calculate_E(const BoutReal &ss_lower, const BoutReal &ss_upper) {
  
// output<<endl<<"l="<<ss_lower<<" u="<<ss_upper<<endl;
// output<<ss_min_array[0]<<" "<<ss_max_array[nfiles-1]<<endl;
  if (!(ss_lower>ss_min_array[0]))
    throw BoutException("ss_lower too small");
  else if (!(ss_lower<ss_max_array[nfiles-1])) {
//     throw BoutException("ss_lower too big");
    output<<"WARNING: ss_lower too big"<<endl;
    interpolation_index_lower = nfiles-1;
  }
  else {
    while (ss_min_array[interpolation_index_lower]>ss_lower)
      interpolation_index_lower--;
    while (ss_max_array[interpolation_index_lower]<ss_lower)
      interpolation_index_lower++;
  }
  if (!(ss_upper>ss_min_array[0]))
    throw BoutException("ss_upper too small");
  else if (!(ss_upper<ss_max_array[nfiles-1])) {
//     throw BoutException("ss_upper too big");
    output<<"WARNING: ss_upper too big"<<endl;
    interpolation_index_upper = nfiles-1;
  }
  else {
    while (ss_min_array[interpolation_index_upper]>ss_upper)
      interpolation_index_upper--;
    while (ss_max_array[interpolation_index_upper]<ss_upper)
      interpolation_index_upper++;
  }
  
  BoutReal delta_ss_lower = ss_lower - ss_min_array[interpolation_index_lower];
//   BoutReal ss_width_lower = ss_max_array[interpolation_index_lower] - ss_min_array[interpolation_index_upper];
  BoutReal delta_ss_upper = ss_upper - ss_min_array[interpolation_index_upper];
//   BoutReal ss_width_upper = ss_max_array[interpolation_index_upper] - ss_min_array[interpolation_index_upper];
  
  for (int i=0; i<n*n; i++) {
    E_lower[i] = E_array[interpolation_index_lower*n*n+i] + delta_ss_lower*Eprime_array[interpolation_index_lower*n*n+i];
    E_upper[i] = E_array[interpolation_index_upper*n*n+i] + delta_ss_upper*Eprime_array[interpolation_index_upper*n*n+i];
  }
  for (int i=0; i<n; i++) {
    E0_lower[i] = E0_array[interpolation_index_lower*n+i] + delta_ss_lower*E0prime_array[interpolation_index_lower*n+i];
    E0_upper[i] = E0_array[interpolation_index_upper*n+i] + delta_ss_upper*E0prime_array[interpolation_index_upper*n+i];
    zero_eigenvalue_E_lower[i] = zero_eigenvalue_E_array[interpolation_index_lower*n+i] + delta_ss_lower*zero_eigenvalue_Eprime_array[interpolation_index_lower*n+i];
    zero_eigenvalue_E_upper[i] = zero_eigenvalue_E_array[interpolation_index_upper*n+i] + delta_ss_upper*zero_eigenvalue_Eprime_array[interpolation_index_upper*n+i];
    nS_E_lower[i] = nS_E_array[interpolation_index_lower*n+i] + delta_ss_lower*nS_Eprime_array[interpolation_index_lower*n+i];
    nS_E_upper[i] = nS_E_array[interpolation_index_upper*n+i] + delta_ss_upper*nS_Eprime_array[interpolation_index_upper*n+i];
  }
  zero_eigenvalue_E0_lower = zero_eigenvalue_E0_array[interpolation_index_lower] + delta_ss_lower*zero_eigenvalue_E0prime_array[interpolation_index_lower];
  zero_eigenvalue_E0_upper = zero_eigenvalue_E0_array[interpolation_index_upper] + delta_ss_upper*zero_eigenvalue_E0prime_array[interpolation_index_upper];
  nS_E0_lower = nS_E0_array[interpolation_index_lower] + delta_ss_lower*nS_E0prime_array[interpolation_index_lower];
  nS_E0_upper = nS_E0_array[interpolation_index_upper] + delta_ss_upper*nS_E0prime_array[interpolation_index_upper];
  
  E_lower_prime = &Eprime_array[interpolation_index_lower*n*n];
  E_upper_prime = &Eprime_array[interpolation_index_upper*n*n];
  E0_lower_prime = &E0prime_array[interpolation_index_lower*n];
  E0_upper_prime = &E0prime_array[interpolation_index_upper*n];
  nS_E_lower_prime = &nS_Eprime_array[interpolation_index_lower*n];
  nS_E_upper_prime = &nS_Eprime_array[interpolation_index_upper*n];
  nS_E0_lower_prime = nS_E0prime_array[interpolation_index_lower];
  nS_E0_upper_prime = nS_E0prime_array[interpolation_index_upper];
  
//   // E_array holds the values of E at the values given by ss_min_array and Eprime_array the (2nd order central finite difference approximation to) the values of ddss_E, etc.
//   // Calculate the value at ss between ss_min and ss_max by a cubic spline (so that the interpolation is continuous and smooth (C1))
//   for (int i=0; i<n*n; i++) {
//     E_lower[i] = E_array[interpolation_index_lower*n*n+i] 
// 		  + delta_ss_lower * ( Eprime_array[interpolation_index_lower*n*n+i]
// 		  + delta_ss_lower * ( (3.*(E_array[(interpolation_index_lower+1)*n*n+i]-E_array[interpolation_index_lower*n*n+i])/ss_width_lower - Eprime_array[(interpolation_index_lower+1)*n*n+i] + 4.*Eprime_array[interpolation_index_lower*n*n+i])/ss_width_lower
// 		  + delta_ss_lower * ( (-2.*E_array[(interpolation_index_lower+1)*n*n+i]+E_array[interpolation_index_lower*n*n+i])/ss_width_lower + Eprime_array[(interpolation_index_lower+1)*n*n+i]-3.*Eprime_array[interpolation_index_lower*n*n+i] )/ss_width_lower/ss_width_lower ) );
//     E_upper[i] = E_array[interpolation_index_upper*n*n+i]
// 		  + delta_ss_upper * ( Eprime_array[interpolation_index_upper*n*n+i]
// 		  + delta_ss_upper * ( (3.*(E_array[(interpolation_index_upper+1)*n*n+i]-E_array[interpolation_index_upper*n*n+i])/ss_width_upper - Eprime_array[(interpolation_index_upper+1)*n*n+i] + 4.*Eprime_array[interpolation_index_upper*n*n+i])/ss_width_upper
// 		  + delta_ss_upper * ( (-2.*E_array[(interpolation_index_upper+1)*n*n+i]+E_array[interpolation_index_upper*n*n+i])/ss_width_upper + Eprime_array[(interpolation_index_upper+1)*n*n+i]-3.*Eprime_array[interpolation_index_upper*n*n+i] )/ss_width_upper/ss_width_upper ) );
//   }
//   for (int i=0; i<n; i++) {
//     E0_lower[i] = E0_array[interpolation_index_lower*n+i] 
// 		  + delta_ss_lower * ( E0prime_array[interpolation_index_lower*n+i]
// 		  + delta_ss_lower * ( (3.*(E0_array[(interpolation_index_lower+1)*n+i]-E0_array[interpolation_index_lower*n+i])/ss_width_lower - E0prime_array[(interpolation_index_lower+1)*n+i] + 4.*E0prime_array[interpolation_index_lower*n+i])/ss_width_lower
// 		  + delta_ss_lower * ( (-2.*E0_array[(interpolation_index_lower+1)*n+i]+E0_array[interpolation_index_lower*n+i])/ss_width_lower + E0prime_array[(interpolation_index_lower+1)*n+i]-3.*E0prime_array[interpolation_index_lower*n+i] )/ss_width_lower/ss_width_lower ) );
//     E0_upper[i] = E0_array[interpolation_index_upper*n+i]
// 		  + delta_ss_upper * ( E0prime_array[interpolation_index_upper*n+i]
// 		  + delta_ss_upper * ( (3.*(E0_array[(interpolation_index_upper+1)*n+i]-E0_array[interpolation_index_upper*n+i])/ss_width_upper - E0prime_array[(interpolation_index_upper+1)*n+i] + 4.*E0prime_array[interpolation_index_upper*n+i])/ss_width_upper
// 		  + delta_ss_upper * ( (-2.*E0_array[(interpolation_index_upper+1)*n+i]+E0_array[interpolation_index_upper*n+i])/ss_width_upper + E0prime_array[(interpolation_index_upper+1)*n+i]-3.*E0prime_array[interpolation_index_upper*n+i] )/ss_width_upper/ss_width_upper ) );
//     zero_eigenvalue_E_lower[i] = zero_eigenvalue_E_array[interpolation_index_lower*n+i] 
// 		  + delta_ss_lower * ( zero_eigenvalue_Eprime_array[interpolation_index_lower*n+i]
// 		  + delta_ss_lower * ( (3.*(zero_eigenvalue_E_array[(interpolation_index_lower+1)*n+i]-zero_eigenvalue_E_array[interpolation_index_lower*n+i])/ss_width_lower - zero_eigenvalue_Eprime_array[(interpolation_index_lower+1)*n+i] + 4.*zero_eigenvalue_Eprime_array[interpolation_index_lower*n+i])/ss_width_lower
// 		  + delta_ss_lower * ( (-2.*zero_eigenvalue_E_array[(interpolation_index_lower+1)*n+i]+zero_eigenvalue_E_array[interpolation_index_lower*n+i])/ss_width_lower + zero_eigenvalue_Eprime_array[(interpolation_index_lower+1)*n+i]-3.*zero_eigenvalue_Eprime_array[interpolation_index_lower*n+i] )/ss_width_lower/ss_width_lower ) );
//     zero_eigenvalue_E_upper[i] = zero_eigenvalue_E_array[interpolation_index_upper*n+i]
// 		  + delta_ss_upper * ( zero_eigenvalue_Eprime_array[interpolation_index_upper*n+i]
// 		  + delta_ss_upper * ( (3.*(zero_eigenvalue_E_array[(interpolation_index_upper+1)*n+i]-zero_eigenvalue_E_array[interpolation_index_upper*n+i])/ss_width_upper - zero_eigenvalue_Eprime_array[(interpolation_index_upper+1)*n+i] + 4.*zero_eigenvalue_Eprime_array[interpolation_index_upper*n+i])/ss_width_upper
// 		  + delta_ss_upper * ( (-2.*zero_eigenvalue_E_array[(interpolation_index_upper+1)*n+i]+zero_eigenvalue_E_array[interpolation_index_upper*n+i])/ss_width_upper + zero_eigenvalue_Eprime_array[(interpolation_index_upper+1)*n+i]-3.*zero_eigenvalue_Eprime_array[interpolation_index_upper*n+i] )/ss_width_upper/ss_width_upper ) );
//     nS_E_lower[i] = nS_E_array[interpolation_index_lower*n+i] 
// 		  + delta_ss_lower * ( nS_Eprime_array[interpolation_index_lower*n+i]
// 		  + delta_ss_lower * ( (3.*(nS_E_array[(interpolation_index_lower+1)*n+i]-nS_E_array[interpolation_index_lower*n+i])/ss_width_lower - nS_Eprime_array[(interpolation_index_lower+1)*n+i] + 4.*nS_Eprime_array[interpolation_index_lower*n+i])/ss_width_lower
// 		  + delta_ss_lower * ( (-2.*nS_E_array[(interpolation_index_lower+1)*n+i]+nS_E_array[interpolation_index_lower*n+i])/ss_width_lower + nS_Eprime_array[(interpolation_index_lower+1)*n+i]-3.*nS_Eprime_array[interpolation_index_lower*n+i] )/ss_width_lower/ss_width_lower ) );
//     nS_E_upper[i] = nS_E_array[interpolation_index_upper*n+i]
// 		  + delta_ss_upper * ( nS_Eprime_array[interpolation_index_upper*n+i]
// 		  + delta_ss_upper * ( (3.*(nS_E_array[(interpolation_index_upper+1)*n+i]-nS_E_array[interpolation_index_upper*n+i])/ss_width_upper - nS_Eprime_array[(interpolation_index_upper+1)*n+i] + 4.*nS_Eprime_array[interpolation_index_upper*n+i])/ss_width_upper
// 		  + delta_ss_upper * ( (-2.*nS_E_array[(interpolation_index_upper+1)*n+i]+nS_E_array[interpolation_index_upper*n+i])/ss_width_upper + nS_Eprime_array[(interpolation_index_upper+1)*n+i]-3.*nS_Eprime_array[interpolation_index_upper*n+i] )/ss_width_upper/ss_width_upper ) );
//   }
//   zero_eigenvalue_E0_lower = zero_eigenvalue_E0_array[interpolation_index_lower] 
// 		+ delta_ss_lower * ( zero_eigenvalue_E0prime_array[interpolation_index_lower]
// 		+ delta_ss_lower * ( (3.*(zero_eigenvalue_E0_array[interpolation_index_lower+1]-zero_eigenvalue_E0_array[interpolation_index_lower])/ss_width_lower - zero_eigenvalue_E0prime_array[interpolation_index_lower+1] + 4.*zero_eigenvalue_E0prime_array[interpolation_index_lower])/ss_width_lower
// 		+ delta_ss_lower * ( (-2.*zero_eigenvalue_E0_array[interpolation_index_lower+1]+zero_eigenvalue_E0_array[interpolation_index_lower])/ss_width_lower + zero_eigenvalue_E0prime_array[interpolation_index_lower+1]-3.*zero_eigenvalue_E0prime_array[interpolation_index_lower] )/ss_width_lower/ss_width_lower ) );
//   zero_eigenvalue_E0_upper = zero_eigenvalue_E0_array[interpolation_index_upper]
// 		+ delta_ss_upper * ( zero_eigenvalue_E0prime_array[interpolation_index_upper]
// 		+ delta_ss_upper * ( (3.*(zero_eigenvalue_E0_array[interpolation_index_upper+1]-zero_eigenvalue_E0_array[interpolation_index_upper])/ss_width_upper - zero_eigenvalue_E0prime_array[interpolation_index_upper+1] + 4.*zero_eigenvalue_E0prime_array[interpolation_index_upper])/ss_width_upper
// 		+ delta_ss_upper * ( (-2.*zero_eigenvalue_E0_array[interpolation_index_upper+1]+zero_eigenvalue_E0_array[interpolation_index_upper])/ss_width_upper + zero_eigenvalue_E0prime_array[interpolation_index_upper+1]-3.*zero_eigenvalue_E0prime_array[interpolation_index_upper] )/ss_width_upper/ss_width_upper ) );
//   nS_E0_lower = nS_E0_array[interpolation_index_lower] 
// 		+ delta_ss_lower * ( nS_E0prime_array[interpolation_index_lower]
// 		+ delta_ss_lower * ( (3.*(nS_E0_array[interpolation_index_lower+1]-nS_E0_array[interpolation_index_lower])/ss_width_lower - nS_E0prime_array[interpolation_index_lower+1] + 4.*nS_E0prime_array[interpolation_index_lower])/ss_width_lower
// 		+ delta_ss_lower * ( (-2.*nS_E0_array[interpolation_index_lower+1]+nS_E0_array[interpolation_index_lower])/ss_width_lower + nS_E0prime_array[interpolation_index_lower+1]-3.*nS_E0prime_array[interpolation_index_lower] )/ss_width_lower/ss_width_lower ) );
//   nS_E0_upper = nS_E0_array[interpolation_index_upper]
// 		+ delta_ss_upper * ( nS_E0prime_array[interpolation_index_upper]
// 		+ delta_ss_upper * ( (3.*(nS_E0_array[interpolation_index_upper+1]-nS_E0_array[interpolation_index_upper])/ss_width_upper - nS_E0prime_array[interpolation_index_upper+1] + 4.*nS_E0prime_array[interpolation_index_upper])/ss_width_upper
// 		+ delta_ss_upper * ( (-2.*nS_E0_array[interpolation_index_upper+1]+nS_E0_array[interpolation_index_upper])/ss_width_upper + nS_E0prime_array[interpolation_index_upper+1]-3.*nS_E0prime_array[interpolation_index_upper] )/ss_width_upper/ss_width_upper ) );
//   for (int i=0; i<n*n; i++) {
//     Eprime_lower[i] = ( Eprime_array[interpolation_index_lower*n*n+i]
// 		  + delta_ss_lower * 2.*( (3.*(E_array[(interpolation_index_lower+1)*n*n+i]-E_array[interpolation_index_lower*n*n+i])/ss_width_lower - Eprime_array[(interpolation_index_lower+1)*n*n+i] + 4.*Eprime_array[interpolation_index_lower*n*n+i])/ss_width_lower
// 		  + delta_ss_lower * 3.*( (-2.*E_array[(interpolation_index_lower+1)*n*n+i]+E_array[interpolation_index_lower*n*n+i])/ss_width_lower + Eprime_array[(interpolation_index_lower+1)*n*n+i]-3.*Eprime_array[interpolation_index_lower*n*n+i] )/ss_width_lower/ss_width_lower );
//     Eprime_upper[i] = ( Eprime_array[interpolation_index_upper*n*n+i]
// 		  + delta_ss_upper * 2.*( (3.*(E_array[(interpolation_index_upper+1)*n*n+i]-E_array[interpolation_index_upper*n*n+i])/ss_width_upper - Eprime_array[(interpolation_index_upper+1)*n*n+i] + 4.*Eprime_array[interpolation_index_upper*n*n+i])/ss_width_upper
// 		  + delta_ss_upper * 3.*( (-2.*E_array[(interpolation_index_upper+1)*n*n+i]+E_array[interpolation_index_upper*n*n+i])/ss_width_upper + Eprime_array[(interpolation_index_upper+1)*n*n+i]-3.*Eprime_array[interpolation_index_upper*n*n+i] )/ss_width_upper/ss_width_upper );
//   }
//   for (int i=0; i<n; i++) {
//     E0prime_lower[i] = ( E0prime_array[interpolation_index_lower*n+i]
// 		  + delta_ss_lower * 2.*( (3.*(E0_array[(interpolation_index_lower+1)*n+i]-E0_array[interpolation_index_lower*n+i])/ss_width_lower - E0prime_array[(interpolation_index_lower+1)*n+i] + 4.*E0prime_array[interpolation_index_lower*n+i])/ss_width_lower
// 		  + delta_ss_lower * 3.*( (-2.*E0_array[(interpolation_index_lower+1)*n+i]+E0_array[interpolation_index_lower*n+i])/ss_width_lower + E0prime_array[(interpolation_index_lower+1)*n+i]-3.*E0prime_array[interpolation_index_lower*n+i] )/ss_width_lower/ss_width_lower );
//     E0prime_upper[i] = ( E0prime_array[interpolation_index_upper*n+i]
// 		  + delta_ss_upper * 2.*( (3.*(E0_array[(interpolation_index_upper+1)*n+i]-E0_array[interpolation_index_upper*n+i])/ss_width_upper - E0prime_array[(interpolation_index_upper+1)*n+i] + 4.*E0prime_array[interpolation_index_upper*n+i])/ss_width_upper
// 		  + delta_ss_upper * 3.*( (-2.*E0_array[(interpolation_index_upper+1)*n+i]+E0_array[interpolation_index_upper*n+i])/ss_width_upper + E0prime_array[(interpolation_index_upper+1)*n+i]-3.*E0prime_array[interpolation_index_upper*n+i] )/ss_width_upper/ss_width_upper );
//     zero_eigenvalue_Eprime_lower[i] = ( zero_eigenvalue_Eprime_array[interpolation_index_lower*n+i]
// 		  + delta_ss_lower * 2.*( (3.*(zero_eigenvalue_E_array[(interpolation_index_lower+1)*n+i]-zero_eigenvalue_E_array[interpolation_index_lower*n+i])/ss_width_lower - zero_eigenvalue_Eprime_array[(interpolation_index_lower+1)*n+i] + 4.*zero_eigenvalue_Eprime_array[interpolation_index_lower*n+i])/ss_width_lower
// 		  + delta_ss_lower * 3.*( (-2.*zero_eigenvalue_E_array[(interpolation_index_lower+1)*n+i]+zero_eigenvalue_E_array[interpolation_index_lower*n+i])/ss_width_lower + zero_eigenvalue_Eprime_array[(interpolation_index_lower+1)*n+i]-3.*zero_eigenvalue_Eprime_array[interpolation_index_lower*n+i] )/ss_width_lower/ss_width_lower );
//     zero_eigenvalue_Eprime_upper[i] = ( zero_eigenvalue_Eprime_array[interpolation_index_upper*n+i]
// 		  + delta_ss_upper * 2.*( (3.*(zero_eigenvalue_E_array[(interpolation_index_upper+1)*n+i]-zero_eigenvalue_E_array[interpolation_index_upper*n+i])/ss_width_upper - zero_eigenvalue_Eprime_array[(interpolation_index_upper+1)*n+i] + 4.*zero_eigenvalue_Eprime_array[interpolation_index_upper*n+i])/ss_width_upper
// 		  + delta_ss_upper * 3.*( (-2.*zero_eigenvalue_E_array[(interpolation_index_upper+1)*n+i]+zero_eigenvalue_E_array[interpolation_index_upper*n+i])/ss_width_upper + zero_eigenvalue_Eprime_array[(interpolation_index_upper+1)*n+i]-3.*zero_eigenvalue_Eprime_array[interpolation_index_upper*n+i] )/ss_width_upper/ss_width_upper );
//     nS_Eprime_lower[i] = ( nS_Eprime_array[interpolation_index_lower*n+i]
// 		  + delta_ss_lower * 2.*( (3.*(nS_E_array[(interpolation_index_lower+1)*n+i]-nS_E_array[interpolation_index_lower*n+i])/ss_width_lower - nS_Eprime_array[(interpolation_index_lower+1)*n+i] + 4.*nS_Eprime_array[interpolation_index_lower*n+i])/ss_width_lower
// 		  + delta_ss_lower * 3.*( (-2.*nS_E_array[(interpolation_index_lower+1)*n+i]+nS_E_array[interpolation_index_lower*n+i])/ss_width_lower + nS_Eprime_array[(interpolation_index_lower+1)*n+i]-3.*nS_Eprime_array[interpolation_index_lower*n+i] )/ss_width_lower/ss_width_lower );
//     nS_Eprime_upper[i] = ( nS_Eprime_array[interpolation_index_upper*n+i]
// 		  + delta_ss_upper * 2.*( (3.*(nS_E_array[(interpolation_index_upper+1)*n+i]-nS_E_array[interpolation_index_upper*n+i])/ss_width_upper - nS_Eprime_array[(interpolation_index_upper+1)*n+i] + 4.*nS_Eprime_array[interpolation_index_upper*n+i])/ss_width_upper
// 		  + delta_ss_upper * 3.*( (-2.*nS_E_array[(interpolation_index_upper+1)*n+i]+nS_E_array[interpolation_index_upper*n+i])/ss_width_upper + nS_Eprime_array[(interpolation_index_upper+1)*n+i]-3.*nS_Eprime_array[interpolation_index_upper*n+i] )/ss_width_upper/ss_width_upper );
//   }
//   zero_eigenvalue_E0prime_lower = ( zero_eigenvalue_E0prime_array[interpolation_index_lower]
// 		+ delta_ss_lower * 2.*( (3.*(zero_eigenvalue_E0_array[interpolation_index_lower+1]-zero_eigenvalue_E0_array[interpolation_index_lower])/ss_width_lower - zero_eigenvalue_E0prime_array[interpolation_index_lower+1] + 4.*zero_eigenvalue_E0prime_array[interpolation_index_lower])/ss_width_lower
// 		+ delta_ss_lower * 3.*( (-2.*zero_eigenvalue_E0_array[interpolation_index_lower+1]+zero_eigenvalue_E0_array[interpolation_index_lower])/ss_width_lower + zero_eigenvalue_E0prime_array[interpolation_index_lower+1]-3.*zero_eigenvalue_E0prime_array[interpolation_index_lower] )/ss_width_lower/ss_width_lower );
//   zero_eigenvalue_E0prime_upper = ( zero_eigenvalue_E0prime_array[interpolation_index_upper]
// 		+ delta_ss_upper * 2.*( (3.*(zero_eigenvalue_E0_array[interpolation_index_upper+1]-zero_eigenvalue_E0_array[interpolation_index_upper])/ss_width_upper - zero_eigenvalue_E0prime_array[interpolation_index_upper+1] + 4.*zero_eigenvalue_E0prime_array[interpolation_index_upper])/ss_width_upper
// 		+ delta_ss_upper * 3.*( (-2.*zero_eigenvalue_E0_array[interpolation_index_upper+1]+zero_eigenvalue_E0_array[interpolation_index_upper])/ss_width_upper + zero_eigenvalue_E0prime_array[interpolation_index_upper+1]-3.*zero_eigenvalue_E0prime_array[interpolation_index_upper] )/ss_width_upper/ss_width_upper );
//   nS_E0prime_lower = ( nS_E0prime_array[interpolation_index_lower]
// 		+ delta_ss_lower * 2.*( (3.*(nS_E0_array[interpolation_index_lower+1]-nS_E0_array[interpolation_index_lower])/ss_width_lower - nS_E0prime_array[interpolation_index_lower+1] + 4.*nS_E0prime_array[interpolation_index_lower])/ss_width_lower
// 		+ delta_ss_lower * 3.*( (-2.*nS_E0_array[interpolation_index_lower+1]+nS_E0_array[interpolation_index_lower])/ss_width_lower + nS_E0prime_array[interpolation_index_lower+1]-3.*nS_E0prime_array[interpolation_index_lower] )/ss_width_lower/ss_width_lower );
//   nS_E0prime_upper = ( nS_E0prime_array[interpolation_index_upper]
// 		+ delta_ss_upper * 2.*( (3.*(nS_E0_array[interpolation_index_upper+1]-nS_E0_array[interpolation_index_upper])/ss_width_upper - nS_E0prime_array[interpolation_index_upper+1] + 4.*nS_E0prime_array[interpolation_index_upper])/ss_width_upper
// 		+ delta_ss_upper * 3.*( (-2.*nS_E0_array[interpolation_index_upper+1]+nS_E0_array[interpolation_index_upper])/ss_width_upper + nS_E0prime_array[interpolation_index_upper+1]-3.*nS_E0prime_array[interpolation_index_upper] )/ss_width_upper/ss_width_upper );
  
}
