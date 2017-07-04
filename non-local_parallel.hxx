/**************************************************************************
 * Calculate non-local electron closures
 *
 **************************************************************************
 * Copyright 2012 J.T.Omotani
 *
 * Contact: John Omotani, john.omotani@york.ac.uk
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

#ifndef __NONLOCALPARALLEL_H__
#define __NONLOCALPARALLEL_H__

// #define NOEDGETERMS

#include <globals.hxx>

#include <sstream>
#include <cmath>

#include <bout_types.hxx>
#include <options.hxx>
#include <boutexception.hxx>
#include <utils.hxx>
#include <bout/sys/timer.hxx>
#include <output.hxx>
#include <msg_stack.hxx>
#include <bout/constants.hxx>
#include <fieldperp.hxx>
#include <field3d.hxx>
#include <field2d.hxx>
#include <stencils.hxx>
#include <interpolation.hxx>

#include "cubic_spline_local.hxx"
#include "non-local_parallel_integration.hxx"
#include "non-local_parallel_toroidal_solver.hxx"

class NonLocalParallel {
public:
  NonLocalParallel(const BoutReal &pass_elementary_charge, const BoutReal &pass_electron_mass/*, const BoutReal &pass_ion_mass*/, const BoutReal &pass_epsilon_0, const BoutReal &pass_logLambda, const bool pass_fluxes_location_is_ylow=false/*, const BoutReal &pass_gamma_factor=5./3.*/, Options* options = NULL);
  ~NonLocalParallel();
  int moments_number;
  Field3D electron_heat_flux;
  Field3D electron_viscosity;
  Field3D electron_friction;
  void set_density_normalization(BoutReal const pass_density_normalization) { density_normalization = pass_density_normalization; normalized_density = true; };
  void set_n_electron(const Field3D &pass_n_electron) { n_electron = pass_n_electron; };
  void set_T_electron(const Field3D &pass_T_electron) { T_electron = pass_T_electron; };
  void set_V_electron(const Field3D &pass_V_electron) { V_electron = pass_V_electron; };
  void set_j_parallel(const Field3D &pass_j_parallel) { j_parallel = pass_j_parallel; };
  void set_maxwellian_source(const Field3D &pass_source_tau, const Field3D &pass_source_strength) { source_tau = pass_source_tau; source_strength = pass_source_strength; };
  void set_V_ion(const Field3D &pass_V_ion) { V_ion = pass_V_ion; };
  void set_potential(const Field3D &pass_potential) { potential = pass_potential; };
  void set_heat_flux_boundary_condition(const Field3D &pass_heat_flux_boundary_condition) { heat_flux_boundary_condition = pass_heat_flux_boundary_condition; };
  void set_viscosity_boundary_condition(const Field3D &pass_viscosity_boundary_condition) { viscosity_boundary_condition = pass_viscosity_boundary_condition; };
  void calculate_nonlocal_closures();
  void set_boundary_gradients();
  void set_neumann_boundary_conditions();
  
  // Bit of a hack: should put the method somewhere more sensible (like in boutmesh.cxx)
  void y_broadcast(BoutReal* input_buffer, const int &size, const int &root_processor, const int &jx); // NB Assumes that the mesh is BoutMesh
//   void y_boundary_broadcast(BoutReal* input_buffer, const int &size, const int &root_processor); // NB Assumes that the mesh is BoutMesh
//   void rms_over_y(const Field3D &input_field, FieldPerp &output_field);
//   void mean_over_y(const Field3D &input_field, FieldPerp &output_field, int exclude_edgecells=0);
  BoutReal interp_to_point_YLOW(const Field3D &input, bindex &position);
  Field3D lambdaC_inverse;	// inverse collision length, i.e. 1/lambdaC
  
private:
  Coordinates *coord; ///< Metric tensor. Set in constructor
  
  Field3D n_electron;
  Field3D T_electron;
  Field3D V_electron;
  Field3D j_parallel;
  Field3D V_ion; // Used for the maxwellian_source_drives_oddp_terms calculation
  Field3D potential; // in same units as T_electron
  Field3D heat_flux_boundary_condition;
  Field3D viscosity_boundary_condition;
  Field3D source_tau; // Used to add kinetic contributions from sources, tau=T_source/T_electron
  Field3D source_strength; // Used to add kinetic contributions from sources, source_strength=particles/lengthunit^3/timeunit
  int number_of_source_drives;
  Field3D * source_driveterms; // Used to add kinetic contributions from sources
  void calculate_source_driveterms(); // Used to add kinetic contributions from sources
  bool gradT_drive;
  bool gradV_drive;
  bool VeminusVi_drive;
  int number_of_drives;
  bool calculate_heatflux;
  bool calculate_viscosity;
  bool calculate_friction;
  bool yperiodic;
  bool bc_heatflux;
  bool bc_viscosity;
  bool bc_apply_reflecting_sheath;
  bool maxwellian_source_drives;
  bool normalized_density;
  BoutReal density_normalization;
//   Field3D lambdaC_inverse;	// inverse collision length, i.e. 1/lambdaC
  Field3D increasing_dimensionless_length;	// d/dl(increasing_dimensionless_length)=1/lambdaC
  Field3D decreasing_dimensionless_length;
  Field3D dimensionless_length_deltas_above; // change in dimensionless_length between jy and jy+1, used in calculating integrals. (More accurate than first adding up and then taking differences)
  Field3D dimensionless_length_deltas_below;
  FieldPerp total_dimensionless_length;
  NonLocalParallelIntegration integration;
  int number_of_negative_eigenvalues;
  BoutReal * eigenvalues;
  BoutReal * exp_total_dimensionless_length_over_eigenvalue;
  BoutReal elementary_charge;
  BoutReal electron_mass;
//   BoutReal ion_mass;
  BoutReal epsilon_0;
  BoutReal logLambda;
  int boundary_gradient_smoothing_length;
  BoutReal boundary_condition_smoothing_range;
  CubicSpline cubic_spline_inverse_lambdaC;
  BoutReal * interp_coefficients;
  Field3D gradT_driveterm;	// drive from Maxwellian moments for heat flux calculation. g^(1,1)/T^1.5 as a function of l.
				      // NB slightly bad notation: drive_term will INCLUDE grad(T) for cell-centred version but EXCLUDE grad(T) (i.e. be only the cell-centred prefactor) for the ylow staggered version, so that grad(T) can be interpolated separately
  Field3D gradT_electron;	// only used for staggered-grids case
  Field3D gradV_driveterm;
  Field3D gradV_electron;
  Field3D VeminusVi_driveterm;
  CubicSpline * cubic_splines_driveterms_centre;
  CubicSpline * cubic_splines_driveterms_ylow;
  BoutReal heatflux_zerocoeff;
  BoutReal * heatflux_coefficients_below;
  BoutReal * heatflux_coefficients_above;
  BoutReal W11_dot_W20;
  BoutReal * WinverseB_11;
  BoutReal W11_dot_W11;
  BoutReal * heatflux_sol_transients_factors;
  FieldPerp pass_interim_upper_boundary_n11;
  FieldPerp upper_boundary_condition_n11;
  BoutReal * viscosity_coefficients;
  BoutReal W20_dot_W11;
  BoutReal * WinverseB_20;
  BoutReal W20_dot_W20;
  BoutReal * viscosity_sol_transients_factors;
  FieldPerp pass_interim_upper_boundary_n20;
  FieldPerp upper_boundary_condition_n20;
  BoutReal friction_zerocoeff;
  BoutReal * friction_coefficients_below;
  BoutReal * friction_coefficients_above;
  BoutReal * C10_1k_dot_W1k_B_times_WinverseB_11;
  BoutReal * C10_1k_dot_W1k_B_times_WinverseB_20;
  BoutReal gradT_zerocoeff;
  BoutReal VeminusVi_zerocoeff;
  BoutReal * driveterm_coefficients_below;
  BoutReal * driveterm_coefficients_above;
  NonLocalParallelToroidalSolver* toroidal_solver_forward;
  NonLocalParallelToroidalSolver* toroidal_solver_reverse;
  int transients_data_size;
  int transients_row_size;
  BoutReal* lower_transients_data;
  BoutReal* upper_transients_data;
  int NONLOCAL_PARALLEL_TAGBASE;
  MPI_Request broadcast_request;
  MPI_Comm* ycomms;
//   MPI_Comm comm_yprocs_minusone;
  bindex* position;
  bool fluxes_location_is_ylow;
  bool* is_lower_boundary;
  bool* is_upper_boundary;
  bool has_lower_boundary;
  bool has_upper_boundary;
  int* ycomms_rank;
  int* ycomms_size;
  int nx_core; // Number of x-gridpoints in core region on this processor
  int nx_sol; // Number of x-gridpoints in sol region on this processor
  int nz; // mesh->ngz-1
  void calculate_nonlocal_closures_cell_centre();
  void calculate_nonlocal_closures_cell_ylow();
  #ifdef CHECK
    bool calculated_before_setting_bcs;
  #endif
};

#endif // __NONLOCALPARALLEL_H__
