/*!
 * \file non-local_parallel.cxx
 *
 * \brief Calculate non-local electron closures
 *
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
 */

#include "non-local_parallel.hxx"

#include <fstream>
#include <sstream>

void savetofile(Field3D &f) {
  std::ofstream file;
  std::stringstream filename;
  filename<<"checknew"<<mesh->getNXPE()<<","<<mesh->getNYPE()<<"_"<<mesh->getXProcIndex()<<","<<mesh->getYProcIndex();
  file.open(filename.str().c_str());
  for (int jx=0; jx<mesh->LocalNx; jx++)
    for (int jy=0; jy<mesh->LocalNy; jy++) {
      file<<jx<<","<<jy<<" ";
      for (int jz=0; jz<mesh->LocalNz; jz++)
	file<<f(jx,jy,jz)<<" ";
      file<<endl;
    }
  file.close();
  MPI_Barrier(BoutComm::get());
  exit(43);
}
void savetofile(FieldPerp &f) {
  std::ofstream file;
  std::stringstream filename;
  filename<<"checknew"<<mesh->getNXPE()<<","<<mesh->getNYPE()<<"_"<<mesh->getXProcIndex()<<","<<mesh->getYProcIndex();
  file.open(filename.str().c_str());
  for (int jx=0; jx<mesh->LocalNx; jx++) {
      file<<jx<<" ";
      for (int jz=0; jz<mesh->LocalNz; jz++)
	file<<f(jx,jz)<<" ";
      file<<endl;
  }
  file.close();
  MPI_Barrier(BoutComm::get());
  exit(43);
}

const FieldPerp exp(const FieldPerp &f) {
  FieldPerp result;
  result.allocate();
  for(auto &i : result) {
    result[i] = ::exp(f[i]);
  }
  return result;
}

// BoutReal NonLocalParallelgamma_factor = 3.;

/**********************************************************************************
 *                    HEAT FLUX INITIALISATION AND CREATION                       *
 **********************************************************************************/

NonLocalParallel::NonLocalParallel(const BoutReal &pass_elementary_charge, const BoutReal &pass_electron_mass/*, const BoutReal &pass_ion_mass*/, const BoutReal &pass_epsilon_0, const BoutReal &pass_logLambda, const bool pass_fluxes_location_is_ylow/*, const BoutReal &pass_gamma_factor*/, Options* options) {

  fluxes_location_is_ylow = pass_fluxes_location_is_ylow;
#ifdef CHECK
  calculated_before_setting_bcs=false;
#endif
  
  // Get the metric tensor from Mesh
  coord = mesh->coordinates();

  normalized_density = false;
  is_lower_boundary = new bool[mesh->LocalNx];
  is_upper_boundary = new bool[mesh->LocalNx];
  has_lower_boundary = false;
  has_upper_boundary = false;
  
  for (int jx=0; jx<mesh->LocalNx; jx++) {
    is_lower_boundary[jx] = mesh->firstY(jx);
    if (is_lower_boundary[jx]) has_lower_boundary = true;
    is_upper_boundary[jx] = mesh->lastY(jx);
    if (is_upper_boundary[jx]) has_upper_boundary = true;
  }
  
  nx_core = 0;
  nx_sol = 0;
  nz = mesh->LocalNz;
  for (int jx=mesh->xstart; jx<=mesh->xend; jx++) {
    if (mesh->periodicY(jx)) nx_core++;
    else nx_sol++;
  }
  
  boundary_gradient_smoothing_length = mesh->GlobalNy/256; // Number of grid points to average over for the gradient used to extrapolate to the guard cells
  boundary_condition_smoothing_range = BoutReal(mesh->GlobalNy)/64.; // Scale length of the exponential decay used to smoothly apply the boundary condition
  if (boundary_gradient_smoothing_length==0)
    boundary_gradient_smoothing_length = 1;
  
  if (fluxes_location_is_ylow && !mesh->StaggerGrids) {
    throw BoutException("Trying to calculate the heat flux at CELL_YLOW while StaggerGrids=false is an error.");
  }
  
  // Get the options for the model
  if (options==NULL) {
    options = Options::getRoot()->getSection("non_local_parallel");
  }
  
  OPTION(options, moments_number, 10);
  OPTION(options, NONLOCAL_PARALLEL_TAGBASE, 12381);
  OPTION(options, gradT_drive, false);
  OPTION(options, gradV_drive, false);
  OPTION(options, VeminusVi_drive, false);
  OPTION(options, calculate_heatflux, false);
  OPTION(options, calculate_viscosity, false);
  OPTION(options, calculate_friction, false);
  OPTION(options, yperiodic, false);
  OPTION(options, bc_heatflux, false);
  OPTION(options, bc_viscosity, false);
  OPTION(options, bc_apply_reflecting_sheath, false);
  OPTION(options, maxwellian_source_drives, false);
  
  if (!gradT_drive && !gradV_drive &&!VeminusVi_drive) {
    throw BoutException("NonLocalParallel has no drive terms switched on");
  }
  if (maxwellian_source_drives) {
    number_of_source_drives = moments_number-2;
  }
  number_of_drives = int(gradT_drive)+int(gradV_drive)+int(VeminusVi_drive)+int(maxwellian_source_drives)*number_of_source_drives;
  
  if (!calculate_heatflux && !calculate_viscosity && !calculate_friction) {
    throw BoutException("NonLocalParallel is not calculating anything");
  }
  
  if (yperiodic) {
    if (bc_heatflux || bc_viscosity || bc_apply_reflecting_sheath) {
      throw BoutException("boundary condition switches in NonLocalParallel are inconsistent: cannot set both periodic and toroidal/sheath boundary conditions");
    }
  }

  cubic_spline_inverse_lambdaC.initialise('y',true,fluxes_location_is_ylow);
  cubic_splines_driveterms_centre = new CubicSpline[number_of_drives];
  int cubic_splines_counter = 0;
  if (yperiodic) {
    if (fluxes_location_is_ylow) {
      cubic_splines_driveterms_ylow = new CubicSpline[number_of_drives];
      if (gradT_drive) {
	cubic_splines_driveterms_centre[cubic_splines_counter].initialise_periodic('y'); // gradT_driveterm coefficient; CELL_CENTRE quantity, so must be adjusted when the grids are staggered.
	cubic_splines_driveterms_ylow[cubic_splines_counter].initialise_periodic('y'); // gradT; will fix up the boundary guard cells with forward/backward derivatives, since at least one guard cell value will be needed
	cubic_splines_counter++;
      }
      if (gradV_drive) {
	cubic_splines_driveterms_centre[cubic_splines_counter].initialise_periodic('y'); // gradV_driveterm; CELL_CENTRE quantity, so must be adjusted when the grids are staggered.
	cubic_splines_driveterms_ylow[cubic_splines_counter].initialise_periodic('y'); // dummy coefficient with value one everywhere as gradV is CELL_CENTRE so there is no CELL_YLOW term needed
	Field3D one = 1.;
	cubic_splines_driveterms_ylow[cubic_splines_counter].calculate(one);
	cubic_splines_counter++;
      }
      if (VeminusVi_drive) {
	cubic_splines_driveterms_centre[cubic_splines_counter].initialise_periodic('y'); // VeminusVi_driveterm coefficient
	cubic_splines_driveterms_ylow[cubic_splines_counter].initialise_periodic('y'); // VeminusVi
	cubic_splines_counter++;
      }
      if (maxwellian_source_drives) {
	for (int source_counter=0; source_counter<number_of_source_drives; source_counter++) {
	  cubic_splines_driveterms_centre[cubic_splines_counter].initialise_periodic('y'); // source_driveterms
	  cubic_splines_driveterms_ylow[cubic_splines_counter].initialise_periodic('y'); // dummy coefficient with value one everywhere as source_driveterms is CELL_CENTRE so there is no CELL_YLOW term needed
	  Field3D one = 1.;
	  cubic_splines_driveterms_ylow[cubic_splines_counter].calculate(one);
	  cubic_splines_counter++;
	}
      }
    } else {
      if (gradT_drive) {
	cubic_splines_driveterms_centre[cubic_splines_counter].initialise_periodic('y'); // gradT_driveterm; false because gradT_driveterm includes a derivative, so we don't know it in the guard cells, hence calculate the interpolation excluding the guard cells at the target boundaries.
	cubic_splines_counter++;
      }
      if (gradV_drive) {
	cubic_splines_driveterms_centre[cubic_splines_counter].initialise_periodic('y'); // gradV_driveterm; false because gradT_driveterm includes a derivative, so we don't know it in the guard cells, hence calculate the interpolation excluding the guard cells at the target boundaries.
	cubic_splines_counter++;
      }
      if (VeminusVi_drive) {
	cubic_splines_driveterms_centre[cubic_splines_counter].initialise_periodic('y'); // VeminusVi_driveterm
	cubic_splines_counter++;
      }
      if (maxwellian_source_drives) {
	for (int source_counter=0; source_counter<number_of_source_drives; source_counter++) {
	  cubic_splines_driveterms_centre[cubic_splines_counter].initialise_periodic('y'); // source_driveterms
	  cubic_splines_counter++;
	}
      }
    }
  } else {
    // NB if cubic_spline_*.periodic_boundary==false, the cubic_spline algorithm will use an algorithm identical to the periodic one on closed flux surfaces, as all the other stuff only gets done where there is a BndryLowerY or BndryUpperY, so we do not need to deal with them separately
    if (fluxes_location_is_ylow) {
      cubic_splines_driveterms_ylow = new CubicSpline[number_of_drives];
      if (gradT_drive) {
	cubic_splines_driveterms_centre[cubic_splines_counter].initialise('y',true,true); // gradT_driveterm coefficient; CELL_CENTRE quantity, so must be adjusted when the grids are staggered.
	cubic_splines_driveterms_ylow[cubic_splines_counter].initialise('y',true,false); // gradT; will fix up the boundary guard cells with forward/backward derivatives, since at least one guard cell value will be needed
	cubic_splines_counter++;
      }
      if (gradV_drive) {
	cubic_splines_driveterms_centre[cubic_splines_counter].initialise('y',false,true); // gradV_driveterm; CELL_CENTRE quantity, so must be adjusted when the grids are staggered.
	cubic_splines_driveterms_ylow[cubic_splines_counter].initialise('y',true,false); // dummy coefficient with value one everywhere as gradV is CELL_CENTRE so there is no CELL_YLOW term needed
	Field3D one = 1.;
	cubic_splines_driveterms_ylow[cubic_splines_counter].calculate(one);
	cubic_splines_counter++;
      }
      if (VeminusVi_drive) {
	cubic_splines_driveterms_centre[cubic_splines_counter].initialise('y',true,true); // VeminusVi_driveterm coefficient
	cubic_splines_driveterms_ylow[cubic_splines_counter].initialise('y',true,false); // VeminusVi
	cubic_splines_counter++;
      }
      if (maxwellian_source_drives) {
	for (int source_counter=0; source_counter<number_of_source_drives; source_counter++) {
	  cubic_splines_driveterms_centre[cubic_splines_counter].initialise('y',true,true); // source_driveterms
	  cubic_splines_driveterms_ylow[cubic_splines_counter].initialise('y',true,false); // dummy coefficient with value one everywhere as source_driveterms is CELL_CENTRE so there is no CELL_YLOW term needed
	  Field3D one = 1.;
	  cubic_splines_driveterms_ylow[cubic_splines_counter].calculate(one);
	  cubic_splines_counter++;
	}
      }
    } else {
      if (gradT_drive) {
	cubic_splines_driveterms_centre[cubic_splines_counter].initialise('y',false,false); // gradT_driveterm; false because gradT_driveterm includes a derivative, so we don't know it in the guard cells, hence calculate the interpolation excluding the guard cells at the target boundaries.
	cubic_splines_counter++;
      }
      if (gradV_drive) {
	cubic_splines_driveterms_centre[cubic_splines_counter].initialise('y',false,false); // gradV_driveterm; false because gradT_driveterm includes a derivative, so we don't know it in the guard cells, hence calculate the interpolation excluding the guard cells at the target boundaries.
	cubic_splines_counter++;
      }
      if (VeminusVi_drive) {
	cubic_splines_driveterms_centre[cubic_splines_counter].initialise('y',true,false); // VeminusVi_driveterm
	cubic_splines_counter++;
      }
      if (maxwellian_source_drives) {
	for (int source_counter=0; source_counter<number_of_source_drives; source_counter++) {
	  cubic_splines_driveterms_centre[cubic_splines_counter].initialise('y',true,false); // source_driveterms
	  cubic_splines_counter++;
	}
      }
    }
  }
  
  position = new bindex;
  broadcast_request = MPI_REQUEST_NULL;

  elementary_charge = pass_elementary_charge;
  electron_mass = pass_electron_mass;
//   ion_mass = pass_ion_mass;
  epsilon_0 = pass_epsilon_0;
  logLambda = pass_logLambda;
//   NonLocalParallelgamma_factor = pass_gamma_factor;
  
  if (fluxes_location_is_ylow) {
    increasing_dimensionless_length.setLocation(CELL_YLOW);
    decreasing_dimensionless_length.setLocation(CELL_YLOW);
    if (calculate_heatflux) {
      electron_heat_flux.setLocation(CELL_YLOW);
    }
    if (calculate_viscosity) {
      electron_viscosity.setLocation(CELL_YLOW); // should really have this as CELL_CENTRE quantity, but then would have to calculate the integrals on both CELL_CENTRE and CELL_YLOW which is perfectly possible but has not been implemented yet
    }
    if (calculate_friction) {
      electron_friction.setLocation(CELL_YLOW);
    }
    if (gradT_drive) {
      gradT_electron.setLocation(CELL_YLOW);
      gradT_electron=0.;
    }
    // No CELL_YLOW drive terms for gradV_drive
    // VeminusVi already CELL_YLOW
  }
  
  if (calculate_heatflux) {
    electron_heat_flux = 0.;
  }
  if (calculate_viscosity) {
    electron_viscosity = 0.;
  }
  if (calculate_friction) {
    electron_friction = 0.;
  }
  if (gradT_drive) {
    gradT_driveterm = 0.;
  }
  if (gradV_drive) {
    gradV_driveterm = 0.;
  }
  if (VeminusVi_drive) {
    VeminusVi_driveterm = 0.;
  }
  lambdaC_inverse = 0.;
  increasing_dimensionless_length = 0.;
  decreasing_dimensionless_length = 0.;
  dimensionless_length_deltas_above = 0.;
  if (fluxes_location_is_ylow) {
    dimensionless_length_deltas_below=0.;
  }
  total_dimensionless_length = 1./0.;
  
  integration.initialise(fluxes_location_is_ylow);
  
  // Get model's eigenvalues and coefficients from files
  BoutReal junk;
  bool allocated_eigenvalues = false;
  if (calculate_heatflux) {
    std::stringstream heatflux_infilename;
    heatflux_infilename<<"nonlocal_coefficients/heatfluxcoeffs"<<moments_number;
    std::ifstream heatflux_infile ( heatflux_infilename.str().c_str() );
    if (!heatflux_infile.is_open()) {
      throw BoutException("Could not open heatfluxcoeffs file");
    }
    heatflux_infile>>number_of_negative_eigenvalues;
    if (!allocated_eigenvalues) {
      eigenvalues = new BoutReal[number_of_negative_eigenvalues];
      allocated_eigenvalues = true;
    }
    heatflux_coefficients_below = new BoutReal[number_of_negative_eigenvalues];
    heatflux_coefficients_above = new BoutReal[number_of_negative_eigenvalues];
    heatflux_infile>>heatflux_zerocoeff;
    for (int i=0; i<number_of_negative_eigenvalues; i++) {
      if (heatflux_infile.eof()) {
	throw BoutException("reached end of heatfluxcoeffs file unexpectedly");
      }
      heatflux_infile>>eigenvalues[i];
      heatflux_infile>>heatflux_coefficients_below[i];
      heatflux_coefficients_above[i] = -heatflux_coefficients_below[i];
    }
    heatflux_infile.close();
    std::stringstream heatfluxbc_infilename;
    heatfluxbc_infilename<<"nonlocal_coefficients/heatfluxbc"<<moments_number;
    std::ifstream heatfluxbc_infile ( heatfluxbc_infilename.str().c_str() );
    if (!heatfluxbc_infile.is_open()) {
      throw BoutException("Could not open heatfluxbc file");
    }
    
    if (bc_heatflux && (nx_sol>0)) {
      WinverseB_11 = new BoutReal[number_of_negative_eigenvalues];
    }
    
    if (bc_heatflux && (nx_sol>0)) {
      heatfluxbc_infile>>W11_dot_W11;
    } else {
      heatfluxbc_infile>>junk;
    }
    
    if (bc_heatflux && (nx_sol>0)) {
      heatfluxbc_infile>>W11_dot_W20;
    } else {
      heatfluxbc_infile>>junk;
    }
    
    for (int i=0; i<number_of_negative_eigenvalues; i++) {
      if (heatfluxbc_infile.eof()) {
	output<<"Error at i="<<i<<endl;
	throw BoutException("reached end of heatfluxbc file unexpectedly");
      }
      
      if (bc_heatflux && (nx_sol>0)) {
	heatfluxbc_infile>>WinverseB_11[i];
      } else {
	heatfluxbc_infile>>junk;
      }
    }
    heatfluxbc_infile.close();
  }
  if (bc_heatflux && (nx_sol>0)) {
    if (has_lower_boundary || has_upper_boundary) {
      pass_interim_upper_boundary_n11 = 0.;
      upper_boundary_condition_n11 = 0.;
    }
  }
  if (calculate_viscosity) {
    std::stringstream viscosity_infilename;
    viscosity_infilename<<"nonlocal_coefficients/viscositycoeffs"<<moments_number;
    std::ifstream viscosity_infile ( viscosity_infilename.str().c_str() );
    if (!viscosity_infile.is_open()) {
      throw BoutException("Could not open viscositycoeffs file");
    }
    viscosity_infile>>number_of_negative_eigenvalues;
    if (!allocated_eigenvalues) {
      eigenvalues = new BoutReal[number_of_negative_eigenvalues];
      allocated_eigenvalues = true;
    }
    viscosity_coefficients = new BoutReal[number_of_negative_eigenvalues];
    viscosity_infile>>junk; // the zero eigenvector has only odd-l components, so does not contribute to the viscosity (2,0)
    for (int i=0; i<number_of_negative_eigenvalues; i++) {
      if (viscosity_infile.eof()) {
	throw BoutException("reached end of viscositycoeffs file unexpectedly");
      }
      viscosity_infile>>eigenvalues[i];
      viscosity_infile>>viscosity_coefficients[i];
    }
    viscosity_infile.close();
    std::stringstream viscositybc_infilename;
    viscositybc_infilename<<"nonlocal_coefficients/viscositybc"<<moments_number;
    std::ifstream viscositybc_infile ( viscositybc_infilename.str().c_str() );
    if (!viscositybc_infile.is_open()) {
      throw BoutException("Could not open viscositybc file");
    }
    if (bc_heatflux && (nx_sol>0)) {
      WinverseB_20 = new BoutReal[number_of_negative_eigenvalues];
    }
    
    if (bc_heatflux && (nx_sol>0)) {
      viscositybc_infile>>W20_dot_W20;
    } else {
      viscositybc_infile>>junk;
    }
    
    if (bc_heatflux && (nx_sol>0)) {
      viscositybc_infile>>W20_dot_W11;
    } else {
      viscositybc_infile>>junk;
    }
    
    for (int i=0; i<number_of_negative_eigenvalues; i++) {
      if (viscositybc_infile.eof()) {
	throw BoutException("reached end of viscositybc file unexpectedly");
      }
      
      if (bc_heatflux && (nx_sol>0)) {
	viscositybc_infile>>WinverseB_20[i];
      } else {
	viscositybc_infile>>junk;
      }
    }
    viscositybc_infile.close();
  }
  if (bc_viscosity && (nx_sol>0)) {
    if (has_lower_boundary || has_upper_boundary) {
      pass_interim_upper_boundary_n20 = 0.;
      upper_boundary_condition_n20 = 0.;
    }
  }
  if (calculate_friction) {
    std::stringstream friction_infilename;
    friction_infilename<<"nonlocal_coefficients/frictioncoeffs"<<moments_number;
    std::ifstream friction_infile ( friction_infilename.str().c_str() );
    if (!friction_infile.is_open()) {
      throw BoutException("Could not open frictioncoeffs file");
    }
    friction_infile>>number_of_negative_eigenvalues;
    if (!allocated_eigenvalues) {
      eigenvalues = new BoutReal[number_of_negative_eigenvalues];
      allocated_eigenvalues = true;
    }
    friction_coefficients_below = new BoutReal[number_of_negative_eigenvalues];
    friction_coefficients_above = new BoutReal[number_of_negative_eigenvalues];
    friction_infile>>friction_zerocoeff;
    for (int i=0; i<number_of_negative_eigenvalues; i++) {
      if (friction_infile.eof()) {
	throw BoutException("reached end of frictioncoeffs file unexpectedly");
      }
      friction_infile>>eigenvalues[i];
      friction_infile>>friction_coefficients_below[i];
      friction_coefficients_above[i] = -friction_coefficients_below[i];
    }
    friction_infile.close();
  }
  
  driveterm_coefficients_below = new BoutReal[number_of_drives*number_of_negative_eigenvalues]; // There must be at least one output, so there will be a valid number_of_negative_eigenvalues by this point
  driveterm_coefficients_above = new BoutReal[number_of_drives*number_of_negative_eigenvalues];
  
  int driveterm_counter = 0;
  if (gradT_drive) {
    std::stringstream gradT_drive_infilename;
    gradT_drive_infilename<<"nonlocal_coefficients/gradTdrivecoeffs"<<moments_number;
    std::ifstream gradT_drive_infile ( gradT_drive_infilename.str().c_str() );
    if (!gradT_drive_infile.is_open()) {
      throw BoutException("Could not open gradTdrivecoeffs file");
    }
    gradT_drive_infile>>number_of_negative_eigenvalues;
    if (!allocated_eigenvalues) {
      eigenvalues = new BoutReal[number_of_negative_eigenvalues];
      allocated_eigenvalues = true;
    }
    gradT_drive_infile>>gradT_zerocoeff;
    for (int i=0; i<number_of_negative_eigenvalues; i++) {
      if (gradT_drive_infile.eof()) {
	throw BoutException("reached end of gradTdrivecoeffs file unexpectedly");
      }
      gradT_drive_infile>>eigenvalues[i];
      gradT_drive_infile>>driveterm_coefficients_below[i*number_of_drives+driveterm_counter];
      driveterm_coefficients_above[i*number_of_drives+driveterm_counter] = -driveterm_coefficients_below[i*number_of_drives+driveterm_counter]; // -ve sign as l=1 is odd for gradT drive
    }
    gradT_drive_infile.close();
    driveterm_counter++;
  }
  if (gradV_drive) {
    std::stringstream gradV_drive_infilename;
    gradV_drive_infilename<<"nonlocal_coefficients/gradVdrivecoeffs"<<moments_number;
    std::ifstream gradV_drive_infile ( gradV_drive_infilename.str().c_str() );
    if (!gradV_drive_infile.is_open()) {
      throw BoutException("Could not open gradVdrivecoeffs file");
    }
    gradV_drive_infile>>number_of_negative_eigenvalues;
    if (!allocated_eigenvalues) {
      eigenvalues = new BoutReal[number_of_negative_eigenvalues];
      allocated_eigenvalues = true;
    }
    gradV_drive_infile>>junk; // gradV drive has l=2, hence zerocoeff is zero
    for (int i=0; i<number_of_negative_eigenvalues; i++) {
      if (gradV_drive_infile.eof())
	throw BoutException("reached end of gradVdrivecoeffs file unexpectedly");
      gradV_drive_infile>>eigenvalues[i];
      gradV_drive_infile>>driveterm_coefficients_below[i*number_of_drives+driveterm_counter];
      driveterm_coefficients_above[i*number_of_drives+driveterm_counter] = driveterm_coefficients_below[i*number_of_drives+driveterm_counter]; // +ve sign as l=2 is even for gradV drive
    }
    gradV_drive_infile.close();
    driveterm_counter++;
  }
  if (VeminusVi_drive) {
    std::stringstream VeminusVi_drive_infilename;
    VeminusVi_drive_infilename<<"nonlocal_coefficients/VeminusVidrivecoeffs"<<moments_number;
    std::ifstream VeminusVi_drive_infile ( VeminusVi_drive_infilename.str().c_str() );
    if (!VeminusVi_drive_infile.is_open())
      throw BoutException("Could not open VeminusVidrivecoeffs file");
    VeminusVi_drive_infile>>number_of_negative_eigenvalues;
    if (!allocated_eigenvalues) {
      eigenvalues = new BoutReal[number_of_negative_eigenvalues];
      allocated_eigenvalues = true;
    }
    VeminusVi_drive_infile>>VeminusVi_zerocoeff;
    for (int i=0; i<number_of_negative_eigenvalues; i++) {
      if (VeminusVi_drive_infile.eof())
	throw BoutException("reached end of VeminusVidrivecoeffs file unexpectedly");
      VeminusVi_drive_infile>>eigenvalues[i];
      VeminusVi_drive_infile>>driveterm_coefficients_below[i*number_of_drives+driveterm_counter];
      driveterm_coefficients_above[i*number_of_drives+driveterm_counter] = -driveterm_coefficients_below[i*number_of_drives+driveterm_counter]; // -ve sign as l=1 is odd for VeminusVi drive
    }
    VeminusVi_drive_infile.close();
    driveterm_counter++;
  }
  if (maxwellian_source_drives) {
    source_driveterms = new Field3D[number_of_source_drives];
    std::stringstream source_drives_infilename;
    source_drives_infilename<<"nonlocal_coefficients/maxwelliansourcedrivescoeffs"<<moments_number;
    std::ifstream source_drives_infile ( source_drives_infilename.str().c_str() );
    if (!source_drives_infile.is_open())
      throw BoutException("Could not open maxwelliansourcedrivescoeffs file");
    for (int source_counter=0; source_counter<number_of_source_drives; source_counter++) {
      // maxwellian source_drives have l=0 which is even, hence zerocoeffs are zero
      for (int i=0; i<number_of_negative_eigenvalues; i++) {
	if (source_drives_infile.eof())
	  throw BoutException("reached end of maxwelliansourcedrivescoeffs file unexpectedly");
	source_drives_infile>>driveterm_coefficients_below[i*number_of_drives+driveterm_counter];
	driveterm_coefficients_above[i*number_of_drives+driveterm_counter] = driveterm_coefficients_below[i*number_of_drives+driveterm_counter]; // +ve sign as l=0 is even for source_drives
      }
      driveterm_counter++;
    }
    source_drives_infile.close();
  }
  
  if (has_lower_boundary) {
    exp_total_dimensionless_length_over_eigenvalue = new BoutReal[number_of_negative_eigenvalues];
  }
  
  if (!yperiodic && (nx_core>0)) {
    toroidal_solver_forward = new NonLocalParallelToroidalSolver(false);
    toroidal_solver_reverse = new NonLocalParallelToroidalSolver(true);
  }
  transients_row_size = number_of_negative_eigenvalues*nz;
  transients_data_size = (mesh->xend-mesh->xstart+1)*transients_row_size;
  lower_transients_data = new BoutReal[transients_data_size];
  upper_transients_data = new BoutReal[transients_data_size];
  if (bc_heatflux && (nx_sol>0)) {
    heatflux_sol_transients_factors = new BoutReal[2*nx_sol*(nz)];
  }
  if (bc_viscosity && (nx_sol>0)) {
    viscosity_sol_transients_factors = new BoutReal[2*nx_sol*(nz)];
  }
  
  ycomms = new MPI_Comm[mesh->xend-mesh->xstart+1];
  ycomms_rank = new int[mesh->xend-mesh->xstart+1];
  ycomms_size = new int[mesh->xend-mesh->xstart+1];
  for (int jx=mesh->xstart,i=0;jx<=mesh->xend;jx++,i++) {
    if (i>=(mesh->xend-mesh->xstart+1)) throw BoutException("Number of points in x exceeded mxsub");
    ycomms[i] = mesh->getYcomm(jx);
    MPI_Comm_rank(ycomms[i],ycomms_rank+i);
    MPI_Comm_size(ycomms[i],ycomms_size+i);
  }
//   // Initialisation for stuff used in y_broadcast functions
//   {
//     MPI_Group group_yprocs;
//     
//     int n_yprocs = mesh->getNYPE();
//     int * indices_yprocs = new int[n_yprocs];
//     for (int i=0; i<n_yprocs; i++)
//       indices_yprocs[i] = i * mesh->getNXPE() + mesh->getXProcIndex();
//     
//     MPI_Group group_world;
//     MPI_Comm_group(BoutComm::get(), &group_world); // Get the entire group
//     MPI_Group_incl(group_world, n_yprocs, indices_yprocs, &group_yprocs);
//     MPI_Group_free(&group_world);
//     delete [] indices_yprocs;
//     
//     MPI_Comm_create(BoutComm::get(), group_yprocs, &comm_yprocs);
//     MPI_Group_free(&group_yprocs);
//   }
//   {
//     MPI_Group group_yprocs;
//     
//     int n_yprocs = mesh->getNYPE()-1;
//     int * indices_yprocs = new int[n_yprocs];
//     for (int i=0; i<n_yprocs; i++)
//       indices_yprocs[i] = (i+1) * mesh->getNXPE() + mesh->getXProcIndex();
//     
//     MPI_Group group_world;
//     MPI_Comm_group(BoutComm::get(), &group_world); // Get the entire group
//     MPI_Group_incl(group_world, n_yprocs, indices_yprocs, &group_yprocs);
//     MPI_Group_free(&group_world);
//     delete [] indices_yprocs;
//     
//     MPI_Comm_create(BoutComm::get(), group_yprocs, &comm_yprocs_minusone);
//     MPI_Group_free(&group_yprocs);
//   }
  
  if (calculate_heatflux) {
    output<<"	NonLocalParallel: calculating heat-flux"<<endl;
  }
  if (calculate_viscosity) {
    output<<"	NonLocalParallel: calculating viscosity"<<endl;
  }
  if (calculate_friction) {
    output<<"	NonLocalParallel: calculating friction"<<endl;
  }
  if (!yperiodic && (nx_core>0)) {
    output<<"	NonLocalParallel: applying toroidal boundary conditions in core"<<endl;
  }
  if (yperiodic) {
    output<<"	NonLocalParallel: applying periodic boundary conditions"<<endl;
  }
  if (bc_heatflux && (nx_sol>0)) {
    output<<"	NonLocalParallel: applying boundary condition to heat-flux in SOL"<<endl;
  }
  if (bc_viscosity && (nx_sol>0)) {
    output<<"	NonLocalParallel: applying boundary condition to viscosity in SOL"<<endl;
  }
  if (bc_apply_reflecting_sheath && (nx_sol>0)) {
    output<<"	NonLocalParallel: applying reflecting sheath boundary conditions in SOL"<<endl;
  }
  if (gradT_drive) {
    output<<"	NonLocalParallel: including gradT_driveterm"<<endl;
  }
  if (gradV_drive) {
    output<<"	NonLocalParallel: including gradV_driveterm"<<endl;
  }
  if (VeminusVi_drive) {
    output<<"	NonLocalParallel: including VeminusVi_driveterm"<<endl;
  }
  if (maxwellian_source_drives) {
    output<<"	NonLocalParallel: including maxwellian source_drives"<<endl;
  }
}


NonLocalParallel::~NonLocalParallel() {
  delete [] eigenvalues;
  delete [] is_lower_boundary;
  delete [] is_upper_boundary;
  delete [] ycomms_rank;
  delete [] ycomms_size;
  delete [] cubic_splines_driveterms_centre;
  if (calculate_heatflux) {
    delete [] heatflux_coefficients_below;
    delete [] heatflux_coefficients_above;
  }
  if (fluxes_location_is_ylow) delete [] cubic_splines_driveterms_ylow;
  if (bc_heatflux && (nx_sol>0)) {
    delete [] heatflux_sol_transients_factors;
    delete [] WinverseB_11;
  }
  if (calculate_viscosity) {
    delete [] viscosity_coefficients;
  }
  if (bc_viscosity && (nx_sol>0)) {
    delete [] viscosity_sol_transients_factors;
    delete [] WinverseB_20;
  }
  if (calculate_friction) {
    delete [] friction_coefficients_below;
    delete [] friction_coefficients_above;
  }
  delete [] driveterm_coefficients_below;
  delete [] driveterm_coefficients_above;
  if (!yperiodic && (nx_core>0)) {
    delete toroidal_solver_forward;
    delete toroidal_solver_reverse;
  }
  delete [] lower_transients_data;
  delete [] upper_transients_data;
  if (has_lower_boundary) {
    delete [] exp_total_dimensionless_length_over_eigenvalue;
  }
//   MPI_Comm_free(&comm_yprocs);
//   MPI_Comm_free(&comm_yprocs_minusone);
}


/**********************************************************************************
 *                    HEAT FLUX CALCULATION ROUTINES
 **********************************************************************************/

void NonLocalParallel::calculate_nonlocal_closures() {
  if (fluxes_location_is_ylow)
    calculate_nonlocal_closures_cell_ylow();
  else
    calculate_nonlocal_closures_cell_centre();
  #ifdef CHECK
    calculated_before_setting_bcs=true;
  #endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NonLocalParallel::calculate_nonlocal_closures_cell_centre() {
  int dataindex=0; // Use for cycling through jx indices of lower/upper_transients_data
  
  lambdaC_inverse = n_electron * pow(elementary_charge,4) * logLambda / 12. / pow(PI,1.5) / pow(epsilon_0,2) / SQ(T_electron);
  if (normalized_density) lambdaC_inverse *= density_normalization; // If set_density_normalization has been called, un-normalize the density used to calculate the collision length so that it has the correct units
  
  if (gradT_drive) {
    gradT_driveterm = 5./4. * n_electron / T_electron * Grad_par(T_electron) / lambdaC_inverse; //g^(1,1)
    mesh->communicate(gradT_driveterm);
  }
  if (gradV_drive) {
    gradV_driveterm = -0.5 * n_electron / sqrt(2.*T_electron/electron_mass) * Grad_par(V_electron) / lambdaC_inverse;
    mesh->communicate(gradV_driveterm);
  }
  if (VeminusVi_drive) {
    VeminusVi_driveterm = -2./sqrt(PI) * (-j_parallel/elementary_charge) / sqrt(2.*T_electron/electron_mass);
//     mesh->communicate(VeminusVi_driveterm);
  }
  
  // Now calculate z and deltaz everywhere
  cubic_spline_inverse_lambdaC.calculate(lambdaC_inverse);
  
  start_index(position);
  
  do {
    
    if (!is_lower_boundary[position->jx]) {
      {
	Timer timer("comms");
	if (position->jx < mesh->DownXSplitIndex())
	  mesh->wait(mesh->irecvYInIndest(&increasing_dimensionless_length(position->jx,mesh->ystart,position->jz),1,
					  NONLOCAL_PARALLEL_TAGBASE + position->jx*mesh->LocalNz+position->jz));
	else
	  mesh->wait(mesh->irecvYInOutdest(&increasing_dimensionless_length(position->jx,mesh->ystart,position->jz),1,
					  NONLOCAL_PARALLEL_TAGBASE + position->jx*mesh->LocalNz+position->jz));
      }
    }
    
    position->jy = mesh->ystart-1;
    calc_index(position);

    interp_coefficients = cubic_spline_inverse_lambdaC.coefficients(position);
    // d/dy(delta) = 1/lambdaC = a + b*t + c*t^2 + d*t^3; t=(ind-jy)=(y-y0)/(sqrt(g_22)*dy); ind is a notional continuous variable equal to jy at the gridpoints so at jy+1 t=1
    dimensionless_length_deltas_above[*position] /* = dy/dt*(a + 1/2*b + 1/3*c + 1/4*d) */
      = coord->dy(position->jx,position->jy)*sqrt(0.5*(coord->g_22(position->jx,position->jy)+coord->g_22(position->jx,position->jyp)))
      *(interp_coefficients[0] + interp_coefficients[1]/2. + interp_coefficients[2]/3. + interp_coefficients[3]/4.);
    next_index_y(position);
    
    do {
      interp_coefficients = cubic_spline_inverse_lambdaC.coefficients(position);
      // d/dy(delta) = 1/lambdaC = a + b*t + c*t^2 + d*t^3; t=(ind-jy)=(y-y0)/(sqrt(g_22)*dy); ind is a notional continuous variable equal to jy at the gridpoints so at jy+1 t=1
      dimensionless_length_deltas_above[*position] /* = dy/dt*(a + 1/2*b + 1/3*c + 1/4*d) */
        = coord->dy(position->jx,position->jy)*sqrt(0.5*(coord->g_22(position->jx,position->jy)+coord->g_22(position->jx,position->jyp)))
        *(interp_coefficients[0] + interp_coefficients[1]/2. + interp_coefficients[2]/3. + interp_coefficients[3]/4.);
      increasing_dimensionless_length(position->jx,position->jyp,position->jz) = increasing_dimensionless_length[*position] + dimensionless_length_deltas_above[*position];
    } while (next_index_y(position));
    
    {
      Timer timer("comms");
      if (position->jx<mesh->UpXSplitIndex())
	mesh->sendYOutIndest(&increasing_dimensionless_length(position->jx,mesh->yend+1,position->jz),1,
			      NONLOCAL_PARALLEL_TAGBASE + position->jx*mesh->LocalNz+position->jz);
      else
	mesh->sendYOutOutdest(&increasing_dimensionless_length(position->jx,mesh->yend+1,position->jz),1,
			      NONLOCAL_PARALLEL_TAGBASE + position->jx*mesh->LocalNz+position->jz);
    }
    
  } while (next_indexperp(position));
  
  
  // Send the total dimensionless_length at the upper boundary back to the other processors.
  if (has_upper_boundary) {
    total_dimensionless_length = sliceXZ(increasing_dimensionless_length,mesh->yend);
  }
  
  for (int jx=mesh->xstart; jx<=mesh->xend; jx++) {
    y_broadcast(&(total_dimensionless_length(jx,0)), // get a pointer to the jx'th row of total_dimensionless_length's data
		nz, ycomms_size[jx-mesh->xstart]-1,jx);
  }
  
  decreasing_dimensionless_length = -increasing_dimensionless_length;
  
  for (auto &i : decreasing_dimensionless_length.region(RGN_NOY)) {
    decreasing_dimensionless_length[i] += total_dimensionless_length[i];
  }

  if (has_lower_boundary) {
    total_dimensionless_length.setIndex(mesh->ystart-1);
    FieldPerp tmp = sliceXZ(dimensionless_length_deltas_above, mesh->ystart-1);
    for(auto &i : total_dimensionless_length) {
      decreasing_dimensionless_length[i] += total_dimensionless_length[i];
      decreasing_dimensionless_length[i] += tmp[i];
    }
  }
  
  if (calculate_heatflux) {
    electron_heat_flux = 0.;
  }
  if (calculate_viscosity) {
    electron_viscosity = 0.;
  }
  if (calculate_friction) {
    electron_friction = 0.;
  }
// BoutReal temp_ef = 0.;
// output<<"first "<<electron_friction[40][35][0]<<endl;
// temp_ef=electron_friction[40][35][0];
  int driveterm_counter = 0;
  if (gradT_drive) {
    if (calculate_heatflux) {
      electron_heat_flux += -heatflux_zerocoeff * gradT_zerocoeff * gradT_driveterm; //zero eigenvalue contribution to n^(1,1) //zero eigenvalue contribution to n^(1,1)/T^1.
    }
    // viscosity gets no zero eigenvalue contribution
    if (calculate_friction) {
      electron_friction += -friction_zerocoeff * gradT_zerocoeff * gradT_driveterm;
    }
    cubic_splines_driveterms_centre[driveterm_counter].calculate(gradT_driveterm);
    driveterm_counter++;
  }
  if (gradV_drive) {
    // gradV drive has no zero eigenvalue contribution
    cubic_splines_driveterms_centre[driveterm_counter].calculate(gradV_driveterm);
    driveterm_counter++;
  }
  if (VeminusVi_drive) {
    if (calculate_heatflux) {
      electron_heat_flux += -heatflux_zerocoeff * VeminusVi_zerocoeff * VeminusVi_driveterm;
    }
    // viscosity gets no zero eigenvalue contribution
    if (calculate_friction) {
      electron_friction += -friction_zerocoeff * VeminusVi_zerocoeff * VeminusVi_driveterm;
    }
    cubic_splines_driveterms_centre[driveterm_counter].calculate(VeminusVi_driveterm);
    driveterm_counter++;
  }
  if (maxwellian_source_drives) {
    calculate_source_driveterms();
    // source_driveterms have no zero eigenvalue contribution
    for (int source_counter=0; source_counter<number_of_source_drives; source_counter++) {
      cubic_splines_driveterms_centre[driveterm_counter].calculate(source_driveterms[source_counter]);
      driveterm_counter++;
    }
  }
  
// BoutReal sumv = 0.;
// BoutReal sumf = 0.;
// BoutReal sumv2 = 0.;
// BoutReal sumf2 = 0.;
// for (int i=0;i<number_of_negative_eigenvalues;i++){
//   sumv += viscosity_gradV_coefficients[i]*eigenvalues[i];
//   sumf += friction_gradV_coefficients[i]*eigenvalues[i];
//   sumv2 += viscosity_gradV_coefficients[i];
//   sumf2 += friction_gradV_coefficients[i];
// }
// output<<sumv<<" "<<sumf<<endl<<sumv2<<" "<<sumf2<<endl;
// exit(17);
// output<<"with zerocoeff "<<electron_friction[40][35][0]<<" "<<electron_friction[40][35][0]-temp_ef<<endl;
// temp_ef=electron_friction[40][35][0];
  {
    BoutReal* temppointer = lower_transients_data;
    for (int i=0; i<transients_data_size; i++) {
      *temppointer=0.;
      temppointer++;
    }
  }
  for (int j=0; j<number_of_negative_eigenvalues; j++) {
    integration.calculateIntegralBelow_cell_centre(eigenvalues[j], dimensionless_length_deltas_above, cubic_spline_inverse_lambdaC, number_of_drives, &driveterm_coefficients_below[j*number_of_drives], cubic_splines_driveterms_centre, j);
    if (calculate_heatflux) {
      electron_heat_flux += heatflux_coefficients_below[j]*integration.integral_below;
    }
    if (calculate_viscosity) {
      electron_viscosity += viscosity_coefficients[j]*integration.integral_below;
    }
    if (calculate_friction) {
      electron_friction += friction_coefficients_below[j]*integration.integral_below;
    }
    dataindex = 0;
    if ( yperiodic || (nx_core>0) || (bc_apply_reflecting_sheath && nx_sol>0) )
      for (int jx=mesh->xstart; jx<=mesh->xend; jx++) {
	if (is_upper_boundary[jx]) {
	  if (mesh->periodicY(jx)) {
	    for (int jz=0; jz<nz; jz++) {
	      lower_transients_data[dataindex*transients_row_size+j*nz+jz] += integration.integral_below(jx,mesh->yend+1,jz); // First store total integrals (at first guard cell) in the transients variables
	      dataindex++;
	    }
	  }
	  else if (bc_apply_reflecting_sheath) {
	    for (int jz=0; jz<nz; jz++) {
	      lower_transients_data[dataindex*transients_row_size+j*nz+jz] += integration.integral_below(jx,mesh->yend,jz); // First store total integrals (at last grid cell) in the transients variables
	      dataindex++;
	    }
	  }
	}
      }
  }
  
// output<<"with IntegralBelow "<<electron_friction[40][35][0]<<" "<<electron_friction[40][35][0]-temp_ef<<endl;
// temp_ef=electron_friction[40][35][0];
  // Signs in the following block: the integral is done 'backwards', so start with -'ve sign. [The following signs are now absorbed into the coefficients] Then multiply by (-1)^(l(A)+l(D)), i.e. if drive and output moment have the same l(mod 2) multiply by +1, if different l(mod 2) by -1
  {
    BoutReal* temppointer = upper_transients_data;
    for (int i=0; i<transients_data_size; i++) {
      *temppointer=0.;
      temppointer++;
    }
  }
  for (int j=0; j<number_of_negative_eigenvalues; j++) {
    integration.calculateIntegralAbove_cell_centre(eigenvalues[j], dimensionless_length_deltas_above, cubic_spline_inverse_lambdaC, number_of_drives, &driveterm_coefficients_above[j*number_of_drives], cubic_splines_driveterms_centre, j);
    
    if (calculate_heatflux) {
      electron_heat_flux += -heatflux_coefficients_above[j]*integration.integral_above;
    }
    if (calculate_viscosity) {
      electron_viscosity += -viscosity_coefficients[j]*integration.integral_above;
    }
    if (calculate_friction) {
      electron_friction += -friction_coefficients_above[j]*integration.integral_above;
    }
    
    dataindex = 0;
    if ( yperiodic || (nx_core>0) || (bc_apply_reflecting_sheath && nx_sol>0) )
      for (int jx=mesh->xstart; jx<=mesh->xend; jx++) {
	if ( is_lower_boundary[jx]) {
	  if (mesh->periodicY(jx)) {
	    for (int jz=0; jz<nz; jz++) {
	      if (calculate_heatflux) {
		upper_transients_data[dataindex*transients_row_size+j*nz+jz] += -heatflux_coefficients_above[j]*integration.integral_above(jx,mesh->ystart-1,jz); // First store total integrals in the transients variables
	      }
	      if (calculate_viscosity) {
		upper_transients_data[dataindex*transients_row_size+j*nz+jz] += -viscosity_coefficients[j]*integration.integral_above(jx,mesh->ystart-1,jz); // First store total integrals in the transients variables
	      }
	      if (calculate_friction) {
		upper_transients_data[dataindex*transients_row_size+j*nz+jz] += -friction_coefficients_above[j]*integration.integral_above(jx,mesh->ystart-1,jz); // First store total integrals in the transients variables
	      }
	      dataindex++;
	    }
	  }
	  else if (bc_apply_reflecting_sheath) {
	    for (int jz=0; jz<nz; jz++) {
	      if (calculate_heatflux) {
		upper_transients_data[dataindex*transients_row_size+j*nz+jz] += -heatflux_coefficients_above[j]*integration.integral_above(jx,mesh->ystart,jz); // First store total integrals in the transients variables
	      }
	      if (calculate_viscosity) {
		upper_transients_data[dataindex*transients_row_size+j*nz+jz] += -viscosity_coefficients[j]*integration.integral_above(jx,mesh->ystart,jz); // First store total integrals in the transients variables
	      }
	      if (calculate_friction) {
		upper_transients_data[dataindex*transients_row_size+j*nz+jz] += -friction_coefficients_above[j]*integration.integral_above(jx,mesh->ystart,jz); // First store total integrals in the transients variables
	      }
	      dataindex++;
	    }
	  }
	}
      }
  }
  
  if (bc_apply_reflecting_sheath && nx_sol>0 && (has_lower_boundary || has_upper_boundary)) {
    dataindex = 0;
    int solindex = 0;
    bool waitformpi = false;
    MPI_Request request1[nx_sol];
    MPI_Request request2[nx_sol];
    for (int jx=mesh->xstart; jx<=mesh->xend; jx++) {
      if (!mesh->periodicY(jx)) {
	if (is_upper_boundary[jx]) {
	  MPI_Isend(&lower_transients_data[dataindex*transients_row_size], transients_row_size, MPI_DOUBLE, 0, NONLOCAL_PARALLEL_TAGBASE+jx,ycomms[jx-mesh->xstart],&request1[solindex]);
	  MPI_Irecv(&upper_transients_data[dataindex*transients_row_size], transients_row_size, MPI_DOUBLE, 0, NONLOCAL_PARALLEL_TAGBASE+mesh->LocalNx+jx,ycomms[jx-mesh->xstart],&request2[solindex]);
	  waitformpi = true;
	}
	if (is_lower_boundary[jx]) {
	  MPI_Isend(&upper_transients_data[dataindex*transients_row_size], transients_row_size, MPI_DOUBLE, ycomms_size[jx-mesh->xstart]-1, NONLOCAL_PARALLEL_TAGBASE+mesh->LocalNx+jx,ycomms[jx-mesh->xstart],&request1[solindex]);
	  MPI_Irecv(&lower_transients_data[dataindex*transients_row_size], transients_row_size, MPI_DOUBLE, ycomms_size[jx-mesh->xstart]-1, NONLOCAL_PARALLEL_TAGBASE+jx,ycomms[jx-mesh->xstart],&request2[solindex]);
	  waitformpi = true;
	}
	solindex++;
      }
      dataindex++;
    }
    if (waitformpi) {
      MPI_Waitall(nx_sol,request1,MPI_STATUSES_IGNORE);
      MPI_Waitall(nx_sol,request2,MPI_STATUSES_IGNORE);
    }
  }
  
  for (int j=0; j<number_of_negative_eigenvalues; j++) {
    if (!yperiodic && ( (nx_core>0) || (bc_apply_reflecting_sheath && nx_sol>0) ) ) {
      BoutReal shiftangle;
      FieldPerp decayfactor_core;
      FieldPerp decayfactor_sol;
      BoutReal* upper_transients_data_copy;
      
      if (has_upper_boundary && nx_core>0) {
        decayfactor_core = exp(sliceXZ(increasing_dimensionless_length, mesh->yend+1)/eigenvalues[j]); 
      }
      if (has_lower_boundary && nx_core>0) {
        decayfactor_core = exp(sliceXZ(decreasing_dimensionless_length, mesh->ystart-1)/eigenvalues[j]);
      }
      
      if (has_upper_boundary && nx_sol>0) {
        decayfactor_sol = exp(sliceXZ(increasing_dimensionless_length, mesh->yend)/eigenvalues[j]);
      }
      if (has_lower_boundary && nx_sol>0) {
        decayfactor_sol = exp(sliceXZ(decreasing_dimensionless_length, mesh->ystart)/eigenvalues[j]);
      }

      dataindex=0;
// {MPI_Barrier(BoutComm::get());Timer timer("toroidal");
      for (int jx=mesh->xstart; jx<=mesh->xend; jx++) {
	if (mesh->periodicY(jx,shiftangle)) {
	  toroidal_solver_reverse->solve(lower_transients_data+dataindex*transients_row_size+j*nz,&decayfactor_core(jx,0),jx,shiftangle);
	  dataindex++;
	}
	else if (bc_apply_reflecting_sheath) {
	  // Calculate lower transients on the upper boundary and vice versa in order to have the same communications to do as for the toroidal/periodic transients, for which this choice is more convenient/efficient
	  if (is_lower_boundary[jx]) {
	    if (is_upper_boundary[jx]) {
	      // Need to copy lower_transients_data before overwriting it to use in the upper_boundary_transients calculation
	      upper_transients_data_copy = new BoutReal[transients_row_size];
	      for (int i=0; i<transients_row_size; i++)
		upper_transients_data_copy[i] = lower_transients_data[(jx-mesh->ystart)*transients_row_size+i];
	    }
	    int i = 0;
	    for (int jz=0; jz<nz; jz++) {
	      i = dataindex*transients_row_size+j*nz+jz;
	      BoutReal thisdecayfactor = decayfactor_sol(jx,jz);
	      upper_transients_data[i] = (upper_transients_data[i]*thisdecayfactor + lower_transients_data[i]) / (1. - pow(thisdecayfactor,2));
	    }
	  }
	  if (is_upper_boundary[jx]) {
	    if (!is_lower_boundary[jx]) {
	      upper_transients_data_copy = upper_transients_data+(jx-mesh->ystart)*transients_row_size;
	    }
	    int i = 0;
	    for (int jz=0; jz<nz; jz++) {
	      // heatflux and friction need an extra minus sign because we do the calculation with the eigenvector coefficients factored in already and they have odd l values
	      i = dataindex*transients_row_size+j*nz+jz;
	      BoutReal thisdecayfactor = decayfactor_sol(jx,jz);
	      lower_transients_data[i] = (lower_transients_data[i]*thisdecayfactor + upper_transients_data_copy[i]) / (1. - pow(thisdecayfactor,2));
	    }
	    if (is_lower_boundary[jx]) {
	      delete [] upper_transients_data_copy;
	    }
	  }
	}
	dataindex++;
      }
    }
// MPI_Barrier(BoutComm::get());}output<<"toroidal solve time is "<<Timer::getTime("toroidal")<<endl;Timer::resetTime("toroidal");
    else if (yperiodic) {
      if (has_upper_boundary) {
	FieldPerp dividebythis = 1.0 - exp(sliceXZ(increasing_dimensionless_length, mesh->yend+1)/eigenvalues[j]); // Now calculate the actual transient using the stored total integral from zero
	for (int jx=mesh->xstart, i=0; jx<=mesh->xend; jx++, i++)
	  for (int jz=0; jz<nz; jz++) {
	    lower_transients_data[i*transients_row_size+j*nz+jz] /= dividebythis[jx][jz];
	  }
      }
    }
    
    if (!yperiodic) {
      BoutReal shiftangle;
      FieldPerp decayfactor;
      if (has_lower_boundary) {
        decayfactor = exp(sliceXZ(decreasing_dimensionless_length, mesh->ystart-1)/eigenvalues[j]);
      }
      dataindex=0;
// {MPI_Barrier(BoutComm::get());Timer timer("toroidal");
      for (int jx=mesh->xstart; jx<=mesh->xend; jx++)
	if (mesh->periodicY(jx,shiftangle)) {
	  toroidal_solver_forward->solve(upper_transients_data+dataindex*transients_row_size+j*nz,&decayfactor(jx,0),jx,shiftangle);
	  dataindex++;
	}
// MPI_Barrier(BoutComm::get());}output<<"toroidal solve time is "<<Timer::getTime("toroidal")<<endl;Timer::resetTime("toroidal");
    }
    else {
      if (has_lower_boundary) {
	FieldPerp dividebythis = 1.0 - exp(sliceXZ(decreasing_dimensionless_length, mesh->ystart-1)/eigenvalues[j]); // Now calculate the actual transient using the stored total integral from zero
	for (int jx=mesh->xstart, i=0; jx<=mesh->xend; jx++, i++)
	  for (int jz=0; jz<nz; jz++) {
	    upper_transients_data[i*transients_row_size+j*nz+jz] /= dividebythis[jx][jz];
	  }
      }
    }
  }
  
  if ( yperiodic || (nx_core>0) || (bc_apply_reflecting_sheath && nx_sol>0) ) {
    dataindex = 0;
    for (int jx=mesh->xstart; jx<=mesh->xend; jx++) {
      if (yperiodic || mesh->periodicY(jx) || bc_apply_reflecting_sheath) {
	y_broadcast(&lower_transients_data[dataindex*transients_row_size], transients_row_size, ycomms_size[jx-mesh->xstart]-1, jx);
	y_broadcast(&upper_transients_data[dataindex*transients_row_size], transients_row_size, 0, jx);
	dataindex++;
      }
    }
  }
  
  if (!yperiodic && (nx_sol>0)) {
    for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
      for (int jz=0; jz<nz; jz++) {
	position->jx=rup.ind;
	position->jy=mesh->yend;
	position->jz=jz;
	calc_index(position);
	BoutReal Te_here = interp_to_point_YLOW(T_electron,*position);
	if (bc_heatflux) {
	  pass_interim_upper_boundary_n11(rup.ind,jz) = electron_heat_flux(rup.ind,mesh->yend,jz);
	  if (bc_apply_reflecting_sheath)
	    for (int j=0; j<number_of_negative_eigenvalues; j++)
	      pass_interim_upper_boundary_n11(rup.ind,jz) += heatflux_coefficients_above[j]*upper_transients_data[(rup.ind-mesh->xstart)*transients_row_size+j*nz+jz]
                + heatflux_coefficients_below[j]*lower_transients_data[(rup.ind-mesh->xstart)*transients_row_size+j*nz+jz]*exp(increasing_dimensionless_length(rup.ind,mesh->yend,jz)/eigenvalues[j]);
	  upper_boundary_condition_n11(rup.ind,jz) = -4./5.*sqrt(electron_mass/2.)/pow(Te_here,1.5)*heat_flux_boundary_condition(rup.ind,mesh->yend,jz);
	}
	if (bc_viscosity) {
	  pass_interim_upper_boundary_n20(rup.ind,jz) = electron_viscosity(rup.ind,mesh->yend,jz);
	  if (bc_apply_reflecting_sheath)
	    for (int j=0; j<number_of_negative_eigenvalues; j++)
	      pass_interim_upper_boundary_n20 += viscosity_coefficients[j]*( upper_transients_data[(rup.ind-mesh->xstart)*transients_row_size+j*nz+jz]
                                                                             + lower_transients_data[(rup.ind-mesh->xstart)*transients_row_size+j*nz+jz]*exp(increasing_dimensionless_length(rup.ind,mesh->yend,jz)/eigenvalues[j]) );
	  upper_boundary_condition_n20(rup.ind,jz) = viscosity_boundary_condition(rup.ind,mesh->yend,jz)/Te_here;
	}
      }
    if (has_upper_boundary && nx_sol>0 && mesh->getNYPE()>1) {
      MPI_Request request1[nx_sol];
      MPI_Request request2[nx_sol];
      MPI_Request request3[nx_sol];
      MPI_Request request4[nx_sol];
      int i=0;
      for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++) {
	if (rup.ind<mesh->xstart || rup.ind>mesh->xend) continue;
	if (bc_heatflux) {
	  MPI_Isend(&pass_interim_upper_boundary_n11(rup.ind,0), // get a pointer to the jx'th row of FieldPerp's data
		    nz,
		    MPI_DOUBLE,
		    0,
		    NONLOCAL_PARALLEL_TAGBASE,
		    ycomms[rup.ind-mesh->xstart],
		    &request1[i]);
	  MPI_Isend(&upper_boundary_condition_n11(rup.ind,0), // get a pointer to the jx'th row of FieldPerp's data
		    nz,
		    MPI_DOUBLE,
		    0,
		    NONLOCAL_PARALLEL_TAGBASE,
		    ycomms[rup.ind-mesh->xstart],
		    &request2[i]);
	}
	if (bc_viscosity) {
	  MPI_Isend(&pass_interim_upper_boundary_n20(rup.ind,0), // get a pointer to the jx'th row of FieldPerp's data
		    nz,
		    MPI_DOUBLE,
		    0,
		    NONLOCAL_PARALLEL_TAGBASE,
		    ycomms[rup.ind-mesh->xstart],
		    &request3[i]);
	  MPI_Isend(&upper_boundary_condition_n20(rup.ind,0), // get a pointer to the jx'th row of FieldPerp's data
		    nz,
		    MPI_DOUBLE,
		    0,
		    NONLOCAL_PARALLEL_TAGBASE,
		    ycomms[rup.ind-mesh->xstart],
		    &request4[i]);
	}
	i++;
      }
      if (bc_heatflux) {
	MPI_Waitall(nx_sol,request1,MPI_STATUSES_IGNORE);
	MPI_Waitall(nx_sol,request2,MPI_STATUSES_IGNORE);
      }
      if (bc_viscosity) {
	MPI_Waitall(nx_sol,request3,MPI_STATUSES_IGNORE);
	MPI_Waitall(nx_sol,request4,MPI_STATUSES_IGNORE);
      }
    }
    if (has_lower_boundary && nx_sol>0 && mesh->getNYPE()>1) {
      MPI_Request request1[nx_sol];
      MPI_Request request2[nx_sol];
      MPI_Request request3[nx_sol];
      MPI_Request request4[nx_sol];
      int i=0;
      for (RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++) {
	if (rlow.ind<mesh->xstart || rlow.ind>mesh->xend) continue;
	if (bc_heatflux) {
	  MPI_Irecv(&pass_interim_upper_boundary_n11(rlow.ind,0), // get a pointer to the jx'th row of FieldPerp's data
		    nz,
		    MPI_DOUBLE,
		    ycomms_size[rlow.ind-mesh->xstart]-1,
		    NONLOCAL_PARALLEL_TAGBASE,
		    ycomms[rlow.ind-mesh->xstart],
		    &request1[i]);
	  MPI_Irecv(&upper_boundary_condition_n11(rlow.ind,0), // get a pointer to the jx'th row of FieldPerp's data
		    nz,
		    MPI_DOUBLE,
		    ycomms_size[rlow.ind-mesh->xstart]-1,
		    NONLOCAL_PARALLEL_TAGBASE,
		    ycomms[rlow.ind-mesh->xstart],
		    &request2[i]);
	}
	if (bc_viscosity) {
	  MPI_Irecv(&pass_interim_upper_boundary_n20(rlow.ind,0), // get a pointer to the jx'th row of FieldPerp's data
		    nz,
		    MPI_DOUBLE,
		    ycomms_size[rlow.ind-mesh->xstart]-1,
		    NONLOCAL_PARALLEL_TAGBASE,
		    ycomms[rlow.ind-mesh->xstart],
		    &request3[i]);
	  MPI_Irecv(&upper_boundary_condition_n20(rlow.ind,0), // get a pointer to the jx'th row of FieldPerp's data
		    nz,
		    MPI_DOUBLE,
		    ycomms_size[rlow.ind-mesh->xstart]-1,
		    NONLOCAL_PARALLEL_TAGBASE,
		    ycomms[rlow.ind-mesh->xstart],
		    &request4[i]);
	}
	i++;
      }
      if (bc_heatflux) {
	MPI_Waitall(nx_sol,request1,MPI_STATUSES_IGNORE);
	MPI_Waitall(nx_sol,request2,MPI_STATUSES_IGNORE);
      }
      if (bc_viscosity) {
	MPI_Waitall(nx_sol,request3,MPI_STATUSES_IGNORE);
	MPI_Waitall(nx_sol,request4,MPI_STATUSES_IGNORE);
      }
    }
    
    if (has_lower_boundary && nx_sol>0) {
      for (int jx=mesh->xstart, i=0; jx<=mesh->xend; jx++) {
	if (!mesh->periodicY(jx)) {
	  for (int jz=0; jz<nz; jz++) {
	    position->jx=jx;
	    position->jy=mesh->ystart;
	    position->jz=jz;
	    calc_index(position);
	    BoutReal Te_here = interp_to_point_YLOW(T_electron,*position);
	    BoutReal lower_boundary_n11;
	    BoutReal upper_boundary_n11;
	    BoutReal interim_lower_boundary_n11;
	    BoutReal interim_upper_boundary_n11;
	    BoutReal lower_boundary_n20;
	    BoutReal upper_boundary_n20;
	    BoutReal interim_lower_boundary_n20;
	    BoutReal interim_upper_boundary_n20;
	    if (bc_heatflux) {
	      lower_boundary_n11 = -4./5.*sqrt(electron_mass/2.)/pow(Te_here,1.5)*heat_flux_boundary_condition(jx,mesh->ystart,jz);
	      upper_boundary_n11 = upper_boundary_condition_n11(jx,jz);
	      interim_lower_boundary_n11 = electron_heat_flux(jx,mesh->ystart,jz);
	      if (bc_apply_reflecting_sheath)
		for (int j=0; j<number_of_negative_eigenvalues; j++)
		  interim_lower_boundary_n11 += heatflux_coefficients_below[j]*lower_transients_data[(jx-mesh->xstart)*transients_row_size+j*nz+jz]
                    + heatflux_coefficients_above[j]*upper_transients_data[(jx-mesh->xstart)*transients_row_size+j*nz+jz]*exp(decreasing_dimensionless_length(jx,mesh->ystart,jz)/eigenvalues[j]);
	      interim_upper_boundary_n11 = pass_interim_upper_boundary_n11(jx,jz);
	    }
	    if (bc_viscosity) {
	      lower_boundary_n20 = viscosity_boundary_condition(jx,mesh->ystart,jz)/Te_here;
	      upper_boundary_n20 = upper_boundary_condition_n20(jx,jz);
	      interim_lower_boundary_n20 = electron_viscosity(jx,mesh->ystart,jz);
	      if (bc_apply_reflecting_sheath)
		for (int j=0; j<number_of_negative_eigenvalues; j++)
		  interim_lower_boundary_n20 += viscosity_coefficients[j]*( lower_transients_data[(jx-mesh->xstart)*transients_row_size+j*nz+jz]
									    + upper_transients_data[(jx-mesh->xstart)*transients_row_size+j*nz+jz]*exp(decreasing_dimensionless_length(jx,mesh->ystart,jz)/eigenvalues[j]) );
	      interim_upper_boundary_n20 = pass_interim_upper_boundary_n20[jx][jz];
	    }
	    if (bc_heatflux && !bc_viscosity) {
	      /*
	      electron_heat_flux is, at this point, the contribution to n11 from nhat_plus.
	      We want the actual heat flux at mesh->ystart to be boundary_heat_flux.
	      Thus the remainder must come from nhat_minus, which we will construct here just to give the right 1,1 component (could set number_of_negative_eigenvalues-1 more components if desired)
	      However, the transients from the other boundary do not necessarily decay to vanishing by the time they get to this boundary, so we must solve for both.
	      Fortunately this reduces to a single algebraic equation which makes it surprisingly easy to do (in the case of just a single condition being imposed at least).
	      */
	      
	      BoutReal sum_decayed_W11_W11_term = 0.;
	      for (int j=0; j<number_of_negative_eigenvalues; j++) {
		exp_total_dimensionless_length_over_eigenvalue[j] = exp(total_dimensionless_length[jx][jz]/eigenvalues[j]);
		sum_decayed_W11_W11_term += heatflux_coefficients_below[j]*WinverseB_11[j]*exp_total_dimensionless_length_over_eigenvalue[j];
	      }
	      heatflux_sol_transients_factors[i*(2*nz)+jz] = ( (lower_boundary_n11 - interim_lower_boundary_n11)*W11_dot_W11
									- sum_decayed_W11_W11_term*(upper_boundary_n11 - interim_upper_boundary_n11) )
									/ ( pow(W11_dot_W11,2) - pow(sum_decayed_W11_W11_term,2) );
	      heatflux_sol_transients_factors[i*(2*nz)+nz+jz] = ( (upper_boundary_n11 - interim_upper_boundary_n11)*W11_dot_W11
														- sum_decayed_W11_W11_term*(lower_boundary_n11 - interim_lower_boundary_n11) )
														/ ( pow(W11_dot_W11,2) - pow(sum_decayed_W11_W11_term,2) );
	    }
	    if (bc_viscosity && !bc_heatflux) {
	      /*
	      electron_viscosity is, at this point, the contribution to n20 from nhat_plus.
	      We want the actual viscosity at mesh->ystart to be boundary_viscosity.
	      Thus the remainder must come from nhat_minus, which we will construct here just to give the right 2,0 component (could set number_of_negative_eigenvalues-1 more components if desired)
	      However, the transients from the other boundary do not necessarily decay to vanishing by the time they get to this boundary, so we must solve for both.
	      Fortunately this reduces to a single algebraic equation which makes it surprisingly easy to do (in the case of just a single condition being imposed at least).
	      */
	      
	      BoutReal sum_decayed_W20_W20_term = 0.;
	      for (int j=0; j<number_of_negative_eigenvalues; j++) {
		exp_total_dimensionless_length_over_eigenvalue[j] = exp(total_dimensionless_length(jx,jz)/eigenvalues[j]);
		sum_decayed_W20_W20_term += viscosity_coefficients[j]*WinverseB_20[j]*exp_total_dimensionless_length_over_eigenvalue[j];
	      }
	      viscosity_sol_transients_factors[i*(2*nz)+jz] = ( (lower_boundary_n20 - interim_lower_boundary_n20)*W20_dot_W20
									- sum_decayed_W20_W20_term*(upper_boundary_n20 - interim_upper_boundary_n20) )
									/ ( pow(W20_dot_W20,2) - pow(sum_decayed_W20_W20_term,2) );
	      viscosity_sol_transients_factors[i*(2*nz)+nz+jz] = ( (upper_boundary_n20 - interim_upper_boundary_n20)*W20_dot_W20
														- sum_decayed_W20_W20_term*(lower_boundary_n20 - interim_lower_boundary_n20) )
														/ ( pow(W20_dot_W20,2) - pow(sum_decayed_W20_W20_term,2));
	    }
	    if (bc_heatflux && bc_viscosity) {
	      BoutReal sum_decayed_W11_W11_term = 0.;
	      BoutReal sum_decayed_W20_W20_term = 0.;
	      BoutReal sum_decayed_W11_W20_term = 0.;
	      BoutReal sum_decayed_W20_W11_term = 0.;
	      for (int j=0; j<number_of_negative_eigenvalues; j++) {
		exp_total_dimensionless_length_over_eigenvalue[j] = exp(total_dimensionless_length(jx,jz)/eigenvalues[j]);
		sum_decayed_W11_W11_term += heatflux_coefficients_below[j]*WinverseB_11[j]*exp_total_dimensionless_length_over_eigenvalue[j];
		sum_decayed_W20_W20_term += viscosity_coefficients[j]*WinverseB_20[j]*exp_total_dimensionless_length_over_eigenvalue[j];
		sum_decayed_W20_W11_term += viscosity_coefficients[j]*WinverseB_11[j]*exp_total_dimensionless_length_over_eigenvalue[j];
		sum_decayed_W11_W20_term += heatflux_coefficients_below[j]*WinverseB_20[j]*exp_total_dimensionless_length_over_eigenvalue[j];
	      }
	      BoutReal det = pow(sum_decayed_W11_W20_term,2)*pow(sum_decayed_W20_W11_term,2)
			    - 2*sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
			    + pow(sum_decayed_W11_W11_term,2)*pow(sum_decayed_W20_W20_term,2)
			    - pow(sum_decayed_W20_W20_term,2)*pow(W11_dot_W11,2)
			    + 2*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W11*W11_dot_W20
			    - pow(sum_decayed_W20_W11_term,2)*pow(W11_dot_W20,2)
			    - 2*sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W11_dot_W11*W20_dot_W11
			    + 2*sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W11_dot_W20*W20_dot_W11
			    - pow(sum_decayed_W11_W20_term,2)*pow(W20_dot_W11,2)
			    + pow(W11_dot_W20,2)*pow(W20_dot_W11,2)
			    + 2*sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W11_dot_W11*W20_dot_W20
			    - 2*sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W11_dot_W20*W20_dot_W20
			    + 2*sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W11*W20_dot_W20
			    - 2*W11_dot_W11*W11_dot_W20*W20_dot_W11*W20_dot_W20
			    - pow(sum_decayed_W11_W11_term,2)*pow(W20_dot_W20,2)
			    + pow(W11_dot_W11,2)*pow(W20_dot_W20,2);
	      heatflux_sol_transients_factors[i*(2*nz)+jz] = ( (upper_boundary_n11-interim_upper_boundary_n11)*( -sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
																    + sum_decayed_W11_W11_term*pow(sum_decayed_W20_W20_term,2)
																    + sum_decayed_W20_W20_term*W11_dot_W20*W20_dot_W11
																    - sum_decayed_W20_W11_term*W11_dot_W20*W20_dot_W20
																    + sum_decayed_W11_W20_term*W20_dot_W11*W20_dot_W20
																    - sum_decayed_W11_W11_term*pow(W20_dot_W20,2) )
										    + (upper_boundary_n20-interim_upper_boundary_n20)*( pow(sum_decayed_W11_W20_term,2)*sum_decayed_W20_W11_term
																	- sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W20_term
																	+ sum_decayed_W20_W20_term*W11_dot_W11*W11_dot_W20
																	- sum_decayed_W20_W11_term*pow(W11_dot_W20,2)
																	+ sum_decayed_W11_W20_term*W11_dot_W11*W20_dot_W20
																	- sum_decayed_W11_W11_term*W11_dot_W20*W20_dot_W20 )
										    + (lower_boundary_n11-interim_lower_boundary_n11)*( -pow(sum_decayed_W20_W20_term,2)*W11_dot_W11
																	+ sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W20
																	- sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W20_dot_W11
																	+ sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W20_dot_W20
																	- W11_dot_W20*W20_dot_W11*W20_dot_W20 + W11_dot_W11*pow(W20_dot_W20,2) )
										    + (lower_boundary_n20-interim_lower_boundary_n20)*( -sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W11_dot_W11
																	+ sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W11_dot_W20
																	- pow(sum_decayed_W11_W20_term,2)*W20_dot_W11
																	+ pow(W11_dot_W20,2)*W20_dot_W11
																	+ sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W20
																	- W11_dot_W11*W11_dot_W20*W20_dot_W20 )
										) / det;
	      heatflux_sol_transients_factors[i*(2*nz)+nz+jz] = ( (upper_boundary_n11-interim_upper_boundary_n11)*( -pow(sum_decayed_W20_W20_term,2)*W11_dot_W11
																					      + sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W20
																					      - sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W20_dot_W11
																					      + sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W20_dot_W20
																					      - W11_dot_W20*W20_dot_W11*W20_dot_W20
																					      + W11_dot_W11*pow(W20_dot_W20,2) )
															    + (upper_boundary_n20-interim_upper_boundary_n20)*( sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W11_dot_W11
																						- sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W11_dot_W20
																						+ pow(sum_decayed_W11_W20_term,2)*W20_dot_W11
																						- pow(W11_dot_W20,2)*W20_dot_W11
																						- sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W20
																						+ W11_dot_W11*W11_dot_W20*W20_dot_W20 )
															    + (lower_boundary_n11-interim_lower_boundary_n11)*( -sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
																						+ sum_decayed_W11_W11_term*pow(sum_decayed_W20_W20_term,2)
																						+ sum_decayed_W20_W20_term*W11_dot_W20*W20_dot_W11
																						- sum_decayed_W20_W11_term*W11_dot_W20*W20_dot_W20
																						+ sum_decayed_W11_W20_term*W20_dot_W11*W20_dot_W20
																						- sum_decayed_W11_W11_term*pow(W20_dot_W20,2) )
															    + (lower_boundary_n20-interim_lower_boundary_n20)*( -pow(sum_decayed_W11_W20_term,2)*sum_decayed_W20_W11_term
																						+ sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W20_term
																						- sum_decayed_W20_W20_term*W11_dot_W11*W11_dot_W20
																						+ sum_decayed_W20_W11_term*pow(W11_dot_W20,2)
																						- sum_decayed_W11_W20_term*W11_dot_W11*W20_dot_W20
																						+ sum_decayed_W11_W11_term*W11_dot_W20*W20_dot_W20 )
															  ) / det;
	      viscosity_sol_transients_factors[i*(2*nz)+jz] = ( (upper_boundary_n11-interim_upper_boundary_n11)*( sum_decayed_W11_W20_term*pow(sum_decayed_W20_W11_term,2)
																    - sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
																    - sum_decayed_W20_W20_term*W11_dot_W11*W20_dot_W11
																    - sum_decayed_W11_W20_term*pow(W20_dot_W11,2)
																    + sum_decayed_W20_W11_term*W11_dot_W11*W20_dot_W20
																    + sum_decayed_W11_W11_term*W20_dot_W11*W20_dot_W20 )
										    + (upper_boundary_n20-interim_upper_boundary_n20)*( -sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W11_term
																	+ pow(sum_decayed_W11_W11_term,2)*sum_decayed_W20_W20_term
																	- sum_decayed_W20_W20_term*pow(W11_dot_W11,2)
																	+ sum_decayed_W20_W11_term*W11_dot_W11*W11_dot_W20
																	- sum_decayed_W11_W20_term*W11_dot_W11*W20_dot_W11
																	+ sum_decayed_W11_W11_term*W11_dot_W20*W20_dot_W11 )
										    + (lower_boundary_n11-interim_lower_boundary_n11)*( sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W11
																	- pow(sum_decayed_W20_W11_term,2)*W11_dot_W20
																	+ sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W20_dot_W11
																	+ W11_dot_W20*pow(W20_dot_W11,2)
																	- sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W20_dot_W20
																	- W11_dot_W11*W20_dot_W11*W20_dot_W20 )
										    + (lower_boundary_n20-interim_lower_boundary_n20)*( sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W11_dot_W11
																	- sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W11_dot_W20
																	+ sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W11
																	- W11_dot_W11*W11_dot_W20*W20_dot_W11
																	- pow(sum_decayed_W11_W11_term,2)*W20_dot_W20
																	+ pow(W11_dot_W11,2)*W20_dot_W20 )
										) / det;
	      viscosity_sol_transients_factors[i*(2*nz)+nz+jz] = ( (upper_boundary_n11-interim_upper_boundary_n11)*( -sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W11
																					      + pow(sum_decayed_W20_W11_term,2)*W11_dot_W20
																					      - sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W20_dot_W11
																					      - W11_dot_W20*pow(W20_dot_W11,2)
																					      + sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W20_dot_W20
																					      + W11_dot_W11*W20_dot_W11*W20_dot_W20 )
															    + (upper_boundary_n20-interim_upper_boundary_n20)*( sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W11_dot_W11
																						- sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W11_dot_W20
																						+ sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W11
																						- W11_dot_W11*W11_dot_W20*W20_dot_W11
																						- pow(sum_decayed_W11_W11_term,2)*W20_dot_W20
																						+ pow(W11_dot_W11,2)*W20_dot_W20 )
															    + (lower_boundary_n11-interim_lower_boundary_n11)*( -sum_decayed_W11_W20_term*pow(sum_decayed_W20_W11_term,2)
																						+ sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
																						+ sum_decayed_W20_W20_term*W11_dot_W11*W20_dot_W11
																						+ sum_decayed_W11_W20_term*pow(W20_dot_W11,2)
																						- sum_decayed_W20_W11_term*W11_dot_W11*W20_dot_W20
																						- sum_decayed_W11_W11_term*W20_dot_W11*W20_dot_W20 )
															    + (lower_boundary_n20-interim_lower_boundary_n20)*( -sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W11_term
																						+ pow(sum_decayed_W11_W11_term,2)*sum_decayed_W20_W20_term
																						- sum_decayed_W20_W20_term*pow(W11_dot_W11,2)
																						+ sum_decayed_W20_W11_term*W11_dot_W11*W11_dot_W20
																						- sum_decayed_W11_W20_term*W11_dot_W11*W20_dot_W11
																						+ sum_decayed_W11_W11_term*W11_dot_W20*W20_dot_W11 )
															  ) / det;
	    }
	  }
	  i++;
	}
      }
    }
  }

  if (bc_heatflux && (nx_sol>0)) {
    for (int jx=mesh->xstart, i=0; jx<=mesh->xend; jx++) {
      if (!mesh->periodicY(jx)) {
	y_broadcast(&heatflux_sol_transients_factors[i*(2*nz)], (nz)*2, 0, jx);
	i++;
      }
    }
  }
  if (bc_viscosity && (nx_sol>0)) {
    for (int jx=mesh->xstart, i=0; jx<=mesh->xend; jx++) {
      if (!mesh->periodicY(jx)) {
	y_broadcast(&viscosity_sol_transients_factors[i*(2*nz)], (nz)*2, 0, jx);
	i++;
      }
    }
  }
  
  if (nx_sol>0) {
    for (int jx=mesh->xstart, i=0; jx<=mesh->xend; jx++, i++) {
      if (!mesh->periodicY(jx)) {
	for (int jz=0; jz<nz; jz++) {
	  if (bc_heatflux) {
	    for (int j=0; j<number_of_negative_eigenvalues; j++) {
	      lower_transients_data[i*transients_row_size+j*nz+jz] += WinverseB_11[j]*heatflux_sol_transients_factors[i*(2*nz)+jz];
	      upper_transients_data[i*transients_row_size+j*nz+jz] += -WinverseB_11[j]*heatflux_sol_transients_factors[i*(2*nz)+nz+jz];
	    }
	  }
	  if (bc_viscosity) {
	    for (int j=0; j<number_of_negative_eigenvalues; j++) {
	      lower_transients_data[i*transients_row_size+j*nz+jz] += WinverseB_20[j]*viscosity_sol_transients_factors[i*(2*nz)+jz];
	      upper_transients_data[i*transients_row_size+j*nz+jz] += WinverseB_20[j]*viscosity_sol_transients_factors[i*(2*nz)+nz+jz];
	    }
	  }
	}
      }
    }
  }

  #ifndef NOEDGETERMS
  start_index(position);
  do {
    position->jy=mesh->ystart;
    calc_index(position);
    do {
      for (int j=0; j<number_of_negative_eigenvalues; j++) {
	BoutReal lower_transient = lower_transients_data[(position->jx-mesh->ystart)*transients_row_size+j*nz+position->jz] * exp(increasing_dimensionless_length[*position]/eigenvalues[j]);
	BoutReal upper_transient = upper_transients_data[(position->jx-mesh->ystart)*transients_row_size+j*nz+position->jz] * exp(decreasing_dimensionless_length[*position]/eigenvalues[j]);
	if (calculate_heatflux) {
	  electron_heat_flux[*position] += heatflux_coefficients_below[j]*lower_transient + heatflux_coefficients_above[j]*upper_transient;
// output<<position->jx<<","<<position->jy<<","<<position->jz<<" "<<heatflux_lower_boundary_transients[j][position->jx][position->jz]<<" "<<heatflux_upper_boundary_transients[j][position->jx][position->jz]<<endl;
	}
	if (calculate_viscosity) {
	  electron_viscosity[*position] += viscosity_coefficients[j]*(lower_transient + upper_transient);
	}
	if (calculate_friction) {
	  electron_friction[*position] += friction_coefficients_below[j]*lower_transient + friction_coefficients_above[j]*upper_transient;
	}
      }
      position->jy++;
      calc_index(position);
    } while (position->jy<mesh->yend+1);
  } while (next_indexperp(position));
  #endif
// output<<"with edgeterms "<<electron_friction[40][35][0]<<" "<<electron_friction[40][35][0]-temp_ef<<endl;
// temp_ef=electron_friction[40][35][0];
  
  if (calculate_heatflux) {
// for (int x=0; x<mesh->LocalNx; x++)for (int y=0; y<mesh->LocalNy; y++){
//   output<<x<<","<<y<<"   ";
//   for (int z=0; z<mesh->LocalNz; z++) output<<electron_heat_flux[x][y][z]<<" ";
//   output<<endl;
// }
    electron_heat_flux *= -5./4.*sqrt(2./electron_mass)*pow(T_electron,1.5); //now we have q=-5/4*v_Telectron*T_electron*n^(1,1)
    mesh->communicate(electron_heat_flux);
  }
  if (calculate_viscosity) {
    electron_viscosity *= T_electron;
    mesh->communicate(electron_viscosity);
  }
  if (calculate_friction) {
    electron_friction *= 2.*T_electron*lambdaC_inverse;
//     electron_friction *= 0.;
// output<<"normalized "<<electron_friction[40][35][0]<<endl;
// temp_ef=electron_friction[40][35][0];

    electron_friction += -sqrt(2.*electron_mass*T_electron)*lambdaC_inverse*(-j_parallel/elementary_charge); // Need to include also the friction due directly to the Maxwellian part of the distribution function
// output<<"with Maxwellian contrib "<<electron_friction[40][35][0]<<" "<<electron_friction[40][35][0]-temp_ef<<endl;
// temp_ef=electron_friction[40][35][0];
    mesh->communicate(electron_friction);
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NonLocalParallel::calculate_nonlocal_closures_cell_ylow() {
  int dataindex=0; // Use for cycling through jx indices of lower/upper_transients_data
  
  lambdaC_inverse = n_electron * pow(elementary_charge,4) * logLambda / 12. / pow(PI,1.5) / pow(epsilon_0,2) / SQ(T_electron);
  if (normalized_density) lambdaC_inverse *= density_normalization; // If set_density_normalization has been called, un-normalize the density used to calculate the collision length so that it has the correct units
  
  if (gradT_drive) {
    gradT_driveterm = 5./4. * n_electron / T_electron / lambdaC_inverse; //g^(1,1)
    gradT_electron = Grad_par(T_electron,CELL_YLOW);
  }
  if (gradV_drive) {
    gradV_driveterm = -0.5 * n_electron / sqrt(2.*T_electron/electron_mass) * Grad_par(V_electron,CELL_CENTRE) / lambdaC_inverse;
//     gradV_electron = Grad_par(V_electron,CELL_YLOW); // would be more consistent to put this on CELL_CENTRE, but that would make imposing a boundary condition a pain.
  }
  if (VeminusVi_drive) {
    VeminusVi_driveterm = -2./sqrt(PI) * (-1./elementary_charge) / sqrt(2.*T_electron/electron_mass);
    // j_parallel is already CELL_YLOW.
  }
  // Calculate target boundary guard cell derivitives (at YLOW) for gradT_electron with 4th order forward/backward differences from T_electron (at CENTRE)
  // Also check for unphysical lambdaC_inverse
  for (RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++)
    for (int jz=0; jz<nz; jz++) {
      for (int jy=mesh->ystart-1; jy>=0; jy--) {
	if (gradT_drive) {
	  gradT_electron(rlow.ind,jy,jz)=(-93.*T_electron(rlow.ind,jy,jz) + 229.*T_electron(rlow.ind,jy+1,jz) - 225.*T_electron(rlow.ind,jy+2,jz) + 111.*T_electron(rlow.ind,jy+3,jz) - 22.*T_electron(rlow.ind,jy+4,jz))/48./coord->dy(rlow.ind,jy)/sqrt((coord->g_22(rlow.ind,jy) + coord->g_22(rlow.ind,jy+1) + coord->g_22(rlow.ind,jy+2) + coord->g_22(rlow.ind,jy+3) + coord->g_22(rlow.ind,jy+4))/5.);
	  if (abs(gradT_driveterm(rlow.ind,jy,jz))>1.e37 || gradT_driveterm(rlow.ind,jy,jz)!=gradT_driveterm(rlow.ind,jy,jz)) gradT_driveterm(rlow.ind,jy,jz) = 1.e37;
	}
	// Nothing to be done for gradV_drive: gradV_driveterm is CELL_CENTRE and the guard cell values are not used
	if (VeminusVi_drive) {
	  if (abs(VeminusVi_driveterm(rlow.ind,jy,jz))>1.e37 || VeminusVi_driveterm(rlow.ind,jy,jz)!=VeminusVi_driveterm(rlow.ind,jy,jz)) VeminusVi_driveterm(rlow.ind,jy,jz) = 1.e37;
	}
	if (abs(lambdaC_inverse(rlow.ind,jy,jz))>1.e37 || lambdaC_inverse(rlow.ind,jy,jz)!=lambdaC_inverse(rlow.ind,jy,jz)) lambdaC_inverse(rlow.ind,jy,jz) = 1.e37;
      }
    }
  for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
    for (int jz=0; jz<nz; jz++) {
      if (gradT_drive) {
	gradT_electron(rup.ind,mesh->yend,jz) = (T_electron(rup.ind,mesh->yend-2,jz)-27.*T_electron(rup.ind,mesh->yend-1,jz)+27.*T_electron(rup.ind,mesh->yend,jz)-T_electron(rup.ind,mesh->yend+1,jz))/24./coord->dy(rup.ind,mesh->yend+1)/sqrt((coord->g_22(rup.ind,mesh->yend-1) + coord->g_22(rup.ind,mesh->yend) + coord->g_22(rup.ind,mesh->yend+1) + coord->g_22(rup.ind,mesh->yend+2))/4.);
      }
      for (int jy=mesh->yend+1; jy<mesh->LocalNy; jy++) {
	if (gradT_drive) {
	  gradT_electron(rup.ind,jy,jz)=(93.*T_electron(rup.ind,jy-1,jz) - 229.*T_electron(rup.ind,jy-2,jz) + 225.*T_electron(rup.ind,jy-3,jz) - 111.*T_electron(rup.ind,jy-4,jz) + 22.*T_electron(rup.ind,jy-5,jz))/48./coord->dy(rup.ind,jy-1)/sqrt((coord->g_22(rup.ind,jy-1) + coord->g_22(rup.ind,jy-2) + coord->g_22(rup.ind,jy-3) + coord->g_22(rup.ind,jy-4) + coord->g_22(rup.ind,jy-5))/5.);
	  if (abs(gradT_driveterm(rup.ind,jy,jz))>1.e37 || gradT_driveterm(rup.ind,jy,jz)!=gradT_driveterm(rup.ind,jy,jz)) gradT_driveterm(rup.ind,jy,jz) = 1.e37;
	}
	// Nothing to be done for gradV_drive: gradV_driveterm is CELL_CENTRE and the guard cell values are not used
	if (VeminusVi_drive) {
	  if (abs(VeminusVi_driveterm(rup.ind,jy,jz))>1.e37 || VeminusVi_driveterm(rup.ind,jy,jz)!=VeminusVi_driveterm(rup.ind,jy,jz)) VeminusVi_driveterm(rup.ind,jy,jz) = 1.e37;
	}
	if (abs(lambdaC_inverse(rup.ind,jy,jz))>1.e37 || lambdaC_inverse(rup.ind,jy,jz)!=lambdaC_inverse(rup.ind,jy,jz)) lambdaC_inverse(rup.ind,jy,jz) = 1.e37;
      }
  }
  
  if (gradT_drive) {
    mesh->communicate(gradT_electron);
  }
  if (gradV_drive) {
    mesh->communicate(gradV_driveterm);
  }
  // No gradient term to communicate for VeminusVi
  
  // Now calculate z and deltaz everywhere
  cubic_spline_inverse_lambdaC.calculate(lambdaC_inverse);
  
  start_index(position);
  
  do {
    
    if (!is_lower_boundary[position->jx]) {
      Timer timer("comms");
      if (position->jx < mesh->DownXSplitIndex()) {
        mesh->wait(mesh->irecvYInIndest(&increasing_dimensionless_length(position->jx,mesh->ystart,position->jz),1,
                                        NONLOCAL_PARALLEL_TAGBASE + position->jx*mesh->LocalNz+position->jz));
      } else {
        mesh->wait(mesh->irecvYInOutdest(&increasing_dimensionless_length(position->jx,mesh->ystart,position->jz),1,
                                         NONLOCAL_PARALLEL_TAGBASE + position->jx*mesh->LocalNz+position->jz));
      }
    } else {
      position->jy = mesh->ystart-2;
      calc_index(position);
      interp_coefficients = cubic_spline_inverse_lambdaC.coefficients(position);
      dimensionless_length_deltas_below(position->jx,position->jyp,position->jz) = coord->dy(position->jx,position->jy)*sqrt(coord->g_22(position->jx,position->jyp))
										      *(interp_coefficients[0]/2. + interp_coefficients[1]/2.*3./4. + interp_coefficients[2]/3.*7./8. + interp_coefficients[3]/4.*15./16.);
    }
    
    position->jy = mesh->ystart-1;
    calc_index(position);
    
    interp_coefficients = cubic_spline_inverse_lambdaC.coefficients(position);
    dimensionless_length_deltas_below(position->jx,position->jyp,position->jz) = coord->dy(position->jx,position->jy)*sqrt(coord->g_22(position->jx,position->jyp))
      *(interp_coefficients[0]/2. + interp_coefficients[1]/2.*3./4. + interp_coefficients[2]/3.*7./8. + interp_coefficients[3]/4.*15./16.);
    interp_coefficients = cubic_spline_inverse_lambdaC.coefficients(position);
    // dimensionless_length_deltas_above[jy] and dimensionless_length_deltas_below[jy] are the deltaz's for the half-step above and below, respectively, the CELL_CENTRE at jy
    // deltaz between position[jy](CELL_CENTRE), where t=0, and position[jyp](CELL_YLOW), where t=0.5
    dimensionless_length_deltas_above(position->jx,position->jy,position->jz) = coord->dy(position->jx,position->jy)*sqrt(coord->g_22(position->jx,position->jy))
      *(interp_coefficients[0]/2. + interp_coefficients[1]/2./4. + interp_coefficients[2]/3./8. + interp_coefficients[3]/4./16.);
    // deltaz between position[jyp](CELL_YLOW), where t=0.5, and position[jyp](CELL_CENTRE), where t=1
    dimensionless_length_deltas_below(position->jx,position->jyp,position->jz) = coord->dy(position->jx,position->jy)*sqrt(coord->g_22(position->jx,position->jyp))
      *(interp_coefficients[0]/2. + interp_coefficients[1]/2.*3./4. + interp_coefficients[2]/3.*7./8. + interp_coefficients[3]/4.*15./16.);
    next_index_y(position);
    
    do{
      interp_coefficients = cubic_spline_inverse_lambdaC.coefficients(position);
      // deltaz between position[jy](CELL_CENTRE), where t=0, and position[jyp](CELL_YLOW), where t=0.5
      dimensionless_length_deltas_above(position->jx,position->jy,position->jz) = coord->dy(position->jx,position->jy)*sqrt(coord->g_22(position->jx,position->jy))
        *(interp_coefficients[0]/2. + interp_coefficients[1]/2./4. + interp_coefficients[2]/3./8. + interp_coefficients[3]/4./16.);
      // deltaz between position[jyp](CELL_YLOW), where t=0.5, and position[jyp](CELL_CENTRE), where t=1
      dimensionless_length_deltas_below(position->jx,position->jyp,position->jz) = coord->dy(position->jx,position->jyp)*sqrt(coord->g_22(position->jx,position->jyp))
        *(interp_coefficients[0]/2. + interp_coefficients[1]/2.*3./4. + interp_coefficients[2]/3.*7./8. + interp_coefficients[3]/4.*15./16.);
      increasing_dimensionless_length(position->jx,position->jyp,position->jz) = increasing_dimensionless_length[*position] + dimensionless_length_deltas_below[*position] + dimensionless_length_deltas_above[*position];
    } while (next_index_y(position));
    
    {
      Timer timer("comms");
      if (position->jx<mesh->UpXSplitIndex())
	mesh->sendYOutIndest(&increasing_dimensionless_length(position->jx,mesh->yend+1,position->jz),1,
			      NONLOCAL_PARALLEL_TAGBASE + position->jx*mesh->LocalNz+position->jz);
      else
	mesh->sendYOutOutdest(&increasing_dimensionless_length(position->jx,mesh->yend+1,position->jz),1,
			      NONLOCAL_PARALLEL_TAGBASE + position->jx*mesh->LocalNz+position->jz);
    }
    
  } while (next_indexperp(position));
  
  
  // Send the total dimensionless_length at the upper boundary back to the other processors.
  if (has_upper_boundary)
    total_dimensionless_length = sliceXZ(increasing_dimensionless_length, mesh->yend);

  for (int jx=mesh->xstart; jx<=mesh->xend; jx++)
    y_broadcast(&total_dimensionless_length(jx,0), // get a pointer to the jx'th row of total_dimensionless_length's data
		nz, ycomms_size[jx-mesh->xstart]-1,jx);
  
  decreasing_dimensionless_length = -increasing_dimensionless_length;
  
  // Add the FieldPerp total_dimensionless_length to all y locations
  // in the Field3D decreasing_dimensionless_length
  // This loop only from ystart to yend (no Y boundaries)
  for (auto &i : decreasing_dimensionless_length.region(RGN_NOY)) {
    decreasing_dimensionless_length[i] += total_dimensionless_length[i];
  }
  
  if (has_lower_boundary) {
    total_dimensionless_length.setIndex(mesh->ystart-1);
    
    FieldPerp tmp = sliceXZ(dimensionless_length_deltas_below, mesh->ystart-1) + sliceXZ(dimensionless_length_deltas_above, mesh->ystart-1);
    for(auto &i : total_dimensionless_length) {
      decreasing_dimensionless_length[i] += total_dimensionless_length[i];
      decreasing_dimensionless_length[i] += tmp[i];
    }
  }
  
  if (calculate_heatflux) {
    electron_heat_flux = 0.;
  }
  if (calculate_viscosity) {
    electron_viscosity = 0.;
  }
  if (calculate_friction) {
    electron_friction = 0.;
  }
  int driveterm_counter = 0;
  if (gradT_drive) {
    if (calculate_heatflux) {
      electron_heat_flux += -heatflux_zerocoeff * gradT_zerocoeff * interp_to(gradT_driveterm,CELL_YLOW) * gradT_electron; //zero eigenvalue contribution to n^(1,1)/T^1.5
    }
    // viscosity gets no zero eigenvalue contribution
    if (calculate_friction) {
      electron_friction += -friction_zerocoeff * gradT_zerocoeff * interp_to(gradT_driveterm,CELL_YLOW) * gradT_electron; //zero eigenvalue contribution
    }
    cubic_splines_driveterms_centre[driveterm_counter].calculate(gradT_driveterm);
    cubic_splines_driveterms_ylow[driveterm_counter].calculate(gradT_electron);
    driveterm_counter++;
  }
  if (gradV_drive) {
    // gradV drive has no zero eigenvalue contribution
    cubic_splines_driveterms_centre[driveterm_counter].calculate(gradV_driveterm);
    driveterm_counter++;
  }
  if (VeminusVi_drive) {
    if (calculate_heatflux) {
      electron_heat_flux += -heatflux_zerocoeff * VeminusVi_zerocoeff * interp_to(VeminusVi_driveterm,CELL_YLOW)*j_parallel;
    }
    // viscosity gets no zero eigenvalue contribution
    if (calculate_friction) {
      electron_friction += -friction_zerocoeff * VeminusVi_zerocoeff * interp_to(VeminusVi_driveterm,CELL_YLOW)*j_parallel;
    }
    cubic_splines_driveterms_centre[driveterm_counter].calculate(VeminusVi_driveterm);
    cubic_splines_driveterms_ylow[driveterm_counter].calculate(j_parallel);
    driveterm_counter++;
  }
  if (maxwellian_source_drives) {
    calculate_source_driveterms();
    // source_driveterms have no zero eigenvalue contribution
    for (int source_counter=0; source_counter<number_of_source_drives; source_counter++) {
      cubic_splines_driveterms_centre[driveterm_counter].calculate(source_driveterms[source_counter]);
      driveterm_counter++;
    }
  }
  
  {
    BoutReal* temppointer = lower_transients_data;
    for (int i=0; i<transients_data_size; i++) {
      *temppointer=0.;
      temppointer++;
    }
  }
  for (int j=0; j<number_of_negative_eigenvalues; j++) {
    integration.calculateIntegralBelow_cell_ylow(eigenvalues[j], dimensionless_length_deltas_below, dimensionless_length_deltas_above, cubic_spline_inverse_lambdaC, number_of_drives, &driveterm_coefficients_below[j*number_of_drives], cubic_splines_driveterms_centre, cubic_splines_driveterms_ylow, j);
// integration.calculateIntegralAbove_cell_ylow(eigenvalues[j], dimensionless_length_deltas_below, dimensionless_length_deltas_above, cubic_spline_inverse_lambdaC, number_of_drives, &driveterm_coefficients_above[j*number_of_drives], cubic_splines_driveterms_centre, cubic_splines_driveterms_ylow, j);
// output<<endl<<j<<endl;for (int y=0;y<mesh->LocalNy;y++)output<<y<<" "<<integration.integral_below[2][y][0]<<" "<<integration.integral_above[2][y][0]<<endl;if(j==number_of_negative_eigenvalues-1){MPI_Barrier(BoutComm::get());exit(7);}
    if (calculate_heatflux) {
      electron_heat_flux += heatflux_coefficients_below[j]*integration.integral_below;
    }
    if (calculate_viscosity) {
      electron_viscosity += viscosity_coefficients[j]*integration.integral_below;
    }
    if (calculate_friction) {
      electron_friction += friction_coefficients_below[j]*integration.integral_below;
    }
    dataindex = 0;
    if ( yperiodic || (nx_core>0) || (bc_apply_reflecting_sheath && nx_sol>0) )
      for (int jx=mesh->xstart; jx<=mesh->xend; jx++) {
	if (is_upper_boundary[jx]) {
	  if (mesh->periodicY(jx)) {
	    for (int jz=0; jz<nz; jz++) {
	      lower_transients_data[dataindex*transients_row_size+j*nz+jz] += integration.integral_below(jx,mesh->yend+1,jz); // First store total integrals (at first guard cell) in the transients variables
	      dataindex++;
	    }
	  }
	  else if (bc_apply_reflecting_sheath) {
	    for (int jz=0; jz<nz; jz++) {
	      lower_transients_data[dataindex*transients_row_size+j*nz+jz] += integration.integral_below(jx,mesh->yend,jz); // First store total integrals (at last grid cell) in the transients variables
	      dataindex++;
	    }
	  }
	}
      }
  }
  
  // Signs in the following block: the integral is done 'backwards', so start with -'ve sign. [The following signs are now absorbed into the coefficients] Then multiply by (-1)^(l(A)+l(D)), i.e. if drive and output moment have the same l(mod 2) multiply by +1, if different l(mod 2) by -1
  {
    BoutReal* temppointer = upper_transients_data;
    for (int i=0; i<transients_data_size; i++) {
      *temppointer=0.;
      temppointer++;
    }
  }
  for (int j=0; j<number_of_negative_eigenvalues; j++) {
    integration.calculateIntegralAbove_cell_ylow(eigenvalues[j], dimensionless_length_deltas_below, dimensionless_length_deltas_above, cubic_spline_inverse_lambdaC, number_of_drives, &driveterm_coefficients_above[j*number_of_drives], cubic_splines_driveterms_centre, cubic_splines_driveterms_ylow, j);
    
    if (calculate_heatflux) {
      electron_heat_flux += -heatflux_coefficients_above[j]*integration.integral_above;
    }
    if (calculate_viscosity) {
      electron_viscosity += -viscosity_coefficients[j]*integration.integral_above;
    }
    if (calculate_friction) {
      electron_friction += -friction_coefficients_above[j]*integration.integral_above;
    }
    
    dataindex = 0;
    if ( yperiodic || (nx_core>0) || (bc_apply_reflecting_sheath && nx_sol>0) )
      for (int jx=mesh->xstart; jx<=mesh->xend; jx++) {
	if ( is_lower_boundary[jx]) {
	  if (mesh->periodicY(jx)) {
	    for (int jz=0; jz<nz; jz++) {
	      upper_transients_data[dataindex*transients_row_size+j*nz+jz] += -integration.integral_above(jx,mesh->ystart-1,jz); // First store total integrals in the transients variables
	      dataindex++;
	    }
	  }
	  else if (bc_apply_reflecting_sheath) {
	    for (int jz=0; jz<nz; jz++) {
	      upper_transients_data[dataindex*transients_row_size+j*nz+jz] += -integration.integral_above(jx,mesh->ystart,jz); // First store total integrals in the transients variables
	      dataindex++;
	    }
	  }
	}
      }
  }
  
// for (int j=0;j<number_of_negative_eigenvalues;j++)output<<j<<" "<<lower_transients_data[j]<<" "<<upper_transients_data[j]<<endl;MPI_Barrier(BoutComm::get());exit(17);
  if (bc_apply_reflecting_sheath && nx_sol>0 && (has_lower_boundary || has_upper_boundary)) {
    dataindex = 0;
    int solindex = 0;
    bool waitformpi = false;
    MPI_Request request1[nx_sol];
    MPI_Request request2[nx_sol];
    for (int jx=mesh->xstart; jx<=mesh->xend; jx++) {
      if (!mesh->periodicY(jx)) {
	if (is_upper_boundary[jx]) {
	  MPI_Isend(&lower_transients_data[dataindex*transients_row_size], transients_row_size, MPI_DOUBLE, 0, NONLOCAL_PARALLEL_TAGBASE+jx,ycomms[jx-mesh->xstart],&request1[solindex]);
	  MPI_Irecv(&upper_transients_data[dataindex*transients_row_size], transients_row_size, MPI_DOUBLE, 0, NONLOCAL_PARALLEL_TAGBASE+mesh->LocalNx+jx,ycomms[jx-mesh->xstart],&request2[solindex]);
	  waitformpi = true;
	}
	if (is_lower_boundary[jx]) {
	  MPI_Isend(&upper_transients_data[dataindex*transients_row_size], transients_row_size, MPI_DOUBLE, ycomms_size[jx-mesh->xstart]-1, NONLOCAL_PARALLEL_TAGBASE+mesh->LocalNx+jx,ycomms[jx-mesh->xstart],&request1[solindex]);
	  MPI_Irecv(&lower_transients_data[dataindex*transients_row_size], transients_row_size, MPI_DOUBLE, ycomms_size[jx-mesh->xstart]-1, NONLOCAL_PARALLEL_TAGBASE+jx,ycomms[jx-mesh->xstart],&request2[solindex]);
	  waitformpi = true;
	}
	solindex++;
      }
      dataindex++;
    }
    if (waitformpi) {
      MPI_Waitall(nx_sol,request1,MPI_STATUSES_IGNORE);
      MPI_Waitall(nx_sol,request2,MPI_STATUSES_IGNORE);
    }
  }
  
  BoutReal reflecting_sheath_eigenvalue_threshold = 10.;
  for (int j=0; j<number_of_negative_eigenvalues; j++) {
    if (!yperiodic && ( (nx_core>0) || (bc_apply_reflecting_sheath && nx_sol>0) ) ) {
      BoutReal shiftangle;
      FieldPerp decayfactor_core;
      FieldPerp decayfactor_sol;
      BoutReal* upper_transients_data_copy;
      if (has_upper_boundary && nx_core>0) {
        decayfactor_core = exp(sliceXZ(increasing_dimensionless_length, mesh->yend+1)/eigenvalues[j]);
      }
      if (has_lower_boundary && nx_core>0) {
        decayfactor_core = exp(sliceXZ(decreasing_dimensionless_length, mesh->ystart-1)/eigenvalues[j]);
      }
      
      if (has_upper_boundary && nx_sol>0) {
        decayfactor_sol = exp(sliceXZ(increasing_dimensionless_length, mesh->yend)/eigenvalues[j]);
      }
      if (has_lower_boundary && nx_sol>0) {
        decayfactor_sol = exp(sliceXZ(decreasing_dimensionless_length, mesh->ystart)/eigenvalues[j]);
      }
      dataindex=0;
// {MPI_Barrier(BoutComm::get());Timer timer("toroidal");
      for (int jx=mesh->xstart; jx<=mesh->xend; jx++) {
	if (mesh->periodicY(jx,shiftangle)) {
	  toroidal_solver_reverse->solve(lower_transients_data+dataindex*transients_row_size+j*nz,&decayfactor_core(jx,0),jx,shiftangle);
	  dataindex++;
	}
	else if (bc_apply_reflecting_sheath) {
	  if (eigenvalues[j]<reflecting_sheath_eigenvalue_threshold) {
	    // Calculate lower transients on the upper boundary and vice versa in order to have the same communications to do as for the toroidal/periodic transients, for which this choice is more convenient/efficient
	    if (is_lower_boundary[jx]) {
	      if (is_upper_boundary[jx]) {
		// Need to copy lower_transients_data before overwriting it to use in the upper_boundary_transients calculation
		upper_transients_data_copy = new BoutReal[transients_row_size];
		for (int i=0; i<transients_row_size; i++)
		  upper_transients_data_copy[i] = upper_transients_data[(jx-mesh->ystart)*transients_row_size+i];
	      }
	      int i = 0;
	      for (int jz=0; jz<nz; jz++) {
		i = dataindex*transients_row_size+j*nz+jz;
		BoutReal thisdecayfactor = decayfactor_sol(jx, jz);
		upper_transients_data[i] = (upper_transients_data[i]*thisdecayfactor + lower_transients_data[i]) / (1. - pow(thisdecayfactor,2));
	      }
	    }
	    if (is_upper_boundary[jx]) {
	      if (!is_lower_boundary[jx]) {
		upper_transients_data_copy = upper_transients_data+(jx-mesh->ystart)*transients_row_size;
	      }
	      int i = 0;
	      for (int jz=0; jz<nz; jz++) {
		// heatflux and friction need an extra minus sign because we do the calculation with the eigenvector coefficients factored in already and they have odd l values
		i = dataindex*transients_row_size+j*nz+jz;
		BoutReal thisdecayfactor = decayfactor_sol(jx,jz);
		lower_transients_data[i] = (lower_transients_data[i]*thisdecayfactor + upper_transients_data_copy[i]) / (1. - pow(thisdecayfactor,2));
	      }
	      if (is_lower_boundary[jx]) {
		delete [] upper_transients_data_copy;
	      }
	    }
	  }
	}
	dataindex++;
      }
    }
// MPI_Barrier(BoutComm::get());}output<<"toroidal solve time is "<<Timer::getTime("toroidal")<<endl;Timer::resetTime("toroidal");
    else if (yperiodic) {
      if (has_upper_boundary) {
	FieldPerp dividebythis = 1.0 - exp(sliceXZ(increasing_dimensionless_length, mesh->yend+1)/eigenvalues[j]); // Now calculate the actual transient using the stored total integral from zero
	for (int jx=mesh->xstart, i=0; jx<=mesh->xend; jx++, i++)
	  for (int jz=0; jz<nz; jz++) {
	    lower_transients_data[i*transients_row_size+j*nz+jz] /= dividebythis(jx,jz);
	  }
      }
    }
    
    if (!yperiodic) {
      BoutReal shiftangle;
      FieldPerp decayfactor;
      if (has_lower_boundary) decayfactor = exp(sliceXZ(decreasing_dimensionless_length, mesh->ystart-1)/eigenvalues[j]);
      dataindex=0;
// {MPI_Barrier(BoutComm::get());Timer timer("toroidal");
      for (int jx=mesh->xstart; jx<=mesh->xend; jx++)
	if (mesh->periodicY(jx,shiftangle)) {
	  toroidal_solver_forward->solve(upper_transients_data+dataindex*transients_row_size+j*nz,&decayfactor(jx,0),jx,shiftangle);
	  dataindex++;
	}
// MPI_Barrier(BoutComm::get());}output<<"toroidal solve time is "<<Timer::getTime("toroidal")<<endl;Timer::resetTime("toroidal");
    }
    else {
      if (has_lower_boundary) {
	FieldPerp dividebythis = 1.0 - exp(sliceXZ(decreasing_dimensionless_length,mesh->ystart-1)/eigenvalues[j]); // Now calculate the actual transient using the stored total integral from zero
	for (int jx=mesh->xstart, i=0; jx<=mesh->xend; jx++, i++)
	  for (int jz=0; jz<nz; jz++) {
	    upper_transients_data[i*transients_row_size+j*nz+jz] /= dividebythis(jx,jz);
	  }
      }
    }
  }
  
// for (int j=0;j<number_of_negative_eigenvalues;j++)output<<j<<" "<<lower_transients_data[j]<<" "<<upper_transients_data[j]<<endl;MPI_Barrier(BoutComm::get());exit(17);
  if ( yperiodic || (nx_core>0) || (bc_apply_reflecting_sheath && nx_sol>0) ) {
    dataindex = 0;
    for (int jx=mesh->xstart; jx<=mesh->xend; jx++) {
      if (yperiodic || mesh->periodicY(jx) || bc_apply_reflecting_sheath) {
	y_broadcast(&lower_transients_data[dataindex*transients_row_size], transients_row_size, ycomms_size[jx-mesh->xstart]-1, jx);
	y_broadcast(&upper_transients_data[dataindex*transients_row_size], transients_row_size, 0, jx);
	dataindex++;
      }
    }
  }
  
  if (!yperiodic && (nx_sol>0)) {
    for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
      for (int jz=0; jz<nz; jz++) {
	position->jx=rup.ind;
	position->jy=mesh->yend;
	position->jz=jz;
	calc_index(position);
	BoutReal Te_here = interp_to_point_YLOW(T_electron,*position);
	if (bc_heatflux) {
	  pass_interim_upper_boundary_n11(rup.ind,jz) = electron_heat_flux(rup.ind,mesh->yend,jz);
	  if (bc_apply_reflecting_sheath)
	    for (int j=0; j<number_of_negative_eigenvalues; j++)
	      pass_interim_upper_boundary_n11(rup.ind,jz) += heatflux_coefficients_above[j]*upper_transients_data[(rup.ind-mesh->xstart)*transients_row_size+j*nz+jz]
                + heatflux_coefficients_below[j]*lower_transients_data[(rup.ind-mesh->xstart)*transients_row_size+j*nz+jz]*exp(increasing_dimensionless_length(rup.ind,mesh->yend,jz)/eigenvalues[j]);
	  upper_boundary_condition_n11(rup.ind,jz) = -4./5.*sqrt(electron_mass/2.)/pow(Te_here,1.5)*heat_flux_boundary_condition(rup.ind,mesh->yend,jz);
	}
	if (bc_viscosity) {
	  pass_interim_upper_boundary_n20(rup.ind,jz) = electron_viscosity(rup.ind,mesh->yend,jz);
	  if (bc_apply_reflecting_sheath)
	    for (int j=0; j<number_of_negative_eigenvalues; j++)
	      pass_interim_upper_boundary_n20 += viscosity_coefficients[j]*( upper_transients_data[(rup.ind-mesh->xstart)*transients_row_size+j*nz+jz]
                                                                             + lower_transients_data[(rup.ind-mesh->xstart)*transients_row_size+j*nz+jz]*exp(increasing_dimensionless_length(rup.ind,mesh->yend,jz)/eigenvalues[j]) );
	  upper_boundary_condition_n20(rup.ind,jz) = viscosity_boundary_condition(rup.ind,mesh->yend,jz)/Te_here;
	}
      }
    
    if (has_upper_boundary && nx_sol>0 && mesh->getNYPE()>1) {
      MPI_Request request1[nx_sol];
      MPI_Request request2[nx_sol];
      MPI_Request request3[nx_sol];
      MPI_Request request4[nx_sol];
      int i=0;
      for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++) {
	if (rup.ind<mesh->xstart || rup.ind>mesh->xend) continue;
	if (bc_heatflux) {
	  MPI_Isend(&pass_interim_upper_boundary_n11(rup.ind,0), // get a pointer to the jx'th row of FieldPerp's data
		    nz,
		    MPI_DOUBLE,
		    0,
		    NONLOCAL_PARALLEL_TAGBASE,
		    ycomms[rup.ind-mesh->xstart],
		    &request1[i]);
	  MPI_Isend(&upper_boundary_condition_n11(rup.ind,0), // get a pointer to the jx'th row of FieldPerp's data
		    nz,
		    MPI_DOUBLE,
		    0,
		    NONLOCAL_PARALLEL_TAGBASE,
		    ycomms[rup.ind-mesh->xstart],
		    &request2[i]);
	}
	if (bc_viscosity) {
	  MPI_Isend(&pass_interim_upper_boundary_n20(rup.ind,0), // get a pointer to the jx'th row of FieldPerp's data
		    nz,
		    MPI_DOUBLE,
		    0,
		    NONLOCAL_PARALLEL_TAGBASE,
		    ycomms[rup.ind-mesh->xstart],
		    &request3[i]);
	  MPI_Isend(&upper_boundary_condition_n20(rup.ind,0), // get a pointer to the jx'th row of FieldPerp's data
		    nz,
		    MPI_DOUBLE,
		    0,
		    NONLOCAL_PARALLEL_TAGBASE,
		    ycomms[rup.ind-mesh->xstart],
		    &request4[i]);
	}
	i++;
      }
      if (bc_heatflux) {
	MPI_Waitall(nx_sol,request1,MPI_STATUSES_IGNORE);
	MPI_Waitall(nx_sol,request2,MPI_STATUSES_IGNORE);
      }
      if (bc_viscosity) {
	MPI_Waitall(nx_sol,request3,MPI_STATUSES_IGNORE);
	MPI_Waitall(nx_sol,request4,MPI_STATUSES_IGNORE);
      }
    }
    if (has_lower_boundary && nx_sol>0 && mesh->getNYPE()>1) {
      MPI_Request request1[nx_sol];
      MPI_Request request2[nx_sol];
      MPI_Request request3[nx_sol];
      MPI_Request request4[nx_sol];
      int i=0;
      for (RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++) {
	if (rlow.ind<mesh->xstart || rlow.ind>mesh->xend) continue;
	if (bc_heatflux) {
	  MPI_Irecv(&pass_interim_upper_boundary_n11(rlow.ind,0), // get a pointer to the jx'th row of FieldPerp's data
		    nz,
		    MPI_DOUBLE,
		    ycomms_size[rlow.ind-mesh->xstart]-1,
		    NONLOCAL_PARALLEL_TAGBASE,
		    ycomms[rlow.ind-mesh->xstart],
		    &request1[i]);
	  MPI_Irecv(&upper_boundary_condition_n11(rlow.ind,0), // get a pointer to the jx'th row of FieldPerp's data
		    nz,
		    MPI_DOUBLE,
		    ycomms_size[rlow.ind-mesh->xstart]-1,
		    NONLOCAL_PARALLEL_TAGBASE,
		    ycomms[rlow.ind-mesh->xstart],
		    &request2[i]);
	}
	if (bc_viscosity) {
	  MPI_Irecv(&pass_interim_upper_boundary_n20(rlow.ind,0), // get a pointer to the jx'th row of FieldPerp's data
		    nz,
		    MPI_DOUBLE,
		    ycomms_size[rlow.ind-mesh->xstart]-1,
		    NONLOCAL_PARALLEL_TAGBASE,
		    ycomms[rlow.ind-mesh->xstart],
		    &request3[i]);
	  MPI_Irecv(&upper_boundary_condition_n20(rlow.ind,0), // get a pointer to the jx'th row of FieldPerp's data
		    nz,
		    MPI_DOUBLE,
		    ycomms_size[rlow.ind-mesh->xstart]-1,
		    NONLOCAL_PARALLEL_TAGBASE,
		    ycomms[rlow.ind-mesh->xstart],
		    &request4[i]);
	}
	i++;
      }
      if (bc_heatflux) {
	MPI_Waitall(nx_sol,request1,MPI_STATUSES_IGNORE);
	MPI_Waitall(nx_sol,request2,MPI_STATUSES_IGNORE);
      }
      if (bc_viscosity) {
	MPI_Waitall(nx_sol,request3,MPI_STATUSES_IGNORE);
	MPI_Waitall(nx_sol,request4,MPI_STATUSES_IGNORE);
      }
    }
    
    if (has_lower_boundary && nx_sol>0) {
      for (int jx=mesh->xstart, i=0; jx<=mesh->xend; jx++) {
	if (!mesh->periodicY(jx)) {
	  for (int jz=0; jz<nz; jz++) {
	    position->jx=jx;
	    position->jy=mesh->ystart;
	    position->jz=jz;
	    calc_index(position);
	    BoutReal Te_here = interp_to_point_YLOW(T_electron,*position);
	    BoutReal lower_boundary_n11;
	    BoutReal upper_boundary_n11;
	    BoutReal interim_lower_boundary_n11;
	    BoutReal interim_upper_boundary_n11;
	    BoutReal lower_boundary_n20;
	    BoutReal upper_boundary_n20;
	    BoutReal interim_lower_boundary_n20;
	    BoutReal interim_upper_boundary_n20;
	    if (bc_heatflux) {
	      lower_boundary_n11 = -4./5.*sqrt(electron_mass/2.)/pow(Te_here,1.5)*heat_flux_boundary_condition(jx,mesh->ystart,jz);
	      upper_boundary_n11 = upper_boundary_condition_n11(jx,jz);
	      interim_lower_boundary_n11 = electron_heat_flux(jx,mesh->ystart,jz);
	      if (bc_apply_reflecting_sheath)
		for (int j=0; j<number_of_negative_eigenvalues; j++)
		  interim_lower_boundary_n11 += heatflux_coefficients_below[j]*lower_transients_data[(jx-mesh->xstart)*transients_row_size+j*nz+jz]
                    + heatflux_coefficients_above[j]*upper_transients_data[(jx-mesh->xstart)*transients_row_size+j*nz+jz]*exp(decreasing_dimensionless_length(jx,mesh->ystart,jz)/eigenvalues[j]);
	      interim_upper_boundary_n11 = pass_interim_upper_boundary_n11(jx,jz);
	    }
	    if (bc_viscosity) {
	      lower_boundary_n20 = viscosity_boundary_condition(jx,mesh->ystart,jz)/Te_here;
	      upper_boundary_n20 = upper_boundary_condition_n20(jx,jz);
	      interim_lower_boundary_n20 = electron_viscosity(jx,mesh->ystart,jz);
	      if (bc_apply_reflecting_sheath)
		for (int j=0; j<number_of_negative_eigenvalues; j++)
		  interim_lower_boundary_n20 += viscosity_coefficients[j]*( lower_transients_data[(jx-mesh->xstart)*transients_row_size+j*nz+jz]
									    + upper_transients_data[(jx-mesh->xstart)*transients_row_size+j*nz+jz]*exp(decreasing_dimensionless_length(jx,mesh->ystart,jz)/eigenvalues[j]) );
	      interim_upper_boundary_n20 = pass_interim_upper_boundary_n20(jx,jz);
	    }
	    if (bc_heatflux && !bc_viscosity) {
	      /*
	      electron_heat_flux is, at this point, the contribution to n11 from nhat_plus.
	      We want the actual heat flux at mesh->ystart to be boundary_heat_flux.
	      Thus the remainder must come from nhat_minus, which we will construct here just to give the right 1,1 component (could set number_of_negative_eigenvalues-1 more components if desired)
	      However, the transients from the other boundary do not necessarily decay to vanishing by the time they get to this boundary, so we must solve for both.
	      Fortunately this reduces to a single algebraic equation which makes it surprisingly easy to do (in the case of just a single condition being imposed at least).
	      */
	      
	      BoutReal sum_decayed_W11_W11_term = 0.;
	      for (int j=0; j<number_of_negative_eigenvalues; j++) {
		exp_total_dimensionless_length_over_eigenvalue[j] = exp(total_dimensionless_length(jx,jz)/eigenvalues[j]);
		sum_decayed_W11_W11_term += heatflux_coefficients_below[j]*WinverseB_11[j]*exp_total_dimensionless_length_over_eigenvalue[j];
	      }
	      heatflux_sol_transients_factors[i*(2*nz)+jz] = ( (lower_boundary_n11 - interim_lower_boundary_n11)*W11_dot_W11
									- sum_decayed_W11_W11_term*(upper_boundary_n11 - interim_upper_boundary_n11) )
									/ ( pow(W11_dot_W11,2) - pow(sum_decayed_W11_W11_term,2) );
	      heatflux_sol_transients_factors[i*(2*nz)+nz+jz] = ( (upper_boundary_n11 - interim_upper_boundary_n11)*W11_dot_W11
														- sum_decayed_W11_W11_term*(lower_boundary_n11 - interim_lower_boundary_n11) )
														/ ( pow(W11_dot_W11,2) - pow(sum_decayed_W11_W11_term,2) );
	    }
	    if (bc_viscosity && !bc_heatflux) {
	      /*
	      electron_viscosity is, at this point, the contribution to n20 from nhat_plus.
	      We want the actual viscosity at mesh->ystart to be boundary_viscosity.
	      Thus the remainder must come from nhat_minus, which we will construct here just to give the right 2,0 component (could set number_of_negative_eigenvalues-1 more components if desired)
	      However, the transients from the other boundary do not necessarily decay to vanishing by the time they get to this boundary, so we must solve for both.
	      Fortunately this reduces to a single algebraic equation which makes it surprisingly easy to do (in the case of just a single condition being imposed at least).
	      */
	      
	      BoutReal sum_decayed_W20_W20_term = 0.;
	      for (int j=0; j<number_of_negative_eigenvalues; j++) {
		exp_total_dimensionless_length_over_eigenvalue[j] = exp(total_dimensionless_length(jx,jz)/eigenvalues[j]);
		sum_decayed_W20_W20_term += viscosity_coefficients[j]*WinverseB_20[j]*exp_total_dimensionless_length_over_eigenvalue[j];
	      }
	      viscosity_sol_transients_factors[i*(2*nz)+jz] = ( (lower_boundary_n20 - interim_lower_boundary_n20)*W20_dot_W20
									- sum_decayed_W20_W20_term*(upper_boundary_n20 - interim_upper_boundary_n20) )
									/ ( pow(W20_dot_W20,2) - pow(sum_decayed_W20_W20_term,2) );
	      viscosity_sol_transients_factors[i*(2*nz)+nz+jz] = ( (upper_boundary_n20 - interim_upper_boundary_n20)*W20_dot_W20
														- sum_decayed_W20_W20_term*(lower_boundary_n20 - interim_lower_boundary_n20) )
														/ ( pow(W20_dot_W20,2) - pow(sum_decayed_W20_W20_term,2));
	    }
	    if (bc_heatflux && bc_viscosity) {
	      BoutReal sum_decayed_W11_W11_term = 0.;
	      BoutReal sum_decayed_W20_W20_term = 0.;
	      BoutReal sum_decayed_W11_W20_term = 0.;
	      BoutReal sum_decayed_W20_W11_term = 0.;
	      for (int j=0; j<number_of_negative_eigenvalues; j++) {
		exp_total_dimensionless_length_over_eigenvalue[j] = exp(total_dimensionless_length(jx,jz)/eigenvalues[j]);
		sum_decayed_W11_W11_term += heatflux_coefficients_below[j]*WinverseB_11[j]*exp_total_dimensionless_length_over_eigenvalue[j];
		sum_decayed_W20_W20_term += viscosity_coefficients[j]*WinverseB_20[j]*exp_total_dimensionless_length_over_eigenvalue[j];
		sum_decayed_W20_W11_term += viscosity_coefficients[j]*WinverseB_11[j]*exp_total_dimensionless_length_over_eigenvalue[j];
		sum_decayed_W11_W20_term += heatflux_coefficients_below[j]*WinverseB_20[j]*exp_total_dimensionless_length_over_eigenvalue[j];
	      }
	      BoutReal det = pow(sum_decayed_W11_W20_term,2)*pow(sum_decayed_W20_W11_term,2)
			    - 2*sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
			    + pow(sum_decayed_W11_W11_term,2)*pow(sum_decayed_W20_W20_term,2)
			    - pow(sum_decayed_W20_W20_term,2)*pow(W11_dot_W11,2)
			    + 2*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W11*W11_dot_W20
			    - pow(sum_decayed_W20_W11_term,2)*pow(W11_dot_W20,2)
			    - 2*sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W11_dot_W11*W20_dot_W11
			    + 2*sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W11_dot_W20*W20_dot_W11
			    - pow(sum_decayed_W11_W20_term,2)*pow(W20_dot_W11,2)
			    + pow(W11_dot_W20,2)*pow(W20_dot_W11,2)
			    + 2*sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W11_dot_W11*W20_dot_W20
			    - 2*sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W11_dot_W20*W20_dot_W20
			    + 2*sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W11*W20_dot_W20
			    - 2*W11_dot_W11*W11_dot_W20*W20_dot_W11*W20_dot_W20
			    - pow(sum_decayed_W11_W11_term,2)*pow(W20_dot_W20,2)
			    + pow(W11_dot_W11,2)*pow(W20_dot_W20,2);
	      heatflux_sol_transients_factors[i*(2*nz)+jz] = ( (upper_boundary_n11-interim_upper_boundary_n11)*( -sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
																    + sum_decayed_W11_W11_term*pow(sum_decayed_W20_W20_term,2)
																    + sum_decayed_W20_W20_term*W11_dot_W20*W20_dot_W11
																    - sum_decayed_W20_W11_term*W11_dot_W20*W20_dot_W20
																    + sum_decayed_W11_W20_term*W20_dot_W11*W20_dot_W20
																    - sum_decayed_W11_W11_term*pow(W20_dot_W20,2) )
										    + (upper_boundary_n20-interim_upper_boundary_n20)*( pow(sum_decayed_W11_W20_term,2)*sum_decayed_W20_W11_term
																	- sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W20_term
																	+ sum_decayed_W20_W20_term*W11_dot_W11*W11_dot_W20
																	- sum_decayed_W20_W11_term*pow(W11_dot_W20,2)
																	+ sum_decayed_W11_W20_term*W11_dot_W11*W20_dot_W20
																	- sum_decayed_W11_W11_term*W11_dot_W20*W20_dot_W20 )
										    + (lower_boundary_n11-interim_lower_boundary_n11)*( -pow(sum_decayed_W20_W20_term,2)*W11_dot_W11
																	+ sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W20
																	- sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W20_dot_W11
																	+ sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W20_dot_W20
																	- W11_dot_W20*W20_dot_W11*W20_dot_W20 + W11_dot_W11*pow(W20_dot_W20,2) )
										    + (lower_boundary_n20-interim_lower_boundary_n20)*( -sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W11_dot_W11
																	+ sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W11_dot_W20
																	- pow(sum_decayed_W11_W20_term,2)*W20_dot_W11
																	+ pow(W11_dot_W20,2)*W20_dot_W11
																	+ sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W20
																	- W11_dot_W11*W11_dot_W20*W20_dot_W20 )
										) / det;
	      heatflux_sol_transients_factors[i*(2*nz)+nz+jz] = ( (upper_boundary_n11-interim_upper_boundary_n11)*( -pow(sum_decayed_W20_W20_term,2)*W11_dot_W11
																					      + sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W20
																					      - sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W20_dot_W11
																					      + sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W20_dot_W20
																					      - W11_dot_W20*W20_dot_W11*W20_dot_W20
																					      + W11_dot_W11*pow(W20_dot_W20,2) )
															    + (upper_boundary_n20-interim_upper_boundary_n20)*( sum_decayed_W11_W20_term*sum_decayed_W20_W20_term*W11_dot_W11
																						- sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W11_dot_W20
																						+ pow(sum_decayed_W11_W20_term,2)*W20_dot_W11
																						- pow(W11_dot_W20,2)*W20_dot_W11
																						- sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W20
																						+ W11_dot_W11*W11_dot_W20*W20_dot_W20 )
															    + (lower_boundary_n11-interim_lower_boundary_n11)*( -sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
																						+ sum_decayed_W11_W11_term*pow(sum_decayed_W20_W20_term,2)
																						+ sum_decayed_W20_W20_term*W11_dot_W20*W20_dot_W11
																						- sum_decayed_W20_W11_term*W11_dot_W20*W20_dot_W20
																						+ sum_decayed_W11_W20_term*W20_dot_W11*W20_dot_W20
																						- sum_decayed_W11_W11_term*pow(W20_dot_W20,2) )
															    + (lower_boundary_n20-interim_lower_boundary_n20)*( -pow(sum_decayed_W11_W20_term,2)*sum_decayed_W20_W11_term
																						+ sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W20_term
																						- sum_decayed_W20_W20_term*W11_dot_W11*W11_dot_W20
																						+ sum_decayed_W20_W11_term*pow(W11_dot_W20,2)
																						- sum_decayed_W11_W20_term*W11_dot_W11*W20_dot_W20
																						+ sum_decayed_W11_W11_term*W11_dot_W20*W20_dot_W20 )
															  ) / det;
	      viscosity_sol_transients_factors[i*(2*nz)+jz] = ( (upper_boundary_n11-interim_upper_boundary_n11)*( sum_decayed_W11_W20_term*pow(sum_decayed_W20_W11_term,2)
																    - sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
																    - sum_decayed_W20_W20_term*W11_dot_W11*W20_dot_W11
																    - sum_decayed_W11_W20_term*pow(W20_dot_W11,2)
																    + sum_decayed_W20_W11_term*W11_dot_W11*W20_dot_W20
																    + sum_decayed_W11_W11_term*W20_dot_W11*W20_dot_W20 )
										    + (upper_boundary_n20-interim_upper_boundary_n20)*( -sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W11_term
																	+ pow(sum_decayed_W11_W11_term,2)*sum_decayed_W20_W20_term
																	- sum_decayed_W20_W20_term*pow(W11_dot_W11,2)
																	+ sum_decayed_W20_W11_term*W11_dot_W11*W11_dot_W20
																	- sum_decayed_W11_W20_term*W11_dot_W11*W20_dot_W11
																	+ sum_decayed_W11_W11_term*W11_dot_W20*W20_dot_W11 )
										    + (lower_boundary_n11-interim_lower_boundary_n11)*( sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W11
																	- pow(sum_decayed_W20_W11_term,2)*W11_dot_W20
																	+ sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W20_dot_W11
																	+ W11_dot_W20*pow(W20_dot_W11,2)
																	- sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W20_dot_W20
																	- W11_dot_W11*W20_dot_W11*W20_dot_W20 )
										    + (lower_boundary_n20-interim_lower_boundary_n20)*( sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W11_dot_W11
																	- sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W11_dot_W20
																	+ sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W11
																	- W11_dot_W11*W11_dot_W20*W20_dot_W11
																	- pow(sum_decayed_W11_W11_term,2)*W20_dot_W20
																	+ pow(W11_dot_W11,2)*W20_dot_W20 )
										) / det;
	      viscosity_sol_transients_factors[i*(2*nz)+nz+jz] = ( (upper_boundary_n11-interim_upper_boundary_n11)*( -sum_decayed_W20_W11_term*sum_decayed_W20_W20_term*W11_dot_W11
																					      + pow(sum_decayed_W20_W11_term,2)*W11_dot_W20
																					      - sum_decayed_W11_W11_term*sum_decayed_W20_W20_term*W20_dot_W11
																					      - W11_dot_W20*pow(W20_dot_W11,2)
																					      + sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W20_dot_W20
																					      + W11_dot_W11*W20_dot_W11*W20_dot_W20 )
															    + (upper_boundary_n20-interim_upper_boundary_n20)*( sum_decayed_W11_W20_term*sum_decayed_W20_W11_term*W11_dot_W11
																						- sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*W11_dot_W20
																						+ sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*W20_dot_W11
																						- W11_dot_W11*W11_dot_W20*W20_dot_W11
																						- pow(sum_decayed_W11_W11_term,2)*W20_dot_W20
																						+ pow(W11_dot_W11,2)*W20_dot_W20 )
															    + (lower_boundary_n11-interim_lower_boundary_n11)*( -sum_decayed_W11_W20_term*pow(sum_decayed_W20_W11_term,2)
																						+ sum_decayed_W11_W11_term*sum_decayed_W20_W11_term*sum_decayed_W20_W20_term
																						+ sum_decayed_W20_W20_term*W11_dot_W11*W20_dot_W11
																						+ sum_decayed_W11_W20_term*pow(W20_dot_W11,2)
																						- sum_decayed_W20_W11_term*W11_dot_W11*W20_dot_W20
																						- sum_decayed_W11_W11_term*W20_dot_W11*W20_dot_W20 )
															    + (lower_boundary_n20-interim_lower_boundary_n20)*( -sum_decayed_W11_W11_term*sum_decayed_W11_W20_term*sum_decayed_W20_W11_term
																						+ pow(sum_decayed_W11_W11_term,2)*sum_decayed_W20_W20_term
																						- sum_decayed_W20_W20_term*pow(W11_dot_W11,2)
																						+ sum_decayed_W20_W11_term*W11_dot_W11*W11_dot_W20
																						- sum_decayed_W11_W20_term*W11_dot_W11*W20_dot_W11
																						+ sum_decayed_W11_W11_term*W11_dot_W20*W20_dot_W11 )
															  ) / det;
	    }
	  }
	  i++;
	}
      }
    }
  }

  if (bc_heatflux && (nx_sol>0)) {
    for (int jx=mesh->xstart, i=0; jx<=mesh->xend; jx++) {
      if (!mesh->periodicY(jx)) {
	y_broadcast(&heatflux_sol_transients_factors[i*(2*nz)], (nz)*2, 0, jx);
	i++;
      }
    }
  }
  if (bc_viscosity && (nx_sol>0)) {
    for (int jx=mesh->xstart, i=0; jx<=mesh->xend; jx++) {
      if (!mesh->periodicY(jx)) {
	y_broadcast(&viscosity_sol_transients_factors[i*(2*nz)], (nz)*2, 0, jx);
	i++;
      }
    }
  }
  
  if (nx_sol>0) {
    for (int jx=mesh->xstart, i=0; jx<=mesh->xend; jx++, i++) {
      if (!mesh->periodicY(jx)) {
	for (int jz=0; jz<nz; jz++) {
	  if (bc_heatflux) {
	    for (int j=0; j<number_of_negative_eigenvalues; j++) {
	      lower_transients_data[i*transients_row_size+j*nz+jz] += WinverseB_11[j]*heatflux_sol_transients_factors[i*(2*nz)+jz];
	      upper_transients_data[i*transients_row_size+j*nz+jz] += -WinverseB_11[j]*heatflux_sol_transients_factors[i*(2*nz)+nz+jz];
	    }
	  }
	  if (bc_viscosity) {
	    for (int j=0; j<number_of_negative_eigenvalues; j++) {
	      lower_transients_data[i*transients_row_size+j*nz+jz] += WinverseB_20[j]*viscosity_sol_transients_factors[i*(2*nz)+jz];
	      upper_transients_data[i*transients_row_size+j*nz+jz] += WinverseB_20[j]*viscosity_sol_transients_factors[i*(2*nz)+nz+jz];
	    }
	  }
	}
      }
    }
  }

  #ifndef NOEDGETERMS
  start_index(position);
  do {
    position->jy=mesh->ystart;
    calc_index(position);
    do {
      for (int j=0; j<number_of_negative_eigenvalues; j++) {
	BoutReal lower_transient = lower_transients_data[(position->jx-mesh->ystart)*transients_row_size+j*nz+position->jz] * exp(increasing_dimensionless_length[*position]/eigenvalues[j]);
	BoutReal upper_transient = upper_transients_data[(position->jx-mesh->ystart)*transients_row_size+j*nz+position->jz] * exp(decreasing_dimensionless_length[*position]/eigenvalues[j]);
	if (calculate_heatflux) {
	  electron_heat_flux[*position] += heatflux_coefficients_below[j]*lower_transient + heatflux_coefficients_above[j]*upper_transient;
	}
	if (calculate_viscosity) {
	  electron_viscosity[*position] += viscosity_coefficients[j]*(lower_transient + upper_transient);
	}
	if (calculate_friction) {
	  electron_friction[*position] += friction_coefficients_below[j]*lower_transient + friction_coefficients_above[j]*upper_transient;
	}
      }
      position->jy++;
      calc_index(position);
    } while (position->jy<mesh->yend+1);
  } while (next_indexperp(position));
  #endif
  
  if (calculate_heatflux) {
    electron_heat_flux *= -5./4.*sqrt(2./electron_mass)*interp_to(pow(T_electron,1.5),CELL_YLOW); //now we have q=-5/4*v_Telectron*T_electron*n^(1,1)
    mesh->communicate(electron_heat_flux);
  }
  if (calculate_viscosity) {
    electron_viscosity *= interp_to(T_electron,CELL_YLOW);
    mesh->communicate(electron_viscosity);
  }
  if (calculate_friction) {
    electron_friction *= 2.*interp_to(T_electron*lambdaC_inverse,CELL_YLOW);
    electron_friction += -interp_to(sqrt(2.*electron_mass*T_electron)*lambdaC_inverse,CELL_YLOW)*(-j_parallel/elementary_charge);
    mesh->communicate(electron_friction);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NonLocalParallel::calculate_source_driveterms() { // NB must be called AFTER calculating lambdaC_inverse
  Field3D one_minus_tau = 1.-source_tau;
  source_driveterms[0] = source_strength*one_minus_tau*one_minus_tau/lambdaC_inverse/sqrt(2*T_electron/electron_mass);
  for (int i=1; i<number_of_source_drives; i++) {
    source_driveterms[i] = source_driveterms[i-1]*one_minus_tau;
  }
}

void NonLocalParallel::set_boundary_gradients() {
  for (RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++)
    for (int jz=0; jz<nz; jz++) {
      BoutReal heat_flux_boundarygradient;
      BoutReal viscosity_boundarygradient;
      BoutReal friction_boundarygradient;
      if (calculate_heatflux) {
//       heat_flux_boundarygradient = (electron_heat_flux[rlow.ind][mesh->ystart][jz]-27.*electron_heat_flux[rlow.ind][mesh->ystart+1][jz]+27.*electron_heat_flux[rlow.ind][mesh->ystart+2][jz]-electron_heat_flux[rlow.ind][mesh->ystart+3][jz])/24.; // NB gradient in index space
//       heat_flux_boundarygradient = (-11.*electron_heat_flux[rlow.ind][mesh->ystart][jz] + 18.*electron_heat_flux[rlow.ind][mesh->ystart+1][jz] - 9.*electron_heat_flux[rlow.ind][mesh->ystart+2][jz] + 2.*electron_heat_flux[rlow.ind][mesh->ystart+3][jz]) / 6. / coord->dy[rlow.ind][mesh->ystart] / sqrt((coord->g_22[rlow.ind][mesh->ystart]+coord->g_22[rlow.ind][mesh->ystart+1]+coord->g_22[rlow.ind][mesh->ystart+2]+coord->g_22[rlow.ind][mesh->ystart+3])/4.);
	heat_flux_boundarygradient = (-electron_heat_flux(rlow.ind,mesh->ystart,jz) + electron_heat_flux(rlow.ind,mesh->ystart+boundary_gradient_smoothing_length,jz))/BoutReal(boundary_gradient_smoothing_length); // NB gradient in index space
      }
      if (calculate_viscosity) {
//       viscosity_boundarygradient = (electron_viscosity[rlow.ind][mesh->ystart][jz]-27.*electron_viscosity[rlow.ind][mesh->ystart+1][jz]+27.*electron_viscosity[rlow.ind][mesh->ystart+2][jz]-electron_viscosity[rlow.ind][mesh->ystart+3][jz])/24.; // NB gradient in index space
//       viscosity_boundarygradient = (-11.*electron_viscosity[rlow.ind][mesh->ystart][jz] + 18.*electron_viscosity[rlow.ind][mesh->ystart+1][jz] - 9.*electron_viscosity[rlow.ind][mesh->ystart+2][jz] + 2.*electron_viscosity[rlow.ind][mesh->ystart+3][jz]) / 6. / coord->dy[rlow.ind][mesh->ystart] / sqrt((coord->g_22[rlow.ind][mesh->ystart]+coord->g_22[rlow.ind][mesh->ystart+1]+coord->g_22[rlow.ind][mesh->ystart+2]+coord->g_22[rlow.ind][mesh->ystart+3])/4.);
	viscosity_boundarygradient = (-electron_viscosity(rlow.ind,mesh->ystart,jz) + electron_viscosity(rlow.ind,mesh->ystart+boundary_gradient_smoothing_length,jz))/BoutReal(boundary_gradient_smoothing_length); // NB gradient in index space
      }
      if (calculate_friction) {
//       friction_boundarygradient = (electron_friction[rlow.ind][mesh->ystart][jz]-27.*electron_friction[rlow.ind][mesh->ystart+1][jz]+27.*electron_friction[rlow.ind][mesh->ystart+2][jz]-electron_friction[rlow.ind][mesh->ystart+3][jz])/24.; // NB gradient in index space
//       friction_boundarygradient = (-11.*electron_friction[rlow.ind][mesh->ystart][jz] + 18.*electron_friction[rlow.ind][mesh->ystart+1][jz] - 9.*electron_friction[rlow.ind][mesh->ystart+2][jz] + 2.*electron_friction[rlow.ind][mesh->ystart+3][jz]) / 6. / coord->dy[rlow.ind][mesh->ystart] / sqrt((coord->g_22[rlow.ind][mesh->ystart]+coord->g_22[rlow.ind][mesh->ystart+1]+coord->g_22[rlow.ind][mesh->ystart+2]+coord->g_22[rlow.ind][mesh->ystart+3])/4.);
	friction_boundarygradient = (-electron_friction(rlow.ind,mesh->ystart,jz) + electron_friction(rlow.ind,mesh->ystart+boundary_gradient_smoothing_length,jz))/BoutReal(boundary_gradient_smoothing_length); // NB gradient in index space
      }
      for (int jy=mesh->ystart-1; jy>=0; jy--) {
	if (calculate_heatflux) {
	  electron_heat_flux(rlow.ind,jy,jz) = electron_heat_flux(rlow.ind,jy+1,jz) - heat_flux_boundarygradient;
	}
	if (calculate_viscosity) {
	  electron_viscosity(rlow.ind,jy,jz) = electron_viscosity(rlow.ind,jy+1,jz) - viscosity_boundarygradient;
	}
	if (calculate_friction) {
	  electron_friction(rlow.ind,jy,jz) = electron_friction(rlow.ind,jy+1,jz) - friction_boundarygradient;
	}
      }
    }
  for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
    for (int jz=0; jz<nz; jz++) {
      BoutReal heat_flux_boundarygradient;
      BoutReal viscosity_boundarygradient;
      BoutReal friction_boundarygradient;
      if (calculate_heatflux) {
// 	heat_flux_boundarygradient = (electron_heat_flux(rup.ind,mesh->yend-3,jz)-27.*electron_heat_flux[rup.ind][mesh->yend-2][jz]+27.*electron_heat_flux[rup.ind][mesh->yend-1][jz]-electron_heat_flux[rup.ind][mesh->yend][jz])/24.; // NB gradient in index space
// 	heat_flux_boundarygradient = (11.*electron_heat_flux[rup.ind][mesh->yend][jz] - 18.*electron_heat_flux[rup.ind][mesh->yend-1][jz] + 9.*electron_heat_flux[rup.ind][mesh->yend-2][jz] - 2.*electron_heat_flux[rup.ind][mesh->yend-3][jz]) / 6. / coord->dy[rup.ind][mesh->yend] / sqrt((coord->g_22[rup.ind][mesh->yend]+coord->g_22[rup.ind][mesh->yend-1]+coord->g_22[rup.ind][mesh->yend-2]+coord->g_22[rup.ind][mesh->yend-3])/4.);
	heat_flux_boundarygradient = (-electron_heat_flux(rup.ind,mesh->yend-boundary_gradient_smoothing_length,jz) + electron_heat_flux(rup.ind,mesh->yend,jz))/BoutReal(boundary_gradient_smoothing_length); // NB gradient in index space
      }
      if (calculate_viscosity) {
// 	viscosity_boundarygradient = (electron_viscosity[rup.ind][mesh->yend-3][jz]-27.*electron_viscosity[rup.ind][mesh->yend-2][jz]+27.*electron_viscosity[rup.ind][mesh->yend-1][jz]-electron_viscosity[rup.ind][mesh->yend][jz])/24.; // NB gradient in index space
//      viscosity_boundarygradient = (11.*electron_viscosity[rup.ind][mesh->yend][jz] - 18.*electron_viscosity[rup.ind][mesh->yend-1][jz] + 9.*electron_viscosity[rup.ind][mesh->yend-2][jz] - 2.*electron_viscosity[rup.ind][mesh->yend-3][jz]) / 6. / coord->dy[rup.ind][mesh->yend] / sqrt((coord->g_22[rup.ind][mesh->yend]+coord->g_22[rup.ind][mesh->yend-1]+coord->g_22[rup.ind][mesh->yend-2]+coord->g_22[rup.ind][mesh->yend-3])/4.);
	viscosity_boundarygradient = (-electron_viscosity(rup.ind,mesh->yend-boundary_gradient_smoothing_length,jz) + electron_viscosity(rup.ind,mesh->yend,jz))/BoutReal(boundary_gradient_smoothing_length); // NB gradient in index space
      }
      if (calculate_friction) {
// 	friction_boundarygradient = (electron_friction[rup.ind][mesh->yend-3][jz]-27.*electron_friction[rup.ind][mesh->yend-2][jz]+27.*electron_friction[rup.ind][mesh->yend-1][jz]-electron_friction[rup.ind][mesh->yend][jz])/24.; // NB gradient in index space
//      friction_boundarygradient = (11.*electron_friction[rup.ind][mesh->yend][jz] - 18.*electron_friction[rup.ind][mesh->yend-1][jz] + 9.*electron_friction[rup.ind][mesh->yend-2][jz] - 2.*electron_friction[rup.ind][mesh->yend-3][jz]) / 6. / coord->dy[rup.ind][mesh->yend] / sqrt((coord->g_22[rup.ind][mesh->yend]+coord->g_22[rup.ind][mesh->yend-1]+coord->g_22[rup.ind][mesh->yend-2]+coord->g_22[rup.ind][mesh->yend-3])/4.);
	friction_boundarygradient = (-electron_friction(rup.ind,mesh->yend-boundary_gradient_smoothing_length,jz) + electron_friction(rup.ind,mesh->yend,jz))/BoutReal(boundary_gradient_smoothing_length); // NB gradient in index space
      }
      for (int jy=mesh->yend+1; jy<mesh->LocalNy; jy++) {
	if (calculate_heatflux) {
	  electron_heat_flux(rup.ind,jy,jz) = electron_heat_flux(rup.ind,jy-1,jz) + heat_flux_boundarygradient;
	}
	if (calculate_viscosity) {
	  electron_viscosity(rup.ind,jy,jz) = electron_viscosity(rup.ind,jy-1,jz) + viscosity_boundarygradient;
	}
	if (calculate_friction) {
	  electron_friction(rup.ind,jy,jz) = electron_friction(rup.ind,jy-1,jz) + friction_boundarygradient;
	}
      }
    }
}

void NonLocalParallel::set_neumann_boundary_conditions() {
  for (RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++)
    for (int jy=0; jy<mesh->ystart; jy++)
      for (int jz=0; jz<nz; jz++) {
	if (calculate_heatflux) {
	  electron_heat_flux(rlow.ind,jy,jz)= electron_heat_flux(rlow.ind,mesh->ystart,jz);
	}
	if (calculate_viscosity) {
	  electron_viscosity(rlow.ind,jy,jz)= electron_viscosity(rlow.ind,mesh->ystart,jz);
	}
	if (calculate_friction) {
	  electron_friction(rlow.ind,jy,jz)= electron_friction(rlow.ind,mesh->ystart,jz);
	}
      }
  for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
    for (int jy=mesh->yend+1; jy<mesh->LocalNy; jy++)
      for (int jz=0; jz<nz; jz++) {
	if (calculate_heatflux) {
	  electron_heat_flux(rup.ind,jy,jz)= electron_heat_flux(rup.ind,mesh->yend,jz);
	}
	if (calculate_viscosity) {
	  electron_viscosity(rup.ind,jy,jz)= electron_viscosity(rup.ind,mesh->yend,jz);
	}
	if (calculate_friction) {
	  electron_friction(rup.ind,jy,jz)= electron_friction(rup.ind,mesh->yend,jz);
	}
      }
}

void NonLocalParallel::y_broadcast(BoutReal* input_buffer, const int &size, const int &root_processor, const int &jx) {
  // NB Assumes that the mesh is BoutMesh
  Timer timer("comms");
  
//  MPI_Bcast(input_buffer, size, PVEC_REAL_MPI_TYPE, root_processor, comm_yprocs);
  MPI_Bcast(input_buffer, size, MPI_DOUBLE, root_processor, ycomms[jx-mesh->xstart]);
// Should use commented out version if method is transferred to boutmesh.cxx
}

// void NonLocalParallel::y_boundary_broadcast(BoutReal* input_buffer, const int &size, const int &root_processor) {
//   // NB Assumes that the mesh is BoutMesh
//   Timer timer("comms");
//   #ifdef CHECK
//     int root_processor_yindex = (root_processor - mesh->getXProcIndex())/mesh->getNXPE();
//     if (root_processor_yindex!=0)
//       throw BoutException("y_boundary_broadcast: it has been assumed that root_processor_yindex==0, but this is not the case here");
//     if (root_processor_yindex==mesh->getNYPE()-1)
//       throw BoutException("y_boundary_broadcast: sends to root_processor+1 but root_processor is the last one");
//   #endif
//   int processor = mesh->getYProcIndex() * mesh->getNXPE() + mesh->getXProcIndex();
//   
//   if (mesh->getNYPE()==1)
//     return;
//   
//   if (mesh->getYProcIndex()==0) {
//     MPI_Wait(&broadcast_request, MPI_STATUS_IGNORE);
//     broadcast_request = mesh->sendToProc(mesh->getXProcIndex(),1,input_buffer,size,NONLOCAL_PARALLEL_TAGBASE);
//   }
//   else {
//     if (mesh->getYProcIndex()==1) {
//       mesh->wait(mesh->receiveFromProc(mesh->getXProcIndex(),0,input_buffer,size,NONLOCAL_PARALLEL_TAGBASE));
//     }
//     
//   //  MPI_Bcast(input_buffer, size, PVEC_REAL_MPI_TYPE, root_processor, comm_yprocs);
//   if (mesh->getNYPE()>2)
//     MPI_Bcast(input_buffer, size, MPI_DOUBLE, root_processor+1, comm_yprocs_minusone);
//   // Should use commented out version if method is transferred to boutmesh.cxx
//   }
// }

// void NonLocalParallel::rms_over_y(const Field3D &input_field, FieldPerp &output_field, int jx) {
//   FieldPerp tempsum;
//   tempsum = 0.;
//   int ye = mesh->yend;
//   if (mesh->StaggerGrids && input_field.getLocation()==CELL_CENTRE) ye--;
//   for (int jx=mesh->xstart; jx<=mesh->xend; jx++)
//     for (int jz=0; jz<nz; jz++)
//       for (int jy=mesh->ystart; jy<=ye; jy++) {
// 	tempsum[jx][jz]+=pow(input_field[jx][jy][jz],2);
//       }
//   MPI_Reduce(*tempsum.getData(),
// 	     *output_field.getData(),
// 	     mesh->LocalNx*mesh->LocalNz,
// 	     MPI_DOUBLE,
// 	     MPI_SUM,
// 	     mesh->getXProcIndex(),
// 	     comm_yprocs);
//   if (mesh->getYProcIndex()==0) {
//     int ny = mesh->GlobalNy;
//     if (mesh->StaggerGrids && input_field.getLocation()==CELL_CENTRE) ny--;
//     for (int jx=mesh->xstart; jx<=mesh->xend; jx++)
//       for (int jz=0; jz<nz;jz++)
// 	output_field[jx][jz] = sqrt(output_field[jx][jz]/ny);
//     mesh->sendToProc(mesh->getXProcIndex(),mesh->getNYPE()-1,*output_field.getData(),mesh->LocalNx*mesh->LocalNz,NONLOCAL_PARALLEL_TAGBASE);
//   }
//   else if (mesh->getYProcIndex()==mesh->getNYPE()-1) {
//     mesh->wait(mesh->receiveFromProc(mesh->getXProcIndex(),0,*output_field.getData(),mesh->LocalNx*mesh->LocalNz,NONLOCAL_PARALLEL_TAGBASE));
//   }
// }

// void NonLocalParallel::mean_over_y(const Field3D &input_field, FieldPerp &output_field, int exclude_edgecells) {
//   FieldPerp tempsum;
//   tempsum = 0.;
//   int ys = mesh->ystart+exclude_edgecells;
//   int ye = mesh->yend-exclude_edgecells;
//   
//   if (mesh->StaggerGrids && input_field.getLocation()==CELL_CENTRE && mesh->lastY()) ye--;
//   
//   for (int jx=mesh->xstart; jx<=mesh->xend; jx++)
//     for (int jz=0; jz<nz; jz++)
//       for (int jy=ys; jy<=ye; jy++) {
// 	tempsum[jx][jz]+=input_field[jx][jy][jz];
//       }
//   MPI_Reduce(*tempsum.getData(),
// 	     *output_field.getData(),
// 	     mesh->LocalNx*mesh->LocalNz,
// 	     MPI_DOUBLE,
// 	     MPI_SUM,
// 	     mesh->getXProcIndex(),
// 	     comm_yprocs);
//   if (mesh->getYProcIndex()==0) {
//     int ny = mesh->GlobalNy;
//     if (mesh->StaggerGrids && input_field.getLocation()==CELL_CENTRE) ny--;
//     ny-=2*exclude_edgecells;
//     for (int jx=mesh->xstart; jx<=mesh->xend; jx++)
//       for (int jz=0; jz<nz;jz++)
// 	output_field[jx][jz] = output_field[jx][jz]/ny;
//     MPI_Request request = mesh->sendToProc(mesh->getXProcIndex(),mesh->getNYPE()-1,*output_field.getData(),mesh->LocalNx*mesh->LocalNz,NONLOCAL_PARALLEL_TAGBASE);
//     MPI_Waitall(1,&request,MPI_STATUSES_IGNORE);
//   }
//   else if (mesh->getYProcIndex()==mesh->getNYPE()-1) {
//     mesh->wait(mesh->receiveFromProc(mesh->getXProcIndex(),0,*output_field.getData(),mesh->LocalNx*mesh->LocalNz,NONLOCAL_PARALLEL_TAGBASE));
//   }
// }

BoutReal NonLocalParallel::interp_to_point_YLOW(const Field3D &input, bindex &position) {
   if(mesh->StaggerGrids)
     return (9.*(input(position.jx,position.jym,position.jz)+input(position.jx,position.jy,position.jz))-(input(position.jx,position.jy2m,position.jz)+input(position.jx,position.jyp,position.jz)))/16.;
   else
     return input[position];
}
