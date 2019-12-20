/*
   SD1D: 1D simulation of plasma-neutral interactions
   ==================================================

     Copyright B.Dudson, University of York, 2016
              email: benjamin.dudson@york.ac.uk
              
    This file is part of SD1D.

    SD1D is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SD1D is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SD1D.  If not, see <http://www.gnu.org/licenses/>.


  Normalisations
  --------------
 
  Ne   (density) normalised to Nnorm [m^-3]
  T    (temperature) normalised to Tnorm [eV]
  B    (magnetic field) normalised to Bnorm [eV]
  
  t    (time) normalised using ion cyclotron frequency Omega_ci [1/s]
  Vi   (velocity) normalised to sound speed Cs [m/s]
  L    (lengths) normalised to hybrid Larmor radius rho_s = Cs/Omega_ci [m]
  
 */

#include <mpi.h>

#include <bout/physicsmodel.hxx>
#include <bout/constants.hxx>
#include <derivs.hxx>
#include <field_factory.hxx>
#include <invert_parderiv.hxx>

#include "div_ops.hxx"
#include "loadmetric.hxx"
#include "radiation.hxx"

#include "non-local_parallel.hxx"  // Non-local heat transport

class SD1D : public PhysicsModel {
protected:
  int init(bool restarting) {
    Options *opt = Options::getRoot()->getSection("sd1d");

    OPTION(opt, cfl_info, false); // Calculate and print CFL information
    
    // Normalisation
    OPTION(opt, Tnorm, 100);  // Reference temperature [eV]
    OPTION(opt, Nnorm, 1e19); // Reference density [m^-3]
    OPTION(opt, Bnorm, 1.0);  // Reference magnetic field [T]
    OPTION(opt, AA, 2.0);     // Ion mass
    SAVE_ONCE4(Tnorm, Nnorm, Bnorm, AA); // Save normalisations
    
    // Model parameters
    OPTION(opt, vwall, 1./3);   // 1/3rd Franck-Condon energy at wall
    OPTION(opt, frecycle, 1.0); // Recycling fraction 100%
    OPTION(opt, fredistribute, 0.0); // Fraction of neutrals redistributed evenly along leg
    OPTION(opt, density_sheath, 0); // Free boundary
    OPTION(opt, pressure_sheath, 0); // Free boundary
    OPTION(opt, gaspuff,  0.0); // Additional gas flux at target
    OPTION(opt, dneut, 1.0);    // Scale neutral gas diffusion
    OPTION(opt, nloss, 0.0);    // Neutral gas loss rate
    OPTION(opt, fimp,  0.01);   // 1% impurity
    OPTION(opt, Eionize,   30);     // Energy loss per ionisation (30eV)
    OPTION(opt, symmetry_plane, true); // Sheath heat transmission
    OPTION(opt, sheath_gamma, 6.5); // Sheath heat transmission
    OPTION(opt, neutral_gamma, 0.25); // Neutral heat transmission
    
    // Plasma anomalous transport
    OPTION(opt, anomalous_D, -1);
    OPTION(opt, anomalous_chi, -1);
    
    if(sheath_gamma < 6) 
      throw BoutException("sheath_gamma < 6 not consistent");

    OPTION(opt, tn_floor, 3.5);  // Minimum neutral gas temperature [eV]

    OPTION(opt, atomic, true);

    OPTION(opt, neutral_f_pn, true);

    OPTION(opt, hyper, -1);            // Numerical hyper-diffusion
    OPTION(opt, ADpar, -1);            // Added Dissipation scheme
    OPTION(opt, viscos, -1);           // Parallel viscosity
    OPTION(opt, ion_viscosity, false); // Braginskii parallel viscosity
    OPTION(opt, heat_conduction, true); // Spitzer-Hahm heat conduction

    OPTION(opt, nonlocal_conduction, false); // Non-local heat conduction
    
    OPTION(opt, elastic_scattering, false); // Include ion-neutral elastic scattering?
    OPTION(opt, excitation, false);    // Include electron impact excitation?

    OPTION(opt, gamma_sound, 5./3);  // Ratio of specific heats
    OPTION(opt, bndry_flux_fix, false);
   
    OPTION(opt, density_form, 4);
    OPTION(opt, momentum_form, 6);
    OPTION(opt, energy_form, 8);
    
    // Read the flux-tube area from input file
    // This goes into the Jacobian.
    string area_string;
    FieldFactory ffact(mesh); 
    
    // Calculate normalisation factors
    
    Cs0      = sqrt(SI::qe*Tnorm / (AA*SI::Mp)); // Reference sound speed [m/s]
    Omega_ci = SI::qe*Bnorm / (AA*SI::Mp);       // Ion cyclotron frequency [1/s]
    rho_s0   = Cs0 / Omega_ci;                   // Length scale [m]
  
    mi_me  = AA*SI::Mp/SI::Me;

    BoutReal Coulomb = 6.6 - 0.5*log(Nnorm * 1e-20) + 1.5*log(Tnorm);
    tau_e0 = 1. / (2.91e-6 * (Nnorm / 1e6) * Coulomb * pow(Tnorm, -3./2));
    
    // Save normalisation factors
    SAVE_ONCE5(Cs0, Omega_ci, rho_s0, tau_e0, mi_me);
    
    OPTION(opt, volume_source, true);
    if(volume_source) {
      // Volume sources of particles and energy
      
      string source_string;
    
      Options* optne = Options::getRoot()->getSection("Ne");
      optne->get("source", source_string, "0.0");
      NeSource = ffact.create2D(source_string, optne);
      //SAVE_ONCE(NeSource);
      
      Options* optpe = Options::getRoot()->getSection("P");
      optpe->get("source", source_string, "0.0");
      PeSource = ffact.create2D(source_string, optpe);
      SAVE_ONCE(PeSource);
      
      // Normalise sources
      NeSource /= Nnorm * Omega_ci;
      PeSource /= SI::qe * Nnorm * Tnorm * Omega_ci;
    }else {
      // Point sources, fixing density and specifying energy flux
      
      
      Options* optpe = Options::getRoot()->getSection("P");
      OPTION(optpe, powerflux, 2e7); // Power flux in W/m^2
      powerflux /= rho_s0 * SI::qe * Tnorm * Nnorm * Omega_ci; // Normalised energy flux
    }

    /////////////////////////
    // Density controller
    OPTION(opt, density_upstream, -1); // Fix upstream density? [m^-3]
    if(density_upstream > 0.0) {
      // Fixing density
      density_upstream /= Nnorm;
      
      // Controller
      OPTION(opt, density_controller_p, 1e-2);
      OPTION(opt, density_controller_i, 1e-3);
      OPTION(opt, density_integral_positive, false);
      OPTION(opt, density_source_positive, true);

      density_error_lasttime = -1.0;  // Signal no value

      if (not symmetry_plane) {
        // Need to control density at the mid-point, so find if the mid-point or
        // either/both mid-points are on this processor

        if (not volume_source) {
          throw BoutException("Using symmetry_plane=false requires volume_source=true");
        }

        // Note GlobalNy includes boundary cells, but YLOCAL takes a global index that
        // does not include boundary cells
        if (mesh->GlobalNy%2 == 0) {
          // even number of grid points, so two symmetric mid-points
          //mid_point1 = mesh->YLOCAL(mesh->GlobalNy/2 - mesh->ystart - 1); // requires BOUT++ v4.2
          mid_point1 = (mesh->GlobalNy/2 - mesh->ystart - 1) % (mesh->LocalNy - 2*mesh->ystart);
          if (mid_point1 < mesh->ystart or mid_point1 > mesh->yend) {
            mid_point1 = -1;
          }
          mid_point_proc1 = (mesh->GlobalNy/2 - mesh->ystart - 1)
                            /(mesh->LocalNy - 2*mesh->ystart);
          //mid_point2 = mesh->YLOCAL(mesh->GlobalNy/2 - mesh->ystart); // requires BOUT++ v4.2
          mid_point2 = (mesh->GlobalNy/2 - mesh->ystart) % (mesh->LocalNy - 2*mesh->ystart);
          if (mid_point2 < mesh->ystart or mid_point2 > mesh->yend) {
            mid_point2 = -1;
          }
          mid_point_proc2 = (mesh->GlobalNy/2 - mesh->ystart)
                            /(mesh->LocalNy - 2*mesh->ystart);
        } else {
          // odd number of grid points, so single mid-point
          //mid_point1 = mesh->YLOCAL(mesh->GlobalNy/2 - mesh->ystart); // requires BOUT++ v4.2
          mid_point1 = (mesh->GlobalNy/2 - mesh->ystart) % (mesh->LocalNy - 2*mesh->ystart);
          if (mid_point1 < mesh->ystart or mid_point1 > mesh->yend) {
            mid_point1 = -1;
          }
          mid_point_proc1 = (mesh->GlobalNy/2 - mesh->ystart)
                            /(mesh->LocalNy - 2*mesh->ystart);
        }
      }

      // Save and load error integral from file, since
      // this determines the source function
      solver->addToRestart(density_error_integral, "density_error_integral");
      solver->addToRestart(density_error_integral2, "density_error_integral2");
      
      if(!restarting) {
        density_error_integral = 0.0;
        density_error_integral2 = 0.0;
        
        if(volume_source) {
          // Set density_error_integral so that
          // the input source is used
          density_error_integral = 1./density_controller_i;
          density_error_integral2 = 1./density_controller_i;
        }
      }
    }

    if(volume_source) {
      if(density_upstream > 0.0) {
        // Evolving NeSource
        SAVE_REPEAT(NeSource);

        NeSource0 = NeSource; // Save initial value
      }else {
        // Fixed NeSource
        SAVE_ONCE(NeSource);
      }
    }
    
    Options::getRoot()->getSection("NVn")->get("evolve", evolve_nvn, true);
    Options::getRoot()->getSection("Pn")->get("evolve", evolve_pn, true);
    
    nloss /= Omega_ci;

    // Specify variables to evolve
    solver->add(Ne, "Ne");
    solver->add(NVi, "NVi");
    solver->add(P, "P");

    if(atomic) {
      solver->add(Nn, "Nn");
      if(evolve_nvn) {
        solver->add(NVn, "NVn");
      }
      if(evolve_pn) {
        solver->add(Pn, "Pn");
      }
    }
    
    // Load the metric tensor
    LoadMetric(rho_s0, Bnorm);
    
    opt->get("area", area_string, "1.0");
    mesh->coordinates()->J = ffact.create2D(area_string, Options::getRoot());
    
    dy4 = SQ(SQ(mesh->coordinates()->dy));
    
    // Use carbon radiation for the impurity
    rad = new HutchinsonCarbonRadiation();

    // Add extra quantities to be saved
    if(atomic) {
      SAVE_REPEAT4(S, R, E, F); // Save net plasma particle source, radiated power, energy transfer, friction
      SAVE_REPEAT2(Dn, kappa_n);  // Neutral diffusion coefficients
      SAVE_REPEAT(flux_ion);   // Flux of ions to target
    }
    if(heat_conduction)
      SAVE_REPEAT(kappa_epar); // Save coefficient of thermal conduction
    
    bool diagnose;
    OPTION(opt, diagnose, true);
    if(diagnose) {
      // Output extra variables
      if(atomic) {
        SAVE_REPEAT2(Srec,Siz);       // Save particle sources
        SAVE_REPEAT3(Frec,Fiz,Fcx);   // Save momentum sources
        SAVE_REPEAT3(Rrec,Riz,Rzrad); // Save radiation sources
        SAVE_REPEAT3(Erec,Eiz,Ecx);   // Save energy transfer

        if(elastic_scattering) {
          SAVE_REPEAT2(Fel, Eel); // Elastic collision transfer channels
        }
        if(excitation) {
          SAVE_REPEAT(Rex); // Electron-neutral excitation
        }
        
        if(evolve_nvn) {
          SAVE_REPEAT(Vn);
        }
      }
      
      SAVE_REPEAT(Vi);
    }
    
    if(ion_viscosity) 
      SAVE_REPEAT(eta_i);
    
    kappa_epar = 0.0;

    Srec = 0.0; Siz = 0.0; S = 0.0;
    Frec = 0.0; Fiz = 0.0; Fcx = 0.0; Fel = 0.0;  F = 0.0;
    Rrec = 0.0; Riz = 0.0; Rzrad = 0.0; Rex = 0.0; R = 0.0;
    Erec = 0.0; Eiz = 0.0; Ecx = 0.0; Eel = 0.0;  E = 0.0;
    
    flux_ion = 0.0;
    
    // Neutral gas diffusion and heat conduction
    Dn = 0.0;
    kappa_n = 0.0;

    // Anomalous transport
    if(anomalous_D > 0.0) {
      // Normalise
      anomalous_D /= rho_s0*rho_s0*Omega_ci; // m^2/s
      output.write("\tnormalised anomalous D_perp = %e\n", anomalous_D);
    }
    if(anomalous_chi > 0.0) {
      // Normalise
      anomalous_chi /= rho_s0*rho_s0*Omega_ci; // m^2/s
      output.write("\tnormalised anomalous chi_perp = %e\n", anomalous_chi);
    }
    
    // Calculate neutral gas redistribution weights over the domain
    string redist_string;
    opt->get("redist_weight", redist_string, "1.0");
    redist_weight = ffact.create2D(redist_string, opt);
    BoutReal localweight = 0.0;
    for(int j=mesh->ystart;j<=mesh->yend;j++) {
      localweight += redist_weight(mesh->xstart, j);
    }

    MPI_Comm ycomm = mesh->getYcomm(mesh->xstart); // MPI communicator

    // Calculate total weight by summing over all processors
    BoutReal totalweight;
    MPI_Allreduce(&localweight, &totalweight, 1, MPI_DOUBLE, MPI_SUM, ycomm);
    // Normalise redist_weight so sum over domain is 1
    redist_weight /= totalweight;
    
    setPrecon( (preconfunc) &SD1D::precon );

    //////////////////////////////////////////
    // Split operator (IMEX) schemes
    // Use combination of explicit and implicit methods
    //
    // Boolean flags rhs_explicit and rhs_implicit
    // turn on explicit and implicit terms 

    bool split_operator;
    OPTION(opt, split_operator, false);
    if(!split_operator) {
      // Turn on all terms in rhs
      rhs_explicit = rhs_implicit = true;
      update_coefficients = true;
    }
    setSplitOperator(split_operator);
    
    //////////////////////////////////////////
    // Non-local heat conduction
    
    if ( heat_conduction && nonlocal_conduction ) {
      // Initialise the non-local heat conduction model
      nonlocal_parallel = new NonLocalParallel(-SI::qe, SI::Me, SI::e0, 
                                               Coulomb, // Coulomb logarithm
                                               false,   // fluxes ylow
                                               Options::getRoot()->getSection("non_local_parallel"));

      // No currents or potentials for now
      nonlocal_parallel->set_potential(0.0);
      nonlocal_parallel->set_j_parallel(0.0);
      
      SAVE_REPEAT2(qe_nonlocal, qe_local);
    }

    return 0;
  }

  /*!
   * This function calculates the time derivatives
   * of all evolving quantities
   *
   */
  int rhs(BoutReal time) {
    fprintf(stderr, "\rTime: %e", time);

    Coordinates *coord = mesh->coordinates();

    mesh->communicate(Ne, NVi, P);
    
    // Floor small values
    P = floor(P, 1e-10);
    Ne = floor(Ne, 1e-10);
    
    Field3D Nelim = floor(Ne, 1e-5);
    
    Vi = NVi / Ne;
    
    Field3D Te = 0.5*P / Ne; // Assuming Te = Ti
    
    mesh->communicate(Vi, Te);
    
    Field3D Nnlim;
    Field3D Tn;
    if(atomic) {
      // Includes atomic processes, neutral gas
      mesh->communicate(Nn);
      if(evolve_nvn) {
        mesh->communicate(NVn);
      }
      if(evolve_pn) {
        mesh->communicate(Pn);
      }
      Nn = floor(Nn, 1e-10);
      Nnlim = floor(Nn, 1e-5);

      
      if(evolve_nvn) {
        Vn = NVn / Nnlim;
      }else {
        Vn = - vwall * sqrt(3.5/Tnorm);
        NVn = Nn * Vn;
      }
      
      if(evolve_pn) {
        Tn = Pn / Nnlim;
        //Tn = floor(Tn, 0.025/Tnorm); // Minimum tn_floor
        Tn = floor(Tn, 1e-12); 
      }else {
        Tn = Te; // Strong CX coupling
        Pn = Tn * floor(Nn, 0.0);
        Tn = floor(Tn, tn_floor/Tnorm); // Minimum of tn_floor
      }
    }
    
    if(update_coefficients) {
      // Update diffusion coefficients
      TRACE("Update coefficients");
      
      tau_e = Omega_ci * tau_e0 * pow(Te,1.5)/Ne;
      
      if(heat_conduction) {
        kappa_epar = 3.2 * mi_me * 0.5*P * tau_e;
        kappa_epar.applyBoundary("neumann");
      }
    
      if(atomic) {
        // Neutral diffusion rate
        
        for(int i=0;i<mesh->LocalNx;i++)
          for(int j=0;j<mesh->LocalNy;j++)
            for(int k=0;k<mesh->LocalNz;k++) {
              // Charge exchange frequency, normalised to ion cyclotron frequency
              BoutReal sigma_cx = Nelim(i,j,k)*Nnorm*hydrogen.chargeExchange(Te(i,j,k)*Tnorm)/Omega_ci;
              
              // Ionisation frequency
              BoutReal sigma_iz = Nelim(i,j,k)*Nnorm*hydrogen.ionisation(Te(i,j,k)*Tnorm)/Omega_ci;
              
              // Neutral thermal velocity
              BoutReal tn = Tn(i,j,k);
              if(tn < tn_floor/Tnorm)
                tn = tn_floor/Tnorm;
              BoutReal vth_n = sqrt(tn); // Normalised to Cs0
              
              // Neutral-neutral mean free path
              BoutReal Lmax = 1.0; // meters
              BoutReal a0 = PI*SQ(5.29e-11);
              BoutReal lambda_nn = 1. / (Nnorm*Nnlim(i,j,k)*a0); // meters
              if(lambda_nn > Lmax) {
                // Limit maximum mean free path
                lambda_nn = Lmax;
              }
              
              lambda_nn /= rho_s0; // Normalised length to rho_s0
              // Neutral-Neutral collision rate
              BoutReal sigma_nn = vth_n / lambda_nn;
              
              // Total neutral collision frequency
              BoutReal sigma = sigma_cx + sigma_iz + sigma_nn;
              
              
              // Neutral gas diffusion
              if(dneut > 0.0) {
                Dn(i,j,k) = dneut * SQ(vth_n) / sigma;
              }
              
              // Neutral gas heat conduction
              kappa_n(i,j,k) = dneut * Nnlim(i,j,k) * SQ(vth_n) / sigma;
            }
        
        kappa_n.applyBoundary("Neumann");
        Dn.applyBoundary("dirichlet_o2");
        mesh->communicate(kappa_n, Dn);
      }
    }
    
    // Set sheath boundary condition on flow
    
    TRACE("Sheath");
    ddt(P) = 0.0; // Need to set heat flux

    if(evolve_pn) {
      ddt(Pn) = 0.0;
    }
    
    for(RangeIterator r=mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      int jz = 0;
      
      // Outward flow velocity to >= Cs
      BoutReal Vout = sqrt(2.0*Te(r.ind, mesh->yend, jz)); // Sound speed outwards
      if(Vi(r.ind, mesh->yend, jz) > Vout)
        Vout = Vi(r.ind, mesh->yend, jz); // If plasma is faster, go to plasma velocity
      
      BoutReal Nout;
      switch( density_sheath ) {
      case 0: {
        // Free boundary on density (constant gradient)
        Nout = 0.5*( 3.*Ne(r.ind, mesh->yend, jz) - Ne(r.ind, mesh->yend-1, jz) );
        break;
      }
      case 1: {
        // Zero gradient
        Nout = Ne(r.ind, mesh->yend, jz); 
        break;
      }
      case 2: {
        // Zero gradient particle flux N*Vi* J*dx*dz
        // Since Vi increases into the sheath, density should drop
        Nout = Ne(r.ind, mesh->yend, jz) * coord->J(r.ind, mesh->yend) * Vi(r.ind, mesh->yend, jz) / (0.5*(coord->J(r.ind, mesh->yend) + coord->J(r.ind, mesh->yend+1)) * Vout);
        break;
      }
      default: throw BoutException("Unrecognised density_sheath option");
      }
      
      if(Nout < 0.0)
        Nout = 0.0; // Prevent Nout from going negative -> Flux is always to the wall
      
      // Flux of particles is Ne*Vout
      BoutReal flux = Nout * Vout;
      
      BoutReal Pout;
      
      switch( pressure_sheath ) {
      case 0: {
        // Free boundary  (constant gradient)
        Pout = 0.5*( 3.*P(r.ind, mesh->yend, jz) - P(r.ind, mesh->yend-1, jz) );
        break;
      }
      case 1: {
        // Zero gradient
        Pout = P(r.ind, mesh->yend, jz); 
        break;
      }
      case 2: {
        // Use energy flux conservation to set pressure
        // (5/2)Pv + (1/2)nv^3 = const
        // 
        Pout = ((5.*P(r.ind,mesh->yend,jz)*Vi(r.ind,mesh->yend,jz) + Ne(r.ind,mesh->yend,jz)*pow(Vi(r.ind,mesh->yend,jz),3))/Vout - Nout*Vout*Vout)/5.;
        break;
      }
      default: throw BoutException("Unrecognised pressure_sheath option");
      }

      if(Pout < 0.0)
        Pout = 0.0;
      
      if(rhs_explicit && (!nonlocal_conduction)) {
        // Additional cooling
        // If nonlocal heat conduction is used then this is
        // already included.
        BoutReal q = (sheath_gamma - 6) * Te(r.ind, mesh->yend, jz) * flux;
        
        // Multiply by cell area to get power
        BoutReal heatflux = q * (coord->J(r.ind, mesh->yend)+coord->J(r.ind, mesh->yend+1))/(sqrt(coord->g_22(r.ind, mesh->yend)) + sqrt(coord->g_22(r.ind, mesh->yend+1)));
        
        // Divide by volume of cell, and 2/3 to get pressure
        ddt(P)(r.ind, mesh->yend, jz) -= (2./3)*heatflux / (coord->dy(r.ind, mesh->yend)*coord->J(r.ind, mesh->yend));
      }
      
      // Set boundary half-way between cells
      for(int jy=mesh->yend+1; jy<mesh->LocalNy; jy++) {
        
        ///// Plasma model
        
        // Vi fixed value (Dirichlet)
        Vi(r.ind, jy, jz)  = 2.*Vout - Vi(r.ind, mesh->yend, jz);
        
        // Ne set from flux (Dirichlet)
        Ne(r.ind, jy, jz)  = 2*Nout - Ne(r.ind, mesh->yend, jz);
        
        // NVi. This can be negative, so set this to the flux
        // going out of the domain (zero gradient)
        NVi(r.ind, jy, jz) = Nout * Vout;
        //NVi(r.ind, jy, jz) = Ne(r.ind, jy, jz)  * Vi(r.ind, jy, jz);        
        //NVi(r.ind, jy, jz) = 2.*Nout * Vout - NVi(r.ind, mesh->yend, jz);

        
        if (nonlocal_conduction) {
          // Constant gradient
          Te(r.ind, jy, jz) = 2.*Te(r.ind, mesh->yend, jz) - Te(r.ind, mesh->yend-1, jz);
        } else {
          // Te zero gradient (Neumann)
          Te(r.ind, jy, jz) = Te(r.ind, mesh->yend, jz);
        }
        
        P(r.ind, jy, jz) = 2.*Pout - P(r.ind, mesh->yend, jz);
        
        if(atomic) {
          ///// Neutral model
          // Flux of neutral particles, momentum, and energy are set later
          // Here the neutral velocity is set to no-flow conditions
          
          // Vn fixed value (Dirichlet)
          Vn(r.ind, jy, jz)  = - Vn(r.ind, mesh->yend, jz);
          
          // Nn free boundary (constant gradient)
          Nn(r.ind, jy, jz) = 2.*Nn(r.ind, mesh->yend, jz) - Nn(r.ind, mesh->yend-1, jz);
          
          if(evolve_pn) {
            // Tn fixed value (Dirichlet)
            //Tn(r.ind, jy, jz) = 3.5/Tnorm - Tn(r.ind, mesh->yend, jz);
            
            // Tn zero gradient. Heat flux set by gamma
            Tn(r.ind, jy, jz) = Tn(r.ind, mesh->yend, jz);
            
            if(rhs_explicit && (neutral_gamma > 0.0)) {
              // Density at the target
              BoutReal Nnout = 0.5*(Nn(r.ind, mesh->yend, jz) + Nn(r.ind, mesh->yend+1, jz));
              // gamma * n * T * cs
              BoutReal q = neutral_gamma * Nnout * Tn(r.ind, jy, jz) * sqrt(Tn(r.ind, jy, jz));
              
              // Multiply by cell area to get power
              BoutReal heatflux = q * (coord->J(r.ind, mesh->yend)+coord->J(r.ind, mesh->yend+1))/(sqrt(coord->g_22(r.ind, mesh->yend)) + sqrt(coord->g_22(r.ind, mesh->yend+1)));
              
              // Divide by volume of cell, and 2/3 to get pressure
              ddt(Pn)(r.ind, mesh->yend, jz) -= (2./3)*heatflux / (coord->dy(r.ind, mesh->yend)*coord->J(r.ind, mesh->yend));
            }
          }else {
            Tn(r.ind, jy, jz) = Te(r.ind, jy, jz);
          }
          Pn(r.ind, jy, jz) = Nn(r.ind, jy, jz) * Tn(r.ind, jy, jz); 
          NVn(r.ind, jy, jz) = -NVn(r.ind, mesh->yend, jz);
        }
      }
    }
    
    if (symmetry_plane) {
      for(RangeIterator r=mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        // No-flow boundary condition on left boundary
        
        for(int jz=0; jz<mesh->LocalNz; jz++) {
          for(int jy=0; jy<mesh->ystart; jy++) {
            Te(r.ind, jy, jz) = Te(r.ind, mesh->ystart, jz);
            Ne(r.ind, jy, jz) = Ne(r.ind, mesh->ystart, jz);
            P(r.ind, jy, jz) = P(r.ind, mesh->ystart, jz);
            Vi(r.ind, jy, jz) = -Vi(r.ind, mesh->ystart, jz);
            NVi(r.ind, jy, jz) = -NVi(r.ind, mesh->ystart, jz);
            
            if(atomic) {
              Vn(r.ind, jy, jz) = -Vn(r.ind, mesh->ystart, jz);
              Nn(r.ind, jy, jz) = Nn(r.ind, jy, jz);
              Pn(r.ind, jy, jz) = Pn(r.ind, jy, jz);
              Tn(r.ind, jy, jz) = Tn(r.ind, jy, jz);
            }
          }
        }
      }
    } else {
      for(RangeIterator r=mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        int jz = 0;
        
        // Outward flow velocity to >= Cs
        BoutReal Vout = -sqrt(2.0*Te(r.ind, mesh->ystart, jz)); // Sound speed outwards
        if(Vi(r.ind, mesh->ystart, jz) > Vout)
          Vout = Vi(r.ind, mesh->ystart, jz); // If plasma is faster, go to plasma velocity
        
        BoutReal Nout;
        switch( density_sheath ) {
        case 0: {
          // Free boundary on density (constant gradient)
          Nout = 0.5*( 3.*Ne(r.ind, mesh->ystart, jz) - Ne(r.ind, mesh->ystart+1, jz) );
          break;
        }
        case 1: {
          // Zero gradient
          Nout = Ne(r.ind, mesh->ystart, jz); 
          break;
        }
        case 2: {
          // Zero gradient particle flux N*Vi* J*dx*dz
          // Since Vi increases into the sheath, density should drop
          Nout = Ne(r.ind, mesh->ystart, jz) * coord->J(r.ind, mesh->ystart) * Vi(r.ind, mesh->ystart, jz) / (0.5*(coord->J(r.ind, mesh->ystart) + coord->J(r.ind, mesh->ystart-1)) * Vout);
          break;
        }
        default: throw BoutException("Unrecognised density_sheath option");
        }
        
        if(Nout < 0.0)
          Nout = 0.0; // Prevent Nout from going negative -> Flux is always to the wall
        
        // Flux of particles is Ne*Vout
        BoutReal flux = Nout * Vout;
        
        BoutReal Pout;
        
        switch( pressure_sheath ) {
        case 0: {
          // Free boundary  (constant gradient)
          Pout = 0.5*( 3.*P(r.ind, mesh->ystart, jz) - P(r.ind, mesh->ystart+1, jz) );
          break;
        }
        case 1: {
          // Zero gradient
          Pout = P(r.ind, mesh->ystart, jz); 
          break;
        }
        case 2: {
          // Use energy flux conservation to set pressure
          // (5/2)Pv + (1/2)nv^3 = const
          // 
          Pout = ((5.*P(r.ind,mesh->ystart,jz)*Vi(r.ind,mesh->ystart,jz) + Ne(r.ind,mesh->ystart,jz)*pow(Vi(r.ind,mesh->ystart,jz),3))/Vout - Nout*Vout*Vout)/5.;
          break;
        }
        default: throw BoutException("Unrecognised pressure_sheath option");
        }

        if(Pout < 0.0)
          Pout = 0.0;
        
        if(rhs_explicit && (!nonlocal_conduction)) {
          // Additional cooling
          // If nonlocal heat conduction is used then this is
          // already included.
          BoutReal q = (sheath_gamma - 6) * Te(r.ind, mesh->ystart, jz) * flux;
          
          // Multiply by cell area to get power
          BoutReal heatflux = q * (coord->J(r.ind, mesh->ystart)+coord->J(r.ind, mesh->ystart-1))/(sqrt(coord->g_22(r.ind, mesh->ystart)) + sqrt(coord->g_22(r.ind, mesh->ystart-1)));
          
          // Divide by volume of cell, and 2/3 to get pressure
          ddt(P)(r.ind, mesh->ystart, jz) -= (2./3)*heatflux / (coord->dy(r.ind, mesh->ystart)*coord->J(r.ind, mesh->ystart));
        }
        
        // Set boundary half-way between cells
        for(int jy=mesh->ystart-1; jy>0; jy--) {
          
          ///// Plasma model
          
          // Vi fixed value (Dirichlet)
          Vi(r.ind, jy, jz)  = 2.*Vout - Vi(r.ind, mesh->ystart, jz);
          
          // Ne set from flux (Dirichlet)
          Ne(r.ind, jy, jz)  = 2*Nout - Ne(r.ind, mesh->ystart, jz);
          
          // NVi. This can be positive, so set this to the flux
          // going out of the domain (zero gradient)
          NVi(r.ind, jy, jz) = Nout * Vout;
          //NVi(r.ind, jy, jz) = Ne(r.ind, jy, jz)  * Vi(r.ind, jy, jz);        
          //NVi(r.ind, jy, jz) = 2.*Nout * Vout - NVi(r.ind, mesh->ystart, jz);

          
          if (nonlocal_conduction) {
            // Constant gradient
            Te(r.ind, jy, jz) = 2.*Te(r.ind, mesh->ystart, jz) - Te(r.ind, mesh->ystart+1, jz);
          } else {
            // Te zero gradient (Neumann)
            Te(r.ind, jy, jz) = Te(r.ind, mesh->ystart, jz);
          }
          
          P(r.ind, jy, jz) = 2.*Pout - P(r.ind, mesh->ystart, jz);
          
          if(atomic) {
            ///// Neutral model
            // Flux of neutral particles, momentum, and energy are set later
            // Here the neutral velocity is set to no-flow conditions
            
            // Vn fixed value (Dirichlet)
            Vn(r.ind, jy, jz)  = - Vn(r.ind, mesh->ystart, jz);
            
            // Nn free boundary (constant gradient)
            Nn(r.ind, jy, jz) = 2.*Nn(r.ind, mesh->ystart, jz) - Nn(r.ind, mesh->ystart+1, jz);
            
            if(evolve_pn) {
              // Tn fixed value (Dirichlet)
              //Tn(r.ind, jy, jz) = 3.5/Tnorm - Tn(r.ind, mesh->ystart, jz);
              
              // Tn zero gradient. Heat flux set by gamma
              Tn(r.ind, jy, jz) = Tn(r.ind, mesh->ystart, jz);
              
              if(rhs_explicit && (neutral_gamma > 0.0)) {
                // Density at the target
                BoutReal Nnout = 0.5*(Nn(r.ind, mesh->ystart, jz) + Nn(r.ind, mesh->ystart-1, jz));
                // gamma * n * T * cs
                BoutReal q = neutral_gamma * Nnout * Tn(r.ind, jy, jz) * sqrt(Tn(r.ind, jy, jz));
                
                // Multiply by cell area to get power
                BoutReal heatflux = q * (coord->J(r.ind, mesh->ystart)+coord->J(r.ind, mesh->ystart-1))/(sqrt(coord->g_22(r.ind, mesh->ystart)) + sqrt(coord->g_22(r.ind, mesh->ystart-1)));
                
                // Divide by volume of cell, and 2/3 to get pressure
                ddt(Pn)(r.ind, mesh->ystart, jz) -= (2./3)*heatflux / (coord->dy(r.ind, mesh->ystart)*coord->J(r.ind, mesh->ystart));
              }
            }else {
              Tn(r.ind, jy, jz) = Te(r.ind, jy, jz);
            }
            Pn(r.ind, jy, jz) = Nn(r.ind, jy, jz) * Tn(r.ind, jy, jz); 
            NVn(r.ind, jy, jz) = -NVn(r.ind, mesh->ystart, jz);
          }
        }
      }
    }
    
    if((density_upstream > 0.0) && rhs_explicit) {
      ///////////////////////////////////////////////
      // Set midplane density
      //
      // This calculates a source needed at the midplane grid cell(s), to relax towards
      // the desired density value. 
      //

      TRACE("Density upstream");
      
      BoutReal source=0., source2=0.;
      if (symmetry_plane) {
        for(RangeIterator r=mesh->iterateBndryLowerY(); !r.isDone(); r++) {
          int jz = 0;
            
          // Density source, so dn/dt = source
          BoutReal error = density_upstream - Ne(r.ind, mesh->ystart, jz);
          
          // PI controller, using crude integral of the error
          if(density_error_lasttime < 0.0) {
            // First time
            density_error_lasttime = time;
            density_error_last = error;
          }
          
          // Integrate using Trapezium rule
          if(time > density_error_lasttime) { // Since time can decrease
            density_error_integral += (time - density_error_lasttime)*
              0.5*(error + density_error_last);
          }
          
          if((density_error_integral < 0.0) && density_integral_positive) {
            // Limit density_error_integral to be >= 0
            density_error_integral = 0.0;
          }
          
          // Calculate source from combination of error and integral
          source = density_controller_p * error + density_controller_i * density_error_integral;
          
          //output.write("\n Source: %e, %e : %e, %e -> %e\n", time, (time - density_error_lasttime), error, density_error_integral, source);

          density_error_last = error;
          density_error_lasttime = time;

          if(!volume_source) {
            // Convert source into a flow velocity
            // through the boundary, based on a zero-gradient boundary on the density.
            // This ensures that the mass and momentum inputs are consistent,
            // but also carries energy through the boundary. This flux
            // of energy is calculated, and subtracted from the pressure equation,
            // so that the density boundary does not contribute to energy balance.
            
            // Calculate needed input velocity
            BoutReal Vin = source * sqrt(coord->g_22(r.ind, mesh->ystart))*coord->dy(r.ind, mesh->ystart) / Ne(r.ind, mesh->ystart, jz);
            
            // Limit at sound speed
            BoutReal cs = sqrt(Te(r.ind, mesh->ystart, jz));
            if( fabs(Vin) > cs ) {
              Vin *= cs / fabs(Vin); // + or - cs
            }
            Vi(r.ind, mesh->ystart-1, jz) = 2.*Vin - Vi(r.ind, mesh->ystart, jz);
            
            // Power flux is v * (5/2 P + 1/2 m n v^2 )
            BoutReal inputflux = Vin * ( 2.5 * P(r.ind, mesh->ystart, jz) + 0.5*Ne(r.ind, mesh->ystart, jz) * Vin*Vin );  // W/m^2 (normalised)
            
            // Subtract input energy flux from P equation
            // so no net power input
            ddt(P)(r.ind, mesh->ystart, jz) -= (2./3)*inputflux / ( coord->dy(r.ind, mesh->ystart) * sqrt(coord->g_22(r.ind, mesh->ystart)));
          }
        }
      } else {
        BoutReal error, error2;
        if (mid_point1 > -1) {
          int jz = 0;
            
          // Density source, so dn/dt = source
          error = density_upstream - Ne(mesh->xstart, mid_point1, jz);
          
          // PI controller, using crude integral of the error
          if(density_error_lasttime < 0.0) {
            // First time
            density_error_lasttime = time;
            density_error_last = error;
          }
          
          // Integrate using Trapezium rule
          if(time > density_error_lasttime) { // Since time can decrease
            density_error_integral += (time - density_error_lasttime)*
              0.5*(error + density_error_last);
          }
          
          if((density_error_integral < 0.0) && density_integral_positive) {
            // Limit density_error_integral to be >= 0
            density_error_integral = 0.0;
          }
          
          // Calculate source from combination of error and integral
          source = density_controller_p * error + density_controller_i * density_error_integral;
          
          //output.write("\n Source: %e, %e : %e, %e -> %e\n", time, (time - density_error_lasttime), error, density_error_integral, source);

          density_error_last = error;
          density_error_lasttime = time;
        }
        if (mid_point2 > -1) {
          int jz = 0;
            
          // Density source, so dn/dt = source
          error2 = density_upstream - Ne(mesh->xstart, mid_point2, jz);
          
          // PI controller, using crude integral of the error
          if(density_error_lasttime2 < 0.0) {
            // First time
            density_error_lasttime2 = time;
            density_error_last2 = error2;
          }
          
          // Integrate using Trapezium rule
          if(time > density_error_lasttime2) { // Since time can decrease
            density_error_integral2 += (time - density_error_lasttime2)*
              0.5*(error2 + density_error_last2);
          }
          
          if((density_error_integral2 < 0.0) && density_integral_positive) {
            // Limit density_error_integral to be >= 0
            density_error_integral2 = 0.0;
          }
          
          // Calculate source from combination of error and integral
          source2 = density_controller_p * error + density_controller_i * density_error_integral;
          
          //output.write("\n Source: %e, %e : %e, %e -> %e\n", time, (time - density_error_lasttime2), error2, density_error_integral2, source);

          density_error_last2 = error2;
          density_error_lasttime2 = time;
        }
      }

      if(volume_source) {
        if((source < 0.0) && density_source_positive) {
          source = 0.0; // Don't remove particles
        }
        
        if (symmetry_plane) {
          // Broadcast the value of source from processor 0
          MPI_Bcast(&source, 1, MPI_DOUBLE, 0, BoutComm::get());
        } else {
          // Broadcast the value of source from midplane processor(s)
          MPI_Bcast(&source, 1, MPI_DOUBLE, mid_point_proc1, BoutComm::get());
          if (mid_point_proc2 > -1) {
            if((source2 < 0.0) && density_source_positive) {
              source2 = 0.0; // Don't remove particles
            }
            MPI_Bcast(&source2, 1, MPI_DOUBLE, mid_point_proc2, BoutComm::get());

            // If there are two mid-plane points, average the sources
            source = (source + source2)/2.;
          }
        }
        
        // Scale NeSource
        NeSource = source * NeSource0;
      }
    }
    
    if(atomic && rhs_explicit) {
      // Atomic physics
      TRACE("Atomic");
      
      if(fimp > 0.0) {
        // Impurity radiation 
        Rzrad = rad->power(Te*Tnorm, Ne*Nnorm, Ne*(Nnorm*fimp)); // J / m^3 / s
        Rzrad /= SI::qe*Tnorm*Nnorm * Omega_ci; // Normalise
      } // else Rzrad = 0.0 set in init()
      
      E = 0.0; // Energy transfer to neutrals

      // Lower floor on Nn for atomic rates
      Field3D Nnlim2 = floor(Nn, 0.0);
      
      for(int i=0;i<mesh->LocalNx;i++)
        for(int j=mesh->ystart;j<=mesh->yend;j++)
          for(int k=0;k<mesh->LocalNz;k++) {
            
            // Integrate rates over each cell using Simpson's rule
            // Calculate cell centre (C), left (L) and right (R) values
            
            BoutReal Te_C = Te(i,j,k), Te_L = 0.5*(Te(i,j-1,k) + Te(i,j,k)), Te_R = 0.5*(Te(i,j,k) + Te(i,j+1,k));
            BoutReal Ne_C = Ne(i,j,k), Ne_L = 0.5*(Ne(i,j-1,k) + Ne(i,j,k)), Ne_R = 0.5*(Ne(i,j,k) + Ne(i,j+1,k));
            BoutReal Vi_C = Vi(i,j,k), Vi_L = 0.5*(Vi(i,j-1,k) + Vi(i,j,k)), Vi_R = 0.5*(Vi(i,j,k) + Vi(i,j+1,k));
            BoutReal Tn_C = Tn(i,j,k), Tn_L = 0.5*(Tn(i,j-1,k) + Tn(i,j,k)), Tn_R = 0.5*(Tn(i,j,k) + Tn(i,j+1,k));
            BoutReal Nn_C = Nnlim2(i,j,k), Nn_L = 0.5*(Nnlim2(i,j-1,k) + Nnlim2(i,j,k)), Nn_R = 0.5*(Nnlim2(i,j,k) + Nnlim2(i,j+1,k));
            BoutReal Vn_C = Vn(i,j,k), Vn_L = 0.5*(Vn(i,j-1,k) + Vn(i,j,k)), Vn_R = 0.5*(Vn(i,j,k) + Vn(i,j+1,k));

            // Jacobian (Cross-sectional area)
            BoutReal J_C = coord->J(i,j), J_L = 0.5*(coord->J(i,j-1) + coord->J(i,j)), J_R = 0.5*(coord->J(i,j) + coord->J(i,j+1));
            
            ///////////////////////////////////////
            // Charge exchange
            
            BoutReal R_cx_L = Ne_L*Nn_L*hydrogen.chargeExchange(Te_L*Tnorm) * (Nnorm / Omega_ci);
            BoutReal R_cx_C = Ne_C*Nn_C*hydrogen.chargeExchange(Te_C*Tnorm) * (Nnorm / Omega_ci);
            BoutReal R_cx_R = Ne_R*Nn_R*hydrogen.chargeExchange(Te_R*Tnorm) * (Nnorm / Omega_ci);
            
            // Ecx is energy transferred to neutrals
            Ecx(i,j,k) = (3./2)* (
                                       J_L * (Te_L - Tn_L)*R_cx_L
                                  + 4.*J_C * (Te_C - Tn_C)*R_cx_C
                                  +    J_R * (Te_R - Tn_R)*R_cx_R
                                  ) / (6. * J_C);
            
            // Fcx is friction between plasma and neutrals 
            Fcx(i,j,k) = (
                               J_L * (Vi_L - Vn_L)*R_cx_L
                          + 4.*J_C * (Vi_C - Vn_C)*R_cx_C
                          +    J_R * (Vi_R - Vn_R)*R_cx_R
                          ) / (6. * J_C);
            
            ///////////////////////////////////////
            // Recombination
            
            BoutReal R_rc_L  = hydrogen.recombination(Ne_L*Nnorm, Te_L*Tnorm)*SQ(Ne_L) * Nnorm / Omega_ci;
            BoutReal R_rc_C  = hydrogen.recombination(Ne_C*Nnorm, Te_C*Tnorm)*SQ(Ne_C) * Nnorm / Omega_ci;
            BoutReal R_rc_R  = hydrogen.recombination(Ne_R*Nnorm, Te_R*Tnorm)*SQ(Ne_R) * Nnorm / Omega_ci;
            
            // Rrec is radiated energy, Erec is energy transferred to neutrals
            // Factor of 1.09 so that recombination becomes an energy source at 5.25eV
            Rrec(i,j,k) = (
                                J_L * (1.09*Te_L - 13.6/Tnorm)*R_rc_L
                           + 4.*J_C * (1.09*Te_C - 13.6/Tnorm)*R_rc_C
                           +    J_R * (1.09*Te_R - 13.6/Tnorm)*R_rc_R
                           ) / (6. * J_C);
            
            Erec(i,j,k) = (3./2) * (
                                         J_L * Te_L * R_rc_L
                                    + 4.*J_C * Te_C * R_rc_C
                                    +    J_R * Te_R * R_rc_R
                                    ) / (6. * J_C);
	    
            Frec(i,j,k) = (
                                 J_L * Vi_L * R_rc_L
                           + 4.* J_C * Vi_C * R_rc_C
                           +     J_R * Vi_R * R_rc_R
                           ) / (6. * J_C);

            Srec(i,j,k) = (
                                 J_L * R_rc_L
                           + 4.* J_C * R_rc_C
                           +     J_R * R_rc_R
                           ) / (6. * J_C);
            
            ///////////////////////////////////////      
            // Ionisation
            BoutReal R_iz_L = Ne_L*Nn_L*hydrogen.ionisation(Te_L*Tnorm) * Nnorm / Omega_ci;
            BoutReal R_iz_C = Ne_C*Nn_C*hydrogen.ionisation(Te_C*Tnorm) * Nnorm / Omega_ci;
            BoutReal R_iz_R = Ne_R*Nn_R*hydrogen.ionisation(Te_R*Tnorm) * Nnorm / Omega_ci;
            
            Riz(i,j,k) = (Eionize/Tnorm) * (    // Energy loss per ionisation
                                                 J_L * R_iz_L
                                            + 4.*J_C * R_iz_C
                                            +    J_R * R_iz_R
                                             ) / (6. * J_C);   
            Eiz(i,j,k) = -(3./2)* (   // Energy from neutral atom temperature
                                          J_L * Tn_L * R_iz_L
                                   + 4. * J_C * Tn_C * R_iz_C
                                   +      J_R * Tn_R * R_iz_R
                                  ) / (6. * J_C);

            // Friction due to ionisation
            Fiz(i,j,k) = - (
                                   J_L * Vn_L * R_iz_L
                            + 4. * J_C * Vn_C * R_iz_C
                            +      J_R * Vn_R * R_iz_R
                            ) / (6. * J_C);
            
            // Plasma sink due to ionisation (negative)
            Siz(i,j,k) = - (
                                 J_L * R_iz_L
                           + 4.* J_C * R_iz_C
                           +     J_R * R_iz_R
                            ) / (6. * J_C);
            
            if(elastic_scattering) {
              /////////////////////////////////////////////////////////
              // Ion-neutral elastic scattering
              //
              // Post "A Review of Recent Developments in Atomic Processes for Divertors and Edge Plasmas" PSI review paper
              //       https://arxiv.org/pdf/plasm-ph/9506003.pdf
              // Relative velocity of two particles in a gas
              // is sqrt(8kT/pi mu) where mu = m_A*m_B/(m_A+m_B)
              // here ions and neutrals have same mass,
              // and the ion temperature is used
              
              BoutReal a0 = 3e-19; // Effective cross-section [m^2]

              // Rates (normalised)
              BoutReal R_el_L = a0*Ne_L*Nn_L*Cs0*sqrt((16./PI)*Te_L) * Nnorm/Omega_ci;
              BoutReal R_el_C = a0*Ne_C*Nn_C*Cs0*sqrt((16./PI)*Te_C) * Nnorm/Omega_ci;
              BoutReal R_el_R = a0*Ne_R*Nn_R*Cs0*sqrt((16./PI)*Te_R) * Nnorm/Omega_ci;

              // Elastic transfer of momentum
              Fel(i,j,k) = (
                               J_L * (Vi_L - Vn_L)*R_el_L
                          + 4.*J_C * (Vi_C - Vn_C)*R_el_C
                          +    J_R * (Vi_R - Vn_R)*R_el_R
                          ) / (6. * J_C);

              // Elastic transfer of thermal energy
              Eel(i,j,k) = (3./2)* (
                                       J_L * (Te_L - Tn_L)*R_el_L
                                  + 4.*J_C * (Te_C - Tn_C)*R_el_C
                                  +    J_R * (Te_R - Tn_R)*R_el_R
                                  ) / (6. * J_C);
            }

            if(excitation) {
              /////////////////////////////////////////////////////////
              // Electron-neutral excitation
              // Note: Rates need checking
              // Currently assuming that quantity calculated is in [eV m^3/s]
              
              BoutReal R_ex_L  = Ne_L*Nn_L*hydrogen.excitation(Te_L*Tnorm) * Nnorm / Omega_ci / Tnorm;
              BoutReal R_ex_C  = Ne_C*Nn_C*hydrogen.excitation(Te_C*Tnorm) * Nnorm / Omega_ci / Tnorm;
              BoutReal R_ex_R  = Ne_R*Nn_R*hydrogen.excitation(Te_R*Tnorm) * Nnorm / Omega_ci / Tnorm;

              Rex(i,j,k) = (
                                 J_L * R_ex_L
                            + 4.*J_C * R_ex_C
                            +    J_R * R_ex_R
                                ) / (6. * J_C); 
            }
            
            
            // Total energy lost from system
            R(i,j,k) = Rzrad(i,j,k)     // Radiated power from impurities
                     + Rrec(i,j,k)      // Recombination
                     + Riz(i,j,k)       // Ionisation
                     + Rex(i,j,k);      // Excitation
            
            // Total energy transferred to neutrals
            E(i,j,k) = Ecx(i,j,k)       // Charge exchange
                     + Erec(i,j,k)      // Recombination
                     + Eiz(i,j,k)       // ionisation
                     + Eel(i,j,k);      // Elastic collisions 

            // Total friction
            F(i,j,k) = Frec(i,j,k)      // Recombination
                     + Fiz(i,j,k)       // Ionisation
                     + Fcx(i,j,k)       // Charge exchange 
                     + Fel(i,j,k);      // Elastic collisions

            // Total sink of plasma, source of neutrals
            S(i,j,k) = Srec(i,j,k) + Siz(i,j,k);
          }
        
      if(!evolve_nvn && neutral_f_pn) {
        // Not evolving neutral momentum
        F = Grad_par(Pn);
      }
    }
    
    ///////////////////////////////////////////////////
    // Plasma model

    /// Density
    TRACE("ddt(Ne)");
    
    if(rhs_explicit) {
      
      if(density_form == 1) {
        // Upwinding fluxes for advection
        ddt(Ne) = 
          - Div_par_FV(Ne, Vi) // Mass flow
          ;
      }else if(density_form == 2) {
        // FE splitting of central difference advection
        ddt(Ne) = 
          - 0.5*( Div_par(NVi) + Vi*Grad_par(Ne) + Ne*Div_par(Vi) )
          ;
      }else if(density_form == 3) {
        // Central differencing for comparison with form 2
        ddt(Ne) = 
          - Div_par(NVi);
      }else if(density_form == 4) {
        // Flux splitting with upwinding
        Field3D a = sqrt( gamma_sound*2.*Te ); // Local sound speed
        ddt(Ne) = 
          - Div_par_FV_FS(Ne, Vi, a, bndry_flux_fix) // Mass flow
          ;
      }else {
        throw BoutException("Unrecognised density_form (%d)", density_form);
      }
      
      if(atomic) {
        ddt(Ne) -= S;                  // Sink to recombination
      }
      
      if(volume_source) {
        ddt(Ne) += NeSource; // External volume source
      }
      
    }else {
      ddt(Ne) = 0.0;
    }

    if(rhs_implicit && (anomalous_D > 0.0)) {
      ddt(Ne) += Div_par_diffusion(anomalous_D, Ne);
    }
    
    if((hyper > 0.0) && (rhs_implicit)) {
      ddt(Ne) += D(Ne, hyper);
    }
    
    if((ADpar > 0.0) && (rhs_implicit)) {
      ddt(Ne) += ADpar * AddedDissipation(1.0, P, Ne, true);
    }
    
    /// Momentum
    TRACE("ddt(NVi)");
    
    if(rhs_explicit) {
      
      if(momentum_form == 1) {
        // Flux form for flow,
        // central differencing of pressure
        
        ddt(NVi) = 
          - Div_par_FV(NVi, Vi) // Momentum flow
          - Grad_par(P)
          ;
      }else if(momentum_form == 3) {
        // FE splitting of flux
        // central differencing pressure
        ddt(NVi) = 
          - 0.5*( Div_par(Vi*NVi) + NVi*Div_par(Vi) + Vi*Grad_par(NVi) )
          - Grad_par(P)
          ;
      }else if(momentum_form == 4) {
        // Flux form for flow,
        // central differencing of pressure
        
        ddt(NVi) = 
          - Div_par_FV3(Ne, Vi, Vi)
          - Grad_par(P)
          ;
      }else if(momentum_form == 5) {
        // Central differencing for comparison
        ddt(NVi) = 
          - Div_par(Vi*NVi)
          - Grad_par(P)
          ;
      }else if(momentum_form == 6) {
        // Flux splitting with upwinding
        Field3D a = sqrt( gamma_sound*2.*Te ); // Local sound speed
        ddt(NVi) = 
          - Div_par_FV_FS(NVi, Vi, a, bndry_flux_fix) // Momentum flow
          - Grad_par(P)
          ;
      }else {
        throw BoutException("Unrecognised momentum_form (%d)", momentum_form);
      }
      
      if(atomic) {
        // Friction with neutrals
        ddt(NVi) -= F;
      }
    }else {
      ddt(NVi) = 0.0; 
    }
    
    if((viscos > 0.) && (rhs_implicit)) {
      ddt(NVi) += viscos*Div_par_diffusion_index(Vi);
    }
    
    if(rhs_implicit && (anomalous_D > 0.0)) {
      ddt(NVi) += Div_par_diffusion(anomalous_D*Vi, Ne);
    }
    
    if(ion_viscosity) {
       // Braginskii ion viscosity
      if(rhs_explicit) {
        // Update viscosity
        
        Field3D tau_i = sqrt(2 * mi_me) * tau_e;
        eta_i = (4./3) * 0.96 * Ne * tau_i * Te;  // Ti = Te
        eta_i.applyBoundary("neumann");
      }
      if(rhs_implicit) {
        ddt(NVi) += Div_par_diffusion(eta_i, Vi);
      }
    }
    
    if((hyper > 0.0) && (rhs_implicit)) {
      ddt(NVi) += D(NVi, hyper);
    }
    
    if((ADpar > 0.0) && (rhs_implicit)) {
      ddt(NVi) += ADpar * AddedDissipation(1.0, P, NVi, true);
    }
    
    /// Pressure

    TRACE("ddt(P)");

    if(rhs_explicit) {
      // Note: ddt(P) set earlier for sheath
      
      if(energy_form == 1) {
        // Upwinding fluxes for advection,
        // central differencing for compression
        ddt(P) += 
          - Div_par_FV(P, Vi)         // Advection
          - (2./3)*P*Div_par(Vi)      // Compression
          ;
      }else if(energy_form == 3) {
        // FE splitting of central difference advection
        // Central differencing compression, grad_par form
        ddt(P) += 
          - (5./3)*0.5*( Div_par(P*Vi) + Vi*Grad_par(P) + P*Div_par(Vi) )
          + (2./3)*Vi*Grad_par(P)
          ;
      }else if(energy_form == 4) {
        // FE splitting, compression in Div_par form 
        ddt(P) += 
          - 0.5*( Div_par(P*Vi) + Vi*Grad_par(P) + P*Div_par(Vi) )
          - (2./3)*P*Div_par(Vi)
          ;
      }else if(energy_form == 5) {
        // Central differencing, symmetric split
        ddt(P) += 
          - (4./3)*Div_par(P*Vi) 
          + (1./3)*(Vi*Grad_par(P) - P*Div_par(Vi))
          ;
      }else if(energy_form == 6) {
        // Upwinding fluxes for advection,
        // central differencing for compression
        // This form uses separate interpolation of Ne and Te 
        ddt(P) += 
          - Div_par_FV3(Ne, 2.*Te, Vi)
          - (2./3)*P*Div_par(Vi)      // Compression
          ;
      }else if(energy_form == 7) {
        // Central differencing for comparison
        ddt(P) += 
          - Div_par(P*Vi)
          - (2./3)*P*Div_par(Vi)
          ;
      }else if(energy_form == 8) {
        Field3D a = sqrt( gamma_sound*2.*Te ); // Local sound speed
        ddt(P) += 
          - Div_par_FV_FS(P, Vi, a, bndry_flux_fix)       // Advection
          - (2./3)*P*Div_par(Vi)          // Compression
          ;
      }else {
        throw BoutException("Unrecognised energy_form (%d)", energy_form);
      }
      
      if(atomic) {
        // Include radiation and neutral interaction
        ddt(P) -= (2./3) * ( 
                            R  // Radiated power
                            +E // Energy transferred to neutrals
                             );
      }
      
      if(volume_source) {
        // Volumetric source
      
        ddt(P) += PeSource;            // External source of energy
      }else {
        // Insert power into the first grid point
        for(RangeIterator r=mesh->iterateBndryLowerY(); !r.isDone(); r++)
          for(int jz=0; jz<mesh->LocalNz; jz++) {
            ddt(P)(r.ind, mesh->ystart, jz) += (2./3)*powerflux / ( coord->dy(r.ind, mesh->ystart) * sqrt(coord->g_22(r.ind, mesh->ystart)));
          }
      }
    }
    
    if (rhs_implicit) {
      if (heat_conduction) {
        if (!nonlocal_conduction) {
          // Standard local heat conduction
          
          //ddt(P) += (2./3)*Div_par_diffusion(kappa_epar, Te);
          ddt(P) += (2./3)*Div_par_diffusion_upwind(kappa_epar, Te);  
          //ddt(P) += (2./3)*Div_par_spitzer( 3.2 * mi_me * Omega_ci * tau_e0, Te);
        } else {
          // Nonlocal heat conduction
          
          // Specify plasma profiles in SI units
          nonlocal_parallel->set_n_electron(Ne*Nnorm);
          nonlocal_parallel->set_T_electron(SI::qe*Te*Tnorm);
          nonlocal_parallel->set_V_electron(Vi*Cs0);
          
          // Set boundary condition on heat flux
          // Calculated from q = gamma e n cs
          Field3D heat_flux_boundary = 0.0;
          for(RangeIterator r=mesh->iterateBndryUpperY(); !r.isDone(); r++) {
            int jz = 0;
            
            // Outward flow velocity to >= Cs
            BoutReal Vout = sqrt(2.0*Te(r.ind, mesh->yend, jz)); // Sound speed outwards
            if(Vi(r.ind, mesh->yend, jz) > Vout)
              Vout = Vi(r.ind, mesh->yend, jz); // If plasma is faster, go to plasma velocity

            BoutReal Nout = 0.5*(Ne(r.ind, mesh->yend, jz) + Ne(r.ind, mesh->yend+1, jz));
            
            BoutReal flux = Nout * Vout;
            BoutReal q = (sheath_gamma - 6) * Te(r.ind, mesh->yend, jz) * flux;
            // Convert to SI units
            heat_flux_boundary(r.ind, mesh->yend, jz) = (Cs0*SI::qe*Tnorm*Nnorm) * q;
          }
          if (not symmetry_plane) {
            for(RangeIterator r=mesh->iterateBndryLowerY(); !r.isDone(); r++) {
              int jz = 0;
              
              // Outward flow velocity to >= Cs
              BoutReal Vout = -sqrt(2.0*Te(r.ind, mesh->ystart, jz)); // Sound speed outwards
              if(Vi(r.ind, mesh->ystart, jz) > Vout)
                Vout = Vi(r.ind, mesh->ystart, jz); // If plasma is faster, go to plasma velocity

              BoutReal Nout = 0.5*(Ne(r.ind, mesh->ystart, jz) + Ne(r.ind, mesh->ystart-1, jz));
              
              BoutReal flux = Nout * Vout;
              BoutReal q = (sheath_gamma - 6) * Te(r.ind, mesh->ystart, jz) * flux;
              // Convert to SI units
              heat_flux_boundary(r.ind, mesh->ystart, jz) = (Cs0*SI::qe*Tnorm*Nnorm) * q;
            }
          } // else: Zero heat flux at lower Y
          
          nonlocal_parallel->set_heat_flux_boundary_condition(heat_flux_boundary);
          
          
          // Specify sources
          //nonlocal_parallel->set_maxwellian_source(0,0);
          
          // Convert parallel metric to SI
          coord->g_22 *= rho_s0*rho_s0; // Convert to m^2
          
          nonlocal_parallel->calculate_nonlocal_closures();

          coord->g_22 /= rho_s0*rho_s0; // Restore g_22

          qe_nonlocal = (1/(Cs0*SI::qe*Tnorm*Nnorm))*nonlocal_parallel->electron_heat_flux;
          qe_local = -kappa_epar*Grad_par(Te);

          // Apply a neumann boundary condition, so that the
          // boundary condition we have imposed on the first and last points
          // is also imposed on the boundaries
          qe_nonlocal.applyBoundary("neumann");
          //qe_nonlocal.applyBoundary("free_o2");
          
          // Note: This is local heat flux
          //ddt(P) += (2./3)*Div_par_diffusion_upwind(kappa_epar, Te);

          // Add divergence of non-local heat flux
          ddt(P) -= (2./3)*Div_par(qe_nonlocal);
        }
      }

      if(anomalous_D > 0.0) {
        ddt(P) += Div_par_diffusion(anomalous_D*2.*Te, Ne);
      }
      if(anomalous_chi > 0.0) {
        ddt(P) += Div_par_diffusion(anomalous_chi, Te);
      }
    }
    
    if((hyper > 0.0) && (rhs_implicit)) {
      ddt(P) += D(P, hyper);
    }
    
    if((ADpar > 0.0) && (rhs_implicit)) {
      ddt(P) += ADpar * AddedDissipation(1.0, P, P, true);
    }

    // Switch off evolution at very low densities
    for(int i=0;i<mesh->LocalNx;i++)
      for(int j=mesh->ystart;j<=mesh->yend;j++)
        for(int k=0;k<mesh->LocalNz;k++) {
          if((Ne(i,j,k) < 1e-5) && (ddt(Ne)(i,j,k) < 0.0)) {
            ddt(Ne)(i,j,k) = 0.0;
            ddt(NVi)(i,j,k) = 0.0;
            ddt(P)(i,j,k) = 0.0;
          }
        }
    
    if(atomic) {
      ///////////////////////////////////////////////////
      // Neutrals model
      //
      // 
      
      TRACE("Neutrals");
      
      Field3D logPn = log(floor(Pn,1e-7));
      logPn.applyBoundary("neumann");
      
      TRACE("ddt(Nn)");
      
      if(rhs_explicit) {
        ddt(Nn) = 
          - Div_par_FV(Nn, Vn)        // Advection
          + S                         // Source from recombining plasma
          - nloss*Nn                  // Loss of neutrals from the system
          ;
      }else {
        ddt(Nn) = 0.0;
      }
      
      if(rhs_implicit) {
        if(dneut > 0.0) {
          ddt(Nn) += Div_par_diffusion(Dn * Nn, logPn); // Diffusion
        }
      }
      
      if((hyper > 0.0) && (rhs_implicit)) {
        ddt(Nn) += D(Nn, hyper);
      }
      
      if(evolve_nvn) {
        // Evolving momentum of the neutral gas
        
        TRACE("ddt(NVn)");
        
        if(rhs_explicit) {
          ddt(NVn) = 
            - Div_par_FV(NVn, Vn)        // Momentum flow
            + F                          // Friction with plasma
            - nloss*NVn                  // Loss of neutrals from the system
            - Grad_par(Pn)               // Pressure gradient
            ;
        }else {
          ddt(NVn) = 0.0;
        }
        
        if(rhs_implicit) {
          if(viscos > 0.) {
            // Note no factor of Nn
            ddt(NVn) += Div_par_diffusion(viscos*SQ(coord->dy), Vn);
          }
          
          if(hyper > 0.) {
            // Numerical dissipation
            ddt(NVn) += D(NVn, hyper);
          }
          
          if(ion_viscosity) {
            // Relationship between heat conduction and viscosity for neutral gas
            // Chapman, Cowling "The Mathematical Theory of Non-Uniform Gases", CUP 1952
            // Ferziger, Kaper "Mathematical Theory of Transport Processes in Gases", 1972
            //
            Field3D eta_n = (2./5) * kappa_n;
            
            ddt(NVn) += Div_par_diffusion(eta_n, Vn);
          }
          
          if(dneut > 0.0)
            ddt(NVn) += Div_par_diffusion(NVn*Dn, logPn); // Diffusion
        }
      }
      
      if(evolve_pn) {
        // Evolving temperature of neutral gas
        // Essentially the same as the plasma equation
        
        TRACE("ddt(Pn)");
        
        if(rhs_explicit) {
          ddt(Pn) +=
            - Div_par_FV(Pn, Vn)                     // Advection
            - (2./3)*Pn*Div_par(Vn)                  // Compression
            + (2./3) * E                             // Energy transferred to neutrals
            - nloss * Pn                             // Loss of neutrals from the system
            ;
        }
        
        if(rhs_implicit) {
          ddt(Pn) += (2./3)*Div_par_diffusion(kappa_n, Tn);  // Parallel heat conduction
          
          if(dneut > 0.0) {
            ddt(Pn) += Div_par_diffusion(Dn*Pn, logPn); // Perpendicular diffusion
          }
        }
        
        if((hyper > 0.0) && (rhs_implicit)) {
          ddt(Pn) += D(Pn, hyper);
        }
        
        // Switch off evolution at very low densities
        // This seems to be necessary to get through initial transients
        
        for(int i=0;i<mesh->LocalNx;i++)
          for(int j=mesh->ystart;j<=mesh->yend;j++)
            for(int k=0;k<mesh->LocalNz;k++) {
              if(Nn(i,j,k) < 1e-5) {
                // Relax to the plasma temperature
                ddt(Pn)(i,j,k) = -1e-2*(Pn(i,j,k) - Te(i,j,k)*Nn(i,j,k));
              }
            }
      }
     
    
      if(rhs_explicit) {
        // Boundary condition on fluxes
        
        TRACE("Fluxes");
        
        BoutReal nredist;
        for(RangeIterator r=mesh->iterateBndryUpperY(); !r.isDone(); r++) {
          int jz = 0; // Z index
          int jy = mesh->yend;
          //flux_ion = 0.0;
          flux_ion = 0.25*(Ne(r.ind, jy, jz) + Ne(r.ind, jy+1, jz)) * (Vi(r.ind, jy, jz) + Vi(r.ind, jy+1, jz)) * (coord->J(r.ind,jy) + coord->J(r.ind,jy+1)) / (sqrt(coord->g_22(r.ind,jy))+ sqrt(coord->g_22(r.ind,jy+1)));
          BoutReal flux_neut = 0.0;
          
          for(int j = mesh->yend+1; j<mesh->LocalNy; j++) {
            //flux_ion += ddt(Ne)(r.ind, j, jz) * coord->J(r.ind,j) * coord->dy(r.ind,j);        
            flux_neut += ddt(Nn)(r.ind, j, jz) * coord->J(r.ind,j) * coord->dy(r.ind,j);
            
            ddt(Ne)(r.ind, j, jz) = 0.0;
            ddt(Nn)(r.ind, j, jz) = 0.0;
          }
          
          // Make sure that mass is conserved
          
          // Total amount of neutral gas to be added
          BoutReal nadd = flux_ion*frecycle + flux_neut + gaspuff;
          
          // Neutral gas arriving at the target
          BoutReal ntarget = (1 - fredistribute) * nadd / ( coord->J(r.ind,mesh->yend) * coord->dy(r.ind,mesh->yend) );
          
          ddt(Nn)(r.ind, mesh->yend, jz) += ntarget;
          
          if(evolve_nvn) {
            // Set velocity of neutrals coming from the wall to a fraction of the 
            // Franck-Condon energy
            BoutReal Vneut = - vwall * sqrt(3.5/Tnorm);
            ddt(NVn)(r.ind, mesh->yend, jz) += ntarget* Vneut;
          }
          
          if(evolve_pn) {
            // Set temperature of the incoming neutrals to F-C
            ddt(Pn)(r.ind, mesh->yend, jz) += ntarget * (3.5/Tnorm);
          }
          
          // Re-distribute neutrals
          nredist = fredistribute * nadd;
          
          // Divide flux_ion by J so that the result in the output file has units of flux per m^2
          flux_ion /= coord->J(mesh->xstart, mesh->yend+1);
        }
        // Now broadcast redistributed neutrals to other processors
        MPI_Comm ycomm = mesh->getYcomm(mesh->xstart); // MPI communicator
        int np; MPI_Comm_size(ycomm, &np); // Number of processors
        
        // Broadcast from final processor (presumably with target)
        // to all other processors
        MPI_Bcast(&nredist, 1, MPI_DOUBLE, np-1, ycomm);
        
        // Distribute along length
        for(int j=mesh->ystart;j<=mesh->yend;j++) {
          // Neutrals into this cell
          BoutReal ncell = nredist * redist_weight(mesh->xstart,j) / ( coord->J(mesh->xstart,j) * coord->dy(mesh->xstart,j) );
          
          ddt(Nn)(mesh->xstart, j, 0) += ncell;
          
          // No momentum
          
          if(evolve_pn) {
            // Set temperature of the incoming neutrals to F-C
            ddt(Pn)(mesh->xstart, j, 0) += ncell * (3.5/Tnorm);
          }
        }
        if (not symmetry_plane) {
          for(RangeIterator r=mesh->iterateBndryLowerY(); !r.isDone(); r++) {
            int jz = 0; // Z index
            int jy = mesh->ystart;
            //flux_ion = 0.0;
            flux_ion = 0.25*(Ne(r.ind, jy, jz) + Ne(r.ind, jy-1, jz)) * (Vi(r.ind, jy, jz) + Vi(r.ind, jy-1, jz)) * (coord->J(r.ind,jy) + coord->J(r.ind,jy-1)) / (sqrt(coord->g_22(r.ind,jy))+ sqrt(coord->g_22(r.ind,jy-1)));
            BoutReal flux_neut = 0.0;
            
            for(int j = mesh->ystart-1; j>0; j--) {
              //flux_ion += ddt(Ne)(r.ind, j, jz) * coord->J(r.ind,j) * coord->dy(r.ind,j);        
              flux_neut += ddt(Nn)(r.ind, j, jz) * coord->J(r.ind,j) * coord->dy(r.ind,j);
              
              ddt(Ne)(r.ind, j, jz) = 0.0;
              ddt(Nn)(r.ind, j, jz) = 0.0;
            }
            
            // Make sure that mass is conserved
            
            // Total amount of neutral gas to be added
            BoutReal nadd = flux_ion*frecycle + flux_neut + gaspuff;
            
            // Neutral gas arriving at the target
            BoutReal ntarget = (1 - fredistribute) * nadd / ( coord->J(r.ind,mesh->ystart) * coord->dy(r.ind,mesh->ystart) );
            
            ddt(Nn)(r.ind, mesh->ystart, jz) += ntarget;
            
            if(evolve_nvn) {
              // Set velocity of neutrals coming from the wall to a fraction of the 
              // Franck-Condon energy
              BoutReal Vneut = vwall * sqrt(3.5/Tnorm);
              ddt(NVn)(r.ind, mesh->ystart, jz) += ntarget* Vneut;
            }
            
            if(evolve_pn) {
              // Set temperature of the incoming neutrals to F-C
              ddt(Pn)(r.ind, mesh->ystart, jz) += ntarget * (3.5/Tnorm);
            }
            
            // Re-distribute neutrals
            nredist = fredistribute * nadd;
            
            // Divide flux_ion by J so that the result in the output file has units of flux per m^2
            flux_ion /= coord->J(mesh->xstart, mesh->ystart-1);
          }
          // Now broadcast redistributed neutrals to other processors
          MPI_Comm ycomm = mesh->getYcomm(mesh->xstart); // MPI communicator
          
          // Broadcast from first processor (presumably with target)
          // to all other processors
          MPI_Bcast(&nredist, 1, MPI_DOUBLE, 0, ycomm);
          
          // Distribute along length
          for(int j=mesh->ystart;j<=mesh->yend;j++) {
            // Neutrals into this cell
            BoutReal ncell = nredist * redist_weight(mesh->xstart,j) / ( coord->J(mesh->xstart,j) * coord->dy(mesh->xstart,j) );
            
            ddt(Nn)(mesh->xstart, j, 0) += ncell;
            
            // No momentum
            
            if(evolve_pn) {
              // Set temperature of the incoming neutrals to F-C
              ddt(Pn)(mesh->xstart, j, 0) += ncell * (3.5/Tnorm);
            }
          }
        }
      }
    }
    return 0;
  }

  /*!
   * Preconditioner. Solves the heat conduction
   * 
   * @param[in] t  The simulation time
   * @param[in] gamma   Factor in front of the Jacobian in (I - gamma*J). Related to timestep
   * @param[in] delta   Not used here
   */
  int precon(BoutReal t, BoutReal gamma, BoutReal delta) {
    
    static InvertPar *inv = NULL;
    if(!inv) {
      // Initialise parallel inversion class
      inv = InvertPar::Create();
      inv->setCoefA(1.0);
    }
    if(heat_conduction) {
      // Set the coefficient in front of Grad2_par2
      inv->setCoefB(-(2./3)*gamma*kappa_epar);
      Field3D dT = ddt(P);
      dT.applyBoundary("neumann");
      ddt(P) = inv->solve(dT);
    }
    
    if (atomic) {
      if (evolve_pn) {
        // Neutral pressure
        inv->setCoefB(-(2./3)*gamma*kappa_n);
        Field3D dT = ddt(Pn);
        dT.applyBoundary("neumann");
        ddt(Pn) = inv->solve(dT);
      }
      
      if (dneut > 0.0) {
        inv->setCoefB(-gamma*Dn);
        Field3D tmp = ddt(Nn);
        tmp.applyBoundary("neumann");
        ddt(Nn) = inv->solve(tmp);
      }
    }

    return 0;
  }

  /*!
   * When split operator is enabled, run only the explicit part
   */
  int convective(BoutReal t) {
    rhs_explicit = true;
    rhs_implicit = false;
    update_coefficients = true;
    return rhs(t);
  }

  /*!
   * When split operator is enabled, run only implicit part
   */
  int diffusive(BoutReal t, bool linear) {
    rhs_explicit = false;
    rhs_implicit = true;
    update_coefficients = !linear; // Don't update coefficients in linear solve
    return rhs(t);
  }
  
  /*!
   * Monitor output solutions
   */
  int outputMonitor(BoutReal simtime, int iter, int NOUT) {
    
    static BoutReal maxinvdt_alltime = 0.0; // Max 1/dt over all output times
    
    ///////////////////////////////////////////////////
    // Check velocities for CFL information
    
    if(cfl_info) {
      // Calculate the maximum velocity, including cell centres
      // and edges.
      
      Coordinates *coord = mesh->coordinates();

      BoutReal maxabsvc = 0.0; // Maximum absolute velocity + sound speed
      BoutReal maxinvdt = 0.0; // Maximum 1/dt
      for(int j=mesh->ystart;j<=mesh->yend;j++) {
        BoutReal g = 5./3;
        
        // cell centre
        BoutReal cs = sqrt( g*P(0,j,0)/Ne(0,j,0) ); // Sound speed
        
        BoutReal vcs = abs( Vi(0,j,0) ) + cs;
        if(vcs > maxabsvc)
          maxabsvc = vcs;
        
        BoutReal dl = coord->dy(0,j) * sqrt(coord->g_22(0,j)); // Length of cell
        if( vcs/dl > maxinvdt )
          maxinvdt = vcs/dl;
        
        // cell left
        BoutReal p = 0.5*(P(0,j-1,0) + P(0,j,0));
        BoutReal n = 0.5*(Ne(0,j-1,0) + Ne(0,j,0));
        cs = sqrt( g*p/n );
        vcs = abs( 0.5*(Vi(0,j-1,0) + Vi(0,j,0)) ) + cs;
        if(vcs > maxabsvc)
          maxabsvc = vcs;
        
        dl = 0.5*(coord->dy(0,j) * sqrt(coord->g_22(0,j)) + coord->dy(0,j-1) * sqrt(coord->g_22(0,j-1)));
        
        if( vcs/dl > maxinvdt )
          maxinvdt = vcs/dl;
        
        // Cell right
        p = 0.5*(P(0,j+1,0) + P(0,j,0));
        n = 0.5*(Ne(0,j+1,0) + Ne(0,j,0));
        cs = sqrt( g*p/n );
        vcs = abs( 0.5*(Vi(0,j+1,0) + Vi(0,j,0)) ) + cs;
        if(vcs > maxabsvc)
          maxabsvc = vcs;

        dl = 0.5*(coord->dy(0,j) * sqrt(coord->g_22(0,j)) + coord->dy(0,j+1) * sqrt(coord->g_22(0,j+1)));
        
        if( vcs/dl > maxinvdt )
          maxinvdt = vcs/dl;
      }
      
      // Get maximum over the domain
      BoutReal maxabsvc_all;
      BoutReal maxinvdt_all;

      MPI_Allreduce(&maxabsvc, &maxabsvc_all, 1, MPI_DOUBLE, MPI_MAX, BoutComm::get());
      MPI_Allreduce(&maxinvdt, &maxinvdt_all, 1, MPI_DOUBLE, MPI_MAX, BoutComm::get());
      
      if(maxinvdt_all > maxinvdt_alltime)
        maxinvdt_alltime = maxinvdt_all;

      output.write("\nLocal max |v|+cs: %e Global max |v|+cs: %e\n", maxabsvc, maxabsvc_all);
      output.write("Local CFL limit: %e Global limit: %e\n", 1./maxinvdt, 1./maxinvdt_all);
      output.write("Minimum global CFL limit %e\n", 1./maxinvdt_alltime);
    }
    return 0;
  }

private:
  bool cfl_info; // Print additional information on CFL limits
  
  // Normalisation parameters
  BoutReal Tnorm, Nnorm, Bnorm, AA;
  BoutReal Cs0, Omega_ci, rho_s0, tau_e0, mi_me;

  /////////////////////////////////////////////////////////////////
  // Evolving quantities
  Field3D Ne, NVi, P; // Plasma (electron) density, momentum, and pressure
  Field3D Nn, NVn, Pn; // Neutral density, momentum, pressure
  
  Field3D Vi, Vn;  // Ion and neutral velocities
  
  bool evolve_nvn; // Evolve neutral momentum?
  bool evolve_pn;  // Evolve neutral pressure?
  
  /////////////////////////////////////////////////////////////////
  // Diffusion and viscosity coefficients
  
  Field3D Dn;           // Neutral gas diffusion
  BoutReal dneut;       // Neutral gas diffusion multiplier
  
  Field3D kappa_n;      // Neutral gas thermal conduction
  Field3D kappa_epar;   // Plasma thermal conduction

  Field3D tau_e;        // Electron collision time
  Field3D eta_i;        // Braginskii ion viscosity
  bool ion_viscosity;   // Braginskii ion viscosity on/off
  bool heat_conduction; // Thermal conduction on/off
  bool elastic_scattering; // Ion-neutral elastic scattering
  bool excitation;      // Include electron-neutral excitation
  
  int density_form;     // Form of the density equation
  int momentum_form;    // Form of the momentum equation
  int energy_form;      // Form of the energy equation
  
  BoutReal nloss;       // Neutral loss rate (1/timescale)

  BoutReal anomalous_D, anomalous_chi; // Anomalous transport
  
  bool nonlocal_conduction; // Turn on non-local heat conduction?
  NonLocalParallel* nonlocal_parallel;
  Field3D qe_nonlocal, qe_local; // Output diagnostics

  /////////////////////////////////////////////////////////////////
  // Atomic physics transfer channels
  
  bool atomic;     // Include atomic physics? This includes neutral gas evolution
  
  Field3D Srec, Siz;        // Plasma particle sinks due to recombination and ionisation
  Field3D Frec, Fiz, Fcx, Fel;   // Plasma momentum sinks due to recombination, ionisation, charge exchange and elastic collisions
  Field3D Rrec, Riz, Rzrad, Rex; // Plasma power sinks due to recombination, ionisation, impurity radiation, and hydrogen excitation
  Field3D Erec, Eiz, Ecx, Eel;   // Transfer of power from plasma to neutrals
  
  Field3D S, F, E; // Exchange of particles, momentum and energy from plasma to neutrals
  Field3D R;       // Radiated power
  
  RadiatedPower *rad;            // Impurity atomic rates
  UpdatedRadiatedPower hydrogen; // Atomic rates
  
  BoutReal fimp;     // Impurity fraction (of Ne)
  BoutReal Eionize;  // Ionisation energy loss
  
  bool neutral_f_pn; // When not evolving NVn, use F = Grad_par(Pn)

  ///////////////////////////////////////////////////////////////
  // Sheath boundary

  bool symmetry_plane;
  int mid_point1 = -1;
  int mid_point_proc1 = -1;
  int mid_point2 = -1;
  int mid_point_proc2 = -1;
  
  BoutReal sheath_gamma;   // Sheath heat transmission factor
  BoutReal neutral_gamma;  // Neutral heat transmission
  
  int density_sheath; // How to handle density boundary?
  int pressure_sheath; // How to handle pressure boundary?

  bool bndry_flux_fix; 

  BoutReal frecycle; // Recycling fraction
  BoutReal gaspuff;  // Additional source of neutral gas at the target plate
  BoutReal vwall;    // Velocity of neutrals coming from the wall
                     // as fraction of Franck-Condon energy

  BoutReal flux_ion; // Flux of ions to target (output)

  // Re-distribution of recycled neutrals
  Field2D redist_weight; // Weighting used to decide redistribution
  BoutReal fredistribute; // Fraction of recycled neutrals re-distributed along length
  
  ///////////////////////////////////////////////////////////////
  // Sources
  
  bool volume_source;         // Include volume sources?
  Field2D NeSource, PeSource; // Volume sources
  Field2D NeSource0; // Used in feedback control
  BoutReal powerflux; // Used if no volume sources
  
  // Upstream density controller
  BoutReal density_upstream;   // The desired density at the midplane (either lower Y (upstream) boundary, or middle of the domain)
  BoutReal density_controller_p, density_controller_i; // Controller settings
  bool density_integral_positive; // Limit the i term to be positive
  bool density_source_positive; // Limit the source to be positive

  BoutReal density_error_lasttime, density_error_last; // Value and time of last error
  BoutReal density_error_integral; // Integral of error
  BoutReal density_error_lasttime2, density_error_last2; // Value and time of last error at second mid-plane point
  BoutReal density_error_integral2; // Integral of error at second mid-plane point
  
  ///////////////////////////////////////////////////////////////
  // Numerical dissipation
  
  BoutReal tn_floor; // Minimum neutral gas temperature [eV]
  
  BoutReal hyper, viscos;     // Numerical dissipation terms
  BoutReal ADpar;             // Added Dissipation numerical term
  
  Field2D dy4;   // SQ(SQ(coord->dy)) cached to avoid recalculating
  
  BoutReal gamma_sound; // Ratio of specific heats in numerical dissipation term

  // Numerical diffusion
  const Field3D D(const Field3D &f, BoutReal d) {
    if(d < 0.0)
      return 0.0;
    return Div_par_diffusion(d*SQ(mesh->coordinates()->dy), f);
    //return -D4DY4_FV(d*dy4,f);
  }
  
  ///////////////////////////////////////////////////////////////
  // Splitting into implicit and explicit
  bool rhs_implicit, rhs_explicit;   // Enable implicit and explicit parts
  bool update_coefficients;          // Re-calculate diffusion coefficients
};

BOUTMAIN(SD1D);
