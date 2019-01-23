/*
   SD1D: 1D simulation of plasma-neutral interactions
   ==================================================

     Copyright B.Dudson, University of York, 2016-2018
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

#include <bout/constants.hxx>
#include <bout/physicsmodel.hxx>
#include <derivs.hxx>
#include "field_factory.hxx"
#include <invert_parderiv.hxx>
#include "unused.hxx"

#include "div_ops.hxx"
#include "loadmetric.hxx"
#include "radiation.hxx"
#include "species.hxx"
#include "reaction.hxx"

class SD1D : public PhysicsModel {
protected:
  int init(bool restarting) {
    Options *opt = Options::getRoot()->getSection("sd1d");

    // Normalisation
    OPTION(opt, Tnorm, 100);             // Reference temperature [eV]
    OPTION(opt, Nnorm, 1e19);            // Reference density [m^-3]
    OPTION(opt, Bnorm, 1.0);             // Reference magnetic field [T]
    SAVE_ONCE3(Tnorm, Nnorm, Bnorm);      // Save normalisations

    // Model parameters
    OPTION(opt, vwall, 1. / 3); // 1/3rd Franck-Condon energy at wall
    OPTION(opt, frecycle, 1.0); // Recycling fraction 100%
    OPTION(opt, fredistribute,
           0.0); // Fraction of neutrals redistributed evenly along leg
    OPTION(opt, density_sheath, 0);   // Free boundary
    OPTION(opt, pressure_sheath, 0);  // Free boundary
    OPTION(opt, gaspuff, 0.0);        // Additional gas flux at target
    OPTION(opt, dneut, 1.0);          // Scale neutral gas diffusion
    OPTION(opt, nloss, 0.0);          // Neutral gas loss rate
    OPTION(opt, sheath_gamma, 6.5);   // Sheath heat transmission
    OPTION(opt, neutral_gamma, 0.25); // Neutral heat transmission

    if (sheath_gamma < 6)
      throw BoutException("sheath_gamma < 6 not consistent");

    OPTION(opt, tn_floor, 3.5); // Minimum neutral gas temperature [eV]

    OPTION(opt, atomic, true);

    OPTION(opt, hyper, -1);             // Numerical hyper-diffusion
    OPTION(opt, ADpar, -1);             // Added Dissipation scheme
    OPTION(opt, viscos, -1);            // Parallel viscosity
    OPTION(opt, ion_viscosity, false);  // Braginskii parallel viscosity
    OPTION(opt, heat_conduction, true); // Spitzer-Hahm heat conduction
    
    if (atomic) {
      // Include atomic rates
      
      Options::getRoot()->getSection("NVn")->get("evolve", evolve_nvn, true);
      Options::getRoot()->getSection("Pn")->get("evolve", evolve_pn, true);
      
      bool elastic_scattering; // Ion-neutral elastic scattering
      OPTION(opt, elastic_scattering, false);
      if (elastic_scattering) {
        // Add elastic scattering to the reactions set
        reactions.push_back(ReactionFactory::getInstance().create("elastic", opt));
      }
      
      bool recombination;
      OPTION(opt, recombination, true);
      if (recombination) {
        reactions.push_back(ReactionFactory::getInstance().create("recombination", opt));
      }
      
      // Ionisation plasma particle source. Doesn't affect neutral diffusion
      bool ionisation; 
      OPTION(opt, ionisation, true);
      if (ionisation) {
        reactions.push_back(ReactionFactory::getInstance().create("ionisation", opt));
      }

      // Include electron-neutral excitation
      bool excitation;
      OPTION(opt, excitation, false);
      if (excitation) {
        reactions.push_back(ReactionFactory::getInstance().create("excitation", opt));
      }

      bool charge_exchange;
      OPTION(opt, charge_exchange, true);
      if (charge_exchange) {
        reactions.push_back(ReactionFactory::getInstance().create("hydrogen_cx", opt));
      }

      bool neutral_f_pn; // When not evolving NVn, use F = Grad_par(Pn)
      OPTION(opt, neutral_f_pn, true);
      if (!evolve_nvn && neutral_f_pn) {
        // Not evolving neutral momentum. Add force calculated from neutral pressure
        reactions.push_back(ReactionFactory::getInstance().create("neutralpressureforce", opt));
      }

    }
    
    OPTION(opt, gamma_sound, 5. / 3); // Ratio of specific heats
    OPTION(opt, bndry_flux_fix, false);

    // Read the flux-tube area from input file
    // This goes into the Jacobian.
    std::string area_string;
    FieldFactory ffact(mesh);

    // Calculate normalisation factors
    // Note: These are calculated for Hydrogen, with
    // factors of atomic mass and charge included in the equations
    
    Cs0 = sqrt(SI::qe * Tnorm / SI::Mp); // Reference sound speed [m/s]
    Omega_ci = SI::qe * Bnorm / SI::Mp;  // Ion cyclotron frequency [1/s]
    rho_s0 = Cs0 / Omega_ci;                    // Length scale [m]

    mi_me = SI::Mp / SI::Me;

    BoutReal Coulomb = 6.6 - 0.5 * log(Nnorm * 1e-20) + 1.5 * log(Tnorm);
    tau_e0 = 1. / (2.91e-6 * (Nnorm / 1e6) * Coulomb * pow(Tnorm, -3. / 2));

    // Save normalisation factors
    SAVE_ONCE5(Cs0, Omega_ci, rho_s0, tau_e0, mi_me);
    
    nloss /= Omega_ci;

    // Electrons are handled differently,
    // so not an evolving species
    species["e"] = new Species(SI::Me / SI::Mp, -1);

    // Save some electron properties
    auto &electrons = *species["e"];
    dump.addRepeat(electrons.N, "Ne");
    dump.addRepeat(electrons.T, "Te");
    
    // Hydrogen ions, evolving
    species["h+"] = new FluidSpecies("h+",
                                     opt,
                                     solver,
                                     restart,
                                     restarting,
                                     Nnorm, Tnorm, Omega_ci, Cs0); 
    
    if (atomic) {

      // Atomic hydrogen species
      species["h"] = new Species();
      
      solver->add(Nn, "Nn");
      if (evolve_nvn) {
        solver->add(NVn, "NVn");
      }
      if (evolve_pn) {
        solver->add(Pn, "Pn");
      }
    }

    // Load the metric tensor
    LoadMetric(rho_s0, Bnorm);

    opt->get("area", area_string, "1.0");
    mesh->getCoordinates()->J = ffact.create2D(area_string, Options::getRoot());

    dy4 = SQ(SQ(mesh->getCoordinates()->dy));

    //////////////////////////////////////////////////
    // Impurities
    OPTION(opt, fimp, 0.0); // Fixed impurity fraction

    if (fimp > 0) {

      OPTION(opt, impurity_adas, false);
      if (impurity_adas) {
        // Use OpenADAS data through Atomicpp
        // Find out which species to model

        OPTION(opt, impurity_species, "c");
        
        reactions.push_back(
                            ReactionFactory::getInstance().create("atomic++coronal", opt));
      } else {
        // Use carbon radiation for the impurity
        impurity_species = "c";
        reactions.push_back(
                            ReactionFactory::getInstance().create("c_hutchinson", opt));
      }

      // Create a Species object for the impurity
      species[impurity_species] = new Species();
    }

    // Add extra quantities to be saved
    if (atomic) {
      SAVE_REPEAT2(Dn, kappa_n); // Neutral diffusion coefficients
      SAVE_REPEAT(flux_ion);     // Flux of ions to target
    }
    if (heat_conduction)
      SAVE_REPEAT(kappa_epar); // Save coefficient of thermal conduction

    bool diagnose;
    OPTION(opt, diagnose, true);
    if (diagnose) {
      // Output extra variables
      if (atomic) {
        if (evolve_nvn) {
          SAVE_REPEAT(Vn);
        }
      }
    }

    kappa_epar = 0.0;
    
    flux_ion = 0.0;
    
    // Neutral gas diffusion and heat conduction
    Dn = 0.0;
    kappa_n = 0.0;

    // Calculate neutral gas redistribution weights over the domain
    std::string redist_string;
    opt->get("redist_weight", redist_string, "1.0");
    redist_weight = ffact.create2D(redist_string, opt);
    BoutReal localweight = 0.0;
    Coordinates *coord = mesh->getCoordinates();
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      localweight += redist_weight(mesh->xstart, j) *
                     coord->J(mesh->xstart, j) * coord->dy(mesh->xstart, j);
    }

    MPI_Comm ycomm = mesh->getYcomm(mesh->xstart); // MPI communicator

    // Calculate total weight by summing over all processors
    BoutReal totalweight;
    MPI_Allreduce(&localweight, &totalweight, 1, MPI_DOUBLE, MPI_SUM, ycomm);
    // Normalise redist_weight so sum over domain:
    //
    // sum ( redist_weight * J * dy ) = 1
    //
    redist_weight /= totalweight;

    setPrecon((preconfunc)&SD1D::precon);

    //////////////////////////////////////////
    // Split operator (IMEX) schemes
    // Use combination of explicit and implicit methods
    //
    // Boolean flags rhs_explicit and rhs_implicit
    // turn on explicit and implicit terms

    bool split_operator;
    OPTION(opt, split_operator, false);
    if (!split_operator) {
      // Turn on all terms in rhs
      rhs_explicit = rhs_implicit = true;
      update_coefficients = true;
    }
    setSplitOperator(split_operator);


    // Molecules
    // species["h2"] = new FluidSpecies("h2",
    //                                  opt,
    //                                  solver,
    //                                  restarting);
    
    output << "\n-----------------------\nAvailable reactions: \n";
    for (auto i :  ReactionFactory::getInstance().listAvailable()) {
      output << "\t" << i << "\n";
    }
    output << "-------------------------\nEnabled reactions: \n";
    if (reactions.empty()) {
      output << "No reactions\n";
    } else {
      for (const auto &r : reactions) {
        output << "\t" << r->str() << "\n";
      }
    }
    output << "-------------------------\n";

    output << "Species:\n";
    for (auto const &s : species) {
      output << "\t" << s.first << "\n";
    }
    output << "-------------------------\n";
    
    return 0;
  }

  /*!
   * This function calculates the time derivatives
   * of all evolving quantities
   *
   */
  int rhs(BoutReal time) {

    // Evolve ion species
    for(auto &s : species) {
      s.second->evolve(time);
    }

    // Electrons handled separately
    auto &ions = *species.at("h+");
    
    Coordinates *coord = mesh->getCoordinates();
    
    // Set electron species properties
    auto &electrons = *species.at("e");
    electrons.T = ions.T;
    electrons.N = ions.N * ions.ZZ;
    electrons.P = ions.P * ions.ZZ;
    electrons.V = ions.V;
    
    Field3D Nnlim;
    Field3D Tn;
    
    if (atomic) {
      // Includes atomic processes, neutral gas
      mesh->communicate(Nn);
      if (evolve_nvn) {
        mesh->communicate(NVn);
      }
      if (evolve_pn) {
        mesh->communicate(Pn);
      }
      Nn = floor(Nn, 1e-10);
      Nnlim = floor(Nn, 1e-5);

      if (evolve_nvn) {
        Vn = NVn / Nnlim;
      } else {
        Vn = -vwall * sqrt(3.5 / Tnorm);
        NVn = Nn * Vn;
      }

      if (evolve_pn) {
        Tn = Pn / Nnlim;
        // Tn = floor(Tn, 0.025/Tnorm); // Minimum tn_floor
        Tn = floor(Tn, 1e-12);
      } else {
        Tn = ions.T; // Strong CX coupling
        Pn = Tn * floor(Nn, 0.0);
        Tn = floor(Tn, tn_floor / Tnorm); // Minimum of tn_floor
      }

      // Neutral atom species
      auto &atoms = *species.at("h");
      atoms.T = Tn;
      atoms.N = Nn;
      atoms.P = Pn;
      atoms.V = Vn;
    }

    if (update_coefficients) {
      // Update diffusion coefficients
      TRACE("Update coefficients");

      tau_e = Omega_ci * tau_e0 * pow(electrons.T, 1.5) / electrons.N;

      if (heat_conduction) {
        kappa_epar = 3.2 * mi_me * 0.5 * electrons.P * tau_e;
        kappa_epar.applyBoundary("neumann");
      }

      if (atomic) {
        // Neutral diffusion rate

        Field3D Nelim = floor(electrons.N, 1e-5);
        
        for (int i = 0; i < mesh->LocalNx; i++)
          for (int j = 0; j < mesh->LocalNy; j++)
            for (int k = 0; k < mesh->LocalNz; k++) {
              // Charge exchange frequency, normalised to ion cyclotron
              // frequency
              BoutReal sigma_cx = Nelim(i, j, k) * Nnorm *
                                  hydrogen.chargeExchange(ions.T(i, j, k) * Tnorm) /
                                  Omega_ci;

              // Ionisation frequency
              BoutReal sigma_iz = Nelim(i, j, k) * Nnorm *
                                  hydrogen.ionisation(electrons.T(i, j, k) * Tnorm) /
                                  Omega_ci;

              // Neutral thermal velocity
              BoutReal tn = Tn(i, j, k);
              if (tn < tn_floor / Tnorm)
                tn = tn_floor / Tnorm;
              BoutReal vth_n = sqrt(tn); // Normalised to Cs0

              // Neutral-neutral mean free path
              BoutReal Lmax = 0.1; // meters
              BoutReal a0 = PI * SQ(5.29e-11);
              BoutReal lambda_nn = 1. / (Nnorm * Nnlim(i, j, k) * a0); // meters
              if (lambda_nn > Lmax) {
                // Limit maximum mean free path
                lambda_nn = Lmax;
              }

              lambda_nn /= rho_s0; // Normalised length to rho_s0
              // Neutral-Neutral collision rate
              BoutReal sigma_nn = vth_n / lambda_nn;

              // Total neutral collision frequency
              BoutReal sigma = sigma_cx + sigma_iz + sigma_nn;

              // Neutral gas diffusion
              if (dneut > 0.0) {
                Dn(i, j, k) = dneut * SQ(vth_n) / sigma;
              }

              // Neutral gas heat conduction
              kappa_n(i, j, k) = dneut * Nnlim(i, j, k) * SQ(vth_n) / sigma;
            }

        kappa_n.applyBoundary("Neumann");
        Dn.applyBoundary("dirichlet_o2");
        mesh->communicate(kappa_n, Dn);
      }
    }

    // Set sheath boundary condition on flow

    TRACE("Sheath");

    if (evolve_pn) {
      ddt(Pn) = 0.0;
    }

    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      int jz = 0;

      // Set boundary half-way between cells
      for (int jy = mesh->yend + 1; jy < mesh->LocalNy; jy++) {
        
        if (atomic) {
          ///// Neutral model
          // Flux of neutral particles, momentum, and energy are set later
          // Here the neutral velocity is set to no-flow conditions

          // Vn fixed value (Dirichlet)
          Vn(r.ind, jy, jz) = -Vn(r.ind, mesh->yend, jz);

          // Nn free boundary (constant gradient)
          Nn(r.ind, jy, jz) =
              2. * Nn(r.ind, mesh->yend, jz) - Nn(r.ind, mesh->yend - 1, jz);

          if (evolve_pn) {
            // Tn fixed value (Dirichlet)
            // Tn(r.ind, jy, jz) = 3.5/Tnorm - Tn(r.ind, mesh->yend, jz);

            // Tn zero gradient. Heat flux set by gamma
            Tn(r.ind, jy, jz) = Tn(r.ind, mesh->yend, jz);

            if (rhs_explicit && (neutral_gamma > 0.0)) {
              // Density at the target
              BoutReal Nnout = 0.5 * (Nn(r.ind, mesh->yend, jz) +
                                      Nn(r.ind, mesh->yend + 1, jz));
              // gamma * n * T * cs
              BoutReal q = neutral_gamma * Nnout * Tn(r.ind, jy, jz) *
                           sqrt(Tn(r.ind, jy, jz));

              // Multiply by cell area to get power
              BoutReal heatflux = q *
                                  (coord->J(r.ind, mesh->yend) +
                                   coord->J(r.ind, mesh->yend + 1)) /
                                  (sqrt(coord->g_22(r.ind, mesh->yend)) +
                                   sqrt(coord->g_22(r.ind, mesh->yend + 1)));

              // Divide by volume of cell, and 2/3 to get pressure
              ddt(Pn)(r.ind, mesh->yend, jz) -=
                  (2. / 3) * heatflux /
                  (coord->dy(r.ind, mesh->yend) * coord->J(r.ind, mesh->yend));
            }
          } else {
            Tn(r.ind, jy, jz) = ions.T(r.ind, jy, jz);
          }
          Pn(r.ind, jy, jz) = Nn(r.ind, jy, jz) * Tn(r.ind, jy, jz);
          NVn(r.ind, jy, jz) = -NVn(r.ind, mesh->yend, jz);
        }
      }
    }

    if (atomic) {
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        // No-flow boundary condition on left boundary
        
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          for (int jy = 0; jy < mesh->ystart; jy++) {
            Vn(r.ind, jy, jz) = -Vn(r.ind, mesh->ystart, jz);
            Nn(r.ind, jy, jz) = Nn(r.ind, jy, jz);
            Pn(r.ind, jy, jz) = Pn(r.ind, jy, jz);
            Tn(r.ind, jy, jz) = Tn(r.ind, jy, jz);
          }
        }
      }
    }
    
    if (atomic && rhs_explicit) {
      // Atomic physics
      TRACE("Atomic");
      
      if (fimp > 0.0) {
        // Fixed fraction impurity radiation 
        species.at(impurity_species)->N = fimp * ions.N;
      }
      
      // Calculate reactions
      for (auto &r : reactions) {
        TRACE("Reaction %s", r->str().c_str());
        
        r->updateSpecies(species, Tnorm, Nnorm, Cs0, Omega_ci);

        // Particle (density) sources
        for (const auto &s : r->densitySources()) {
          // s.first contains the species label (std::string)
          // s.second is a Field3D with the source in normalised units
          try {
            ddt(species.at(s.first)->N) += s.second;
          } catch (const std::out_of_range &e) {
            throw BoutException("Unhandled density source for species '%s'", s.first.c_str());
          }
        }
        
        // Momentum sources
        for (const auto &s : r->momentumSources()) {
          // Note: Mass should be accounted for somewhere
          try {
            ddt(species.at(s.first)->NV) += s.second;
          } catch (const std::out_of_range &e) {
            throw BoutException("Unhandled momentum source for species '%s'", s.first.c_str());
          }
        }   
        
        // Energy sources
        for (const auto &s : r->energySources()) {
          // Add power to pressure equation, with 2/3 factor
          try {
            ddt(species.at(s.first)->P) += (2./3) * s.second;
          } catch (const std::out_of_range &e) {
            throw BoutException("Unhandled energy source for species '%s'", s.first.c_str());
          }
        }
      }
    }

    ///////////////////////////////////////////////////
    // electron model
    
    {
      TRACE("Electron pressure");
      
      if (heat_conduction) {
        ddt(electrons.P) += (2. / 3) * Div_par_diffusion_upwind(kappa_epar, electrons.T);
      }

      // Electron pressure acts on ions
      ddt(ions.NV) -= Grad_par(ions.P);
    }
    
    if (atomic) {
      ///////////////////////////////////////////////////
      // Neutrals model
      //
      //

      TRACE("Neutrals");

      Field3D logPn = log(floor(Pn, 1e-7));
      logPn.applyBoundary("neumann");

      TRACE("ddt(Nn)");

      if (rhs_explicit) {
        ddt(Nn) = -Div_par_FV(Nn, Vn) // Advection
                  - nloss * Nn        // Loss of neutrals from the system
            ;

      } else {
        ddt(Nn) = 0.0;
      }

      if (rhs_implicit) {
        if (dneut > 0.0) {
          ddt(Nn) += Div_par_diffusion(Dn * Nn, logPn); // Diffusion
        }
      }

      if ((hyper > 0.0) && (rhs_implicit)) {
        ddt(Nn) += D(Nn, hyper);
      }

      if (evolve_nvn) {
        // Evolving momentum of the neutral gas

        TRACE("ddt(NVn)");

        if (rhs_explicit) {
          ddt(NVn) = -Div_par_FV(NVn, Vn) // Momentum flow
                     - nloss * NVn        // Loss of neutrals from the system
                     - Grad_par(Pn)       // Pressure gradient
              ;
        } else {
          ddt(NVn) = 0.0;
        }

        if (rhs_implicit) {
          if (viscos > 0.) {
            // Note no factor of Nn
            ddt(NVn) += Div_par_diffusion(viscos * SQ(coord->dy), Vn);
          }

          if (hyper > 0.) {
            // Numerical dissipation
            ddt(NVn) += D(NVn, hyper);
          }

          if (ion_viscosity) {
            // Relationship between heat conduction and viscosity for neutral
            // gas Chapman, Cowling "The Mathematical Theory of Non-Uniform
            // Gases", CUP 1952 Ferziger, Kaper "Mathematical Theory of
            // Transport Processes in Gases", 1972
            //
            Field3D eta_n = (2. / 5) * kappa_n;

            ddt(NVn) += Div_par_diffusion(eta_n, Vn);
          }

          if (dneut > 0.0)
            ddt(NVn) += Div_par_diffusion(NVn * Dn, logPn); // Diffusion
        }
      }

      if (evolve_pn) {
        // Evolving temperature of neutral gas
        // Essentially the same as the plasma equation

        TRACE("ddt(Pn)");

        if (rhs_explicit) {
          ddt(Pn) += -Div_par_FV(Pn, Vn)           // Advection
                     - (2. / 3) * Pn * Div_par(Vn) // Compression
                     - nloss * Pn   // Loss of neutrals from the system
              ;
        }

        if (rhs_implicit) {
          ddt(Pn) += (2. / 3) *
                     Div_par_diffusion(kappa_n, Tn); // Parallel heat conduction

          if (dneut > 0.0) {
            ddt(Pn) +=
                Div_par_diffusion(Dn * Pn, logPn); // Perpendicular diffusion
          }
        }

        if ((hyper > 0.0) && (rhs_implicit)) {
          ddt(Pn) += D(Pn, hyper);
        }

        // Switch off evolution at very low densities
        // This seems to be necessary to get through initial transients

        for (auto i : ddt(Nn).getRegion(RGN_NOBNDRY)) {
          if (Nn[i] < 1e-5) {
            // Relax to the plasma temperature
            ddt(Pn)[i] = -1e-2 * (Pn[i] - ions.T[i] * Nn[i]);
          }
        }
      }

      if (rhs_explicit) {
        // Boundary condition on fluxes

        TRACE("Fluxes");

        BoutReal nredist;
        for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
          int jz = 0; // Z index
          int jy = mesh->yend;

          BoutReal ion_density = 0.5 * (ions.N(r.ind, jy, jz) + ions.N(r.ind, jy + 1, jz));
          BoutReal ion_velocity = 0.5 * (ions.V(r.ind, jy, jz) + ions.V(r.ind, jy + 1, jz));

          // Ion flux to the target
          flux_ion = ion_density * ion_velocity *
            (coord->J(r.ind, jy) + coord->J(r.ind, jy + 1)) /
            (sqrt(coord->g_22(r.ind, jy)) + sqrt(coord->g_22(r.ind, jy + 1)));
          
          BoutReal flux_neut = 0.0;

          for (int j = mesh->yend + 1; j < mesh->LocalNy; j++) {
            // flux_ion += ddt(Ne)(r.ind, j, jz) * coord->J(r.ind,j) *
            // coord->dy(r.ind,j);
            flux_neut += ddt(Nn)(r.ind, j, jz) * coord->J(r.ind, j) *
                         coord->dy(r.ind, j);

            //ddt(Ne)(r.ind, j, jz) = 0.0;
            ddt(Nn)(r.ind, j, jz) = 0.0;
          }

          // Make sure that mass is conserved

          // Total amount of neutral gas to be added
          BoutReal nadd = flux_ion * frecycle + flux_neut + gaspuff;

          // Neutral gas arriving at the target
          BoutReal ntarget =
              (1 - fredistribute) * nadd /
              (coord->J(r.ind, mesh->yend) * coord->dy(r.ind, mesh->yend));

          ddt(Nn)(r.ind, mesh->yend, jz) += ntarget;

          if (evolve_nvn) {
            // Set velocity of neutrals coming from the wall to a fraction of
            // the Franck-Condon energy
            BoutReal Vneut = -vwall * sqrt(3.5 / Tnorm);
            ddt(NVn)(r.ind, mesh->yend, jz) += ntarget * Vneut;
          }

          if (evolve_pn) {
            // Set temperature of the incoming neutrals to F-C
            ddt(Pn)(r.ind, mesh->yend, jz) += ntarget * (3.5 / Tnorm);
          }

          // Re-distribute neutrals
          nredist = fredistribute * nadd;

          // Divide flux_ion by J so that the result in the output file has
          // units of flux per m^2
          flux_ion /= coord->J(mesh->xstart, mesh->yend + 1);
        }

        // Now broadcast redistributed neutrals to other processors
        MPI_Comm ycomm = mesh->getYcomm(mesh->xstart); // MPI communicator
        int np;
        MPI_Comm_size(ycomm, &np); // Number of processors

        // Broadcast from final processor (presumably with target)
        // to all other processors
        MPI_Bcast(&nredist, 1, MPI_DOUBLE, np - 1, ycomm);

        // Distribute along length
        for (int j = mesh->ystart; j <= mesh->yend; j++) {
          // Neutrals into this cell
          // Note: from earlier normalisation the sum ( redist_weight * J * dy )
          // = 1 This ensures that if redist_weight is constant then the source
          // of particles per volume is also constant.
          BoutReal ncell = nredist * redist_weight(mesh->xstart, j);

          ddt(Nn)(mesh->xstart, j, 0) += ncell;

          // No momentum

          if (evolve_pn) {
            // Set temperature of the incoming neutrals to F-C
            ddt(Pn)(mesh->xstart, j, 0) += ncell * (3.5 / Tnorm);
          }
        }
      }
    }

    // Equal electron and ion temperatures
    ddt(ions.P)   += 0.5*ddt(electrons.P);
    
    // Add reaction sources
    if (atomic && rhs_explicit) {
      // Plasma equations sum electron and ion contributions
      try {
        ddt(Nn) += ddt(species.at("h")->N);
        if (evolve_nvn) {
          ddt(NVn) += ddt(species.at("h")->NV);
        }
        if (evolve_pn) {
          ddt(Pn) += ddt(species.at("h")->P);
        }
      } catch (const std::out_of_range &e) {
        throw BoutException("Failed while adding sources");
      }
    }
    
    return 0;
  }

  /*!
   * Preconditioner. Solves the heat conduction
   *
   * @param[in] t  The simulation time
   * @param[in] gamma   Factor in front of the Jacobian in (I - gamma*J).
   * Related to timestep
   * @param[in] delta   Not used here
   */
  int precon(BoutReal UNUSED(t), BoutReal gamma, BoutReal UNUSED(delta)) {
    
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

private:

  // Normalisation parameters
  BoutReal Tnorm, Nnorm, Bnorm;
  BoutReal Cs0, Omega_ci, rho_s0, tau_e0, mi_me;

  /////////////////////////////////////////////////////////////////
  // Evolving quantities
  Field3D Nn, NVn, Pn; // Neutral density, momentum, pressure

  Field3D Vn; // Ion and neutral velocities

  bool evolve_nvn; // Evolve neutral momentum?
  bool evolve_pn;  // Evolve neutral pressure?

  SpeciesMap species; // Map of species, indexed by strings
  
  /////////////////////////////////////////////////////////////////
  // Diffusion and viscosity coefficients

  Field3D Dn;     // Neutral gas diffusion
  BoutReal dneut; // Neutral gas diffusion multiplier

  Field3D kappa_n;    // Neutral gas thermal conduction
  Field3D kappa_epar; // Plasma thermal conduction

  Field3D tau_e;        // Electron collision time
  bool ion_viscosity;   // Braginskii ion viscosity on/off
  bool heat_conduction; // Thermal conduction on/off

  BoutReal nloss; // Neutral loss rate (1/timescale)

  /////////////////////////////////////////////////////////////////
  // Atomic physics transfer channels

  bool atomic; // Include atomic physics? This includes neutral gas evolution
  
  UpdatedRadiatedPower hydrogen; // Atomic rates

  BoutReal fimp;             // Impurity fraction (of Ne)
  bool impurity_adas;        // True if using Atomic++ library
  std::string impurity_species;   // Name of impurity species to use
  
  std::vector<Reaction*> reactions; // Reaction set to include
  
  ///////////////////////////////////////////////////////////////
  // Sheath boundary

  BoutReal sheath_gamma;  // Sheath heat transmission factor
  BoutReal neutral_gamma; // Neutral heat transmission

  int density_sheath;  // How to handle density boundary?
  int pressure_sheath; // How to handle pressure boundary?

  bool bndry_flux_fix;

  BoutReal frecycle; // Recycling fraction
  BoutReal gaspuff;  // Additional source of neutral gas at the target plate
  BoutReal vwall;    // Velocity of neutrals coming from the wall
                     // as fraction of Franck-Condon energy

  BoutReal flux_ion; // Flux of ions to target (output)

  // Re-distribution of recycled neutrals
  Field2D redist_weight;  // Weighting used to decide redistribution
  BoutReal fredistribute; // Fraction of recycled neutrals re-distributed along
                          // length

  ///////////////////////////////////////////////////////////////
  // Sources

  bool volume_source;         // Include volume sources?
  BoutReal powerflux;         // Used if no volume sources

  

  ///////////////////////////////////////////////////////////////
  // Numerical dissipation

  BoutReal tn_floor; // Minimum neutral gas temperature [eV]

  BoutReal hyper, viscos; // Numerical dissipation terms
  BoutReal ADpar;         // Added Dissipation numerical term

  Field2D dy4; // SQ(SQ(coord->dy)) cached to avoid recalculating

  BoutReal gamma_sound; // Ratio of specific heats in numerical dissipation term

  // Numerical diffusion
  const Field3D D(const Field3D &f, BoutReal d) {
    if (d < 0.0)
      return 0.0;
    return Div_par_diffusion(d * SQ(mesh->getCoordinates()->dy), f);
    // return -D4DY4_FV(d*dy4,f);
  }

  ///////////////////////////////////////////////////////////////
  // Splitting into implicit and explicit
  bool rhs_implicit, rhs_explicit; // Enable implicit and explicit parts
  bool update_coefficients;        // Re-calculate diffusion coefficients
};

BOUTMAIN(SD1D);
