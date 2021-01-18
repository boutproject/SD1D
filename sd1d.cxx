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
#include <bout/snb.hxx>
#include <derivs.hxx>
#include "field_factory.hxx"
#include <invert_parderiv.hxx>
#include "unused.hxx"

#include "div_ops.hxx"
#include "loadmetric.hxx"
#include "radiation.hxx"
#include "species.hxx"
#include "reaction.hxx"

using bout::HeatFluxSNB;

class SD1D : public PhysicsModel {
protected:
  int init(bool restarting) {
    auto &opt = Options::root()["sd1d"];

    // Normalisation
    OPTION(opt, Tnorm, 100);             // Reference temperature [eV]
    OPTION(opt, Nnorm, 1e19);            // Reference density [m^-3]
    OPTION(opt, Bnorm, 1.0);             // Reference magnetic field [T]
    SAVE_ONCE3(Tnorm, Nnorm, Bnorm);      // Save normalisations

    OPTION(opt, heat_conduction, true); // Spitzer-Hahm heat conduction

    snb_model = opt["snb_model"]
                    .doc("Use SNB non-local heat flux model")
                    .withDefault<bool>(false);	// SNB non-locality
    if (snb_model) {
      // Create a solver to calculate the SNB heat flux
      snb = new HeatFluxSNB();
    }


    // Get a list of reactions to include, separated by commas
    
    std::string reaction_names = opt["reactions"].withDefault<std::string>("");
    for (const auto& name : strsplit(reaction_names, ',')) {
      auto name_trimmed = trim(name);
      reactions.push_back(ReactionFactory::getInstance().create(
          name_trimmed, &Options::root()[name_trimmed]));
    }

    // Calculate normalisation factors
    // Note: These are calculated for Hydrogen, with
    // factors of atomic mass and charge included in the equations
    
    Cs0 = sqrt(SI::qe * Tnorm / SI::Mp); // Reference sound speed [m/s]
    Omega_ci = SI::qe * Bnorm / SI::Mp;  // Ion cyclotron frequency [1/s]
    rho_s0 = Cs0 / Omega_ci;             // Length scale [m]

    mi_me = SI::Mp / SI::Me;

    BoutReal Coulomb = 6.6 - 0.5 * log(Nnorm * 1e-20) + 1.5 * log(Tnorm);
    tau_e0 = 1. / (2.91e-6 * (Nnorm / 1e6) * Coulomb * pow(Tnorm, -3. / 2));

    // Save normalisation factors
    SAVE_ONCE5(Cs0, Omega_ci, rho_s0, tau_e0, mi_me);
    
    // Electrons are handled differently,
    // so not an evolving species
    species["e"] = new Species(SI::Me / SI::Mp, -1);

    // Save some electron properties
    auto &electrons = *species["e"];
    dump.addRepeat(electrons.N, "Ne");
    dump.addRepeat(electrons.T, "Te");

    // A list of species to evolve
    std::string species_names = opt["species"].withDefault<std::string>("h+");
    for (const auto& name : strsplit(species_names, ',')) {
      auto name_trimmed = trim(name);

      species[name_trimmed] = new FluidSpecies(name_trimmed,
                                               &Options::root()[name_trimmed],
                                               solver,
                                               restart,
                                               restarting,
                                               Nnorm, Tnorm, Omega_ci, Cs0); 
    }

    // Load the metric tensor
    LoadMetric(rho_s0, Bnorm);
    
    // Read the flux-tube area from input file
    // This goes into the Jacobian.
    std::string area_string = opt["area"].withDefault<std::string>("1.0");
    
    mesh->getCoordinates()->J = FieldFactory(mesh).create2D(area_string, Options::getRoot());
    
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
                            ReactionFactory::getInstance().create("atomic++coronal", &opt));
      } else {
        // Use carbon radiation for the impurity
        impurity_species = "c";
        reactions.push_back(
                            ReactionFactory::getInstance().create("c_hutchinson", &opt));
      }

      // Create a Species object for the impurity
      species[impurity_species] = new Species();
    }

    if (heat_conduction) {
      SAVE_REPEAT(kappa_epar); // Save coefficient of thermal conduction
    }

    if (snb_model) {
      SAVE_REPEAT(Div_Q_SH, Div_Q_SNB);
    }
    
    kappa_epar = 0.0;
    
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

    ///////////////////////////////////////////////////////////
    // Evolve all ion and neutral species
    
    for(auto &s : species) {
      s.second->evolve(time);
    }

    //////////////////////////////////////////////////////////
    // Electrons
    
    auto &electrons = *species.at("e");
    
    // Set electron species properties
    auto &ions = *species.at("h+");
    electrons.T = ions.T;
    electrons.N = ions.N * ions.ZZ;
    electrons.P = ions.P * ions.ZZ;
    electrons.V = ions.V;

    //////////////////////////////////////////////////////////
    // Fixed fraction inpurity species
    
    if (fimp > 0.0) {
      // Fixed fraction impurity radiation 
      species.at(impurity_species)->N = fimp * ions.N;
    }

    //////////////////////////////////////////////////////////
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

    ///////////////////////////////////////////////////
    // electron model
    
    
    {
      TRACE("Electron pressure");
      
      ddt(electrons.P) = 0.0;

      if (heat_conduction) {

        if (update_coefficients) {
          // Update diffusion coefficients
          TRACE("Update coefficients");
          
          tau_e = Omega_ci * tau_e0 * pow(electrons.T, 1.5) / electrons.N;
          
          kappa_epar = 3.2 * mi_me * 0.5 * electrons.P * tau_e;
          kappa_epar.applyBoundary("neumann");
          
        }

        if (snb_model) {
            // SNB non-local heat flux. Also returns the Spitzer-Harm value for comparison
            // Note: Te in eV, Ne in Nnorm
            Field2D dy_orig = mesh->getCoordinates()->dy;
            mesh->getCoordinates()->dy *= rho_s0; // Convert distances to m
            Div_Q_SNB = snb->divHeatFlux(electrons.T * Tnorm, electrons.N * Nnorm, &Div_Q_SH);
            mesh->getCoordinates()->dy = dy_orig;
            
            // Normalise from eV/m^3/s
            Div_Q_SNB /= Tnorm * Nnorm * Omega_ci;
            Div_Q_SH /= Tnorm * Nnorm * Omega_ci;

            // Add to pressure equation
            ddt(electrons.P) -= (2. / 3) * Div_Q_SNB * ions.ZZ / (1+ions.ZZ);
          } else {
            // The standard Spitzer-Harm model
            // NOTE: Factor of Z/(1+Z) accounts for instant equipartition of energy between electrons and ions in this Te=Ti assumption,

            ddt(electrons.P) += (2. / 3) * Div_par_diffusion_upwind(kappa_epar, electrons.T) * ions.ZZ / (1+ions.ZZ);
        }
      }

      // Electron pressure acts on ions
      ddt(ions.NV) -= Grad_par(ions.P);
    }
    
    // Equal electron and ion temperatures
    ddt(ions.P)   += ddt(electrons.P) / ions.ZZ;
    
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
    static InvertPar *inv = NULL;
    if (!inv) {
      // Initialise parallel inversion class
      inv = InvertPar::Create();
      inv->setCoefA(1.0);
    }
    if (heat_conduction) {
      auto &ions = *species.at("h+");
      
      // Set the coefficient in front of Grad2_par2
      inv->setCoefB(-(2. / 3) * gamma * kappa_epar);
      Field3D dT = ddt(ions.P);
      dT.applyBoundary("neumann");
      ddt(ions.P) = inv->solve(dT);
    }

    // if (atomic) {
    //   if (evolve_pn && (dneut > 0.0)) {
    //     // Neutral pressure
    //     inv->setCoefB(-(2. / 3) * gamma * kappa_n);
    //     Field3D dT = ddt(Pn);
    //     dT.applyBoundary("neumann");
    //     ddt(Pn) = inv->solve(dT);
    //   }

    //   if (dneut > 0.0) {
    //     inv->setCoefB(-gamma * Dn);
    //     Field3D tmp = ddt(Nn);
    //     tmp.applyBoundary("neumann");
    //     ddt(Nn) = inv->solve(tmp);
    //   }
    // }
    
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

  SpeciesMap species; // Map of species, indexed by strings
  
  /////////////////////////////////////////////////////////////////
  // Electron model

  Field3D kappa_epar; // Plasma thermal conduction
  Field3D tau_e;        // Electron collision time
  bool heat_conduction; // Thermal conduction on/off
  bool snb_model;       // Use the SNB model for heat conduction?
  HeatFluxSNB *snb;
  Field3D Div_Q_SH, Div_Q_SNB; // Divergence of heat flux from Spitzer-Harm and SNB

  /////////////////////////////////////////////////////////////////
  // Atomic physics transfer channels

  BoutReal fimp;             // Impurity fraction (of Ne)
  bool impurity_adas;        // True if using Atomic++ library
  std::string impurity_species;   // Name of impurity species to use
  
  std::vector<Reaction*> reactions; // Reaction set to include
    
  ///////////////////////////////////////////////////////////////
  // Splitting into implicit and explicit
  bool rhs_implicit, rhs_explicit; // Enable implicit and explicit parts
  bool update_coefficients;        // Re-calculate diffusion coefficients
};

BOUTMAIN(SD1D);
