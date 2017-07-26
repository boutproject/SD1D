#pragma once

class FluidFull;

#ifndef __SD1D_FLUID_FULL_H__
#define __SD1D_FLUID_FULL_H__

#include "fluid.hxx"

/// Full fluid model
///
/// Implements a 1D model which evolves density, pressure and momentum
/// 
class FluidFull : public Fluid {
public:
  ///
  /// @param[in] gamma_    Ratio of specific heats: 3/2 for monatomic gas
  /// @param[in] AA_       Atomic mass number : 2 for Deuterium
  FluidFull(BoutReal gamma_, BoutReal AA_) : gamma(gamma_), AA(AA_) {}
  
  void init(Solver *solver,           // Time integration solver
            bool restarting,          // True if simulation is restarting
            const string &suffix,     // Suffix for all evolving variables
            BoutReal Tnorm, BoutReal Nnorm, // Temperature, density normalisation [eV], [m^-3]
            BoutReal Cs0, BoutReal rho_s0) override;  // Velocity, length normalisation [m/s], [m]

  /// Calculate time derivatives
  void rhs(BoutReal time) override;

  BoutReal density(const DataIterator &i) override; // Return normalised density
  BoutReal temperature(const DataIterator &i) override; // Return normalised temperature
  BoutReal velocity(const DataIterator &i) override;  // normalised velocity

  void addParticles(const DataIterator &i, BoutReal Sn) override;
  void addEnergy(const DataIterator &i, BoutReal Se) override;
  void addMomentum(const DataIterator &i, BoutReal Sm) override;
private:
  BoutReal gamma; ///< Ratio of specific heats (adiabatic index, 5/3 for ideal monatomic gas)
  BoutReal AA;    ///< Atomic mass number (2 = Deuterium)
  
  Field3D N, P, NV; ///< density, pressure, parallel momentum
  
  Field3D T, V; // Derived temperature, parallel velocity
};

#endif //__SD1D_FLUID_FULL_H__
