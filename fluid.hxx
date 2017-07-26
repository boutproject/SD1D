
#pragma once

#ifndef __SD1D_FLUID_H__
#define __SD1D_FLUID_H__

#include <bout.hxx>


/// Fluid class
///
/// Defines the interface for fluid species
class Fluid {
public:
  Fluid() {}
  virtual ~Fluid() {}

  virtual void init(Solver *solver,           // Time integration solver
                    bool restarting,          // True if simulation is restarting
                    const string &suffix,     // Suffix for all evolving variables
                    BoutReal Tnorm, BoutReal Nnorm, // Temperature, density normalisation [eV], [m^-3]
                    BoutReal Cs0, BoutReal rho_s0) = 0;  // Velocity, length normalisation [m/s], [m]

  /// Calculate time derivatives
  virtual void rhs(BoutReal time) = 0;

  virtual BoutReal density(const DataIterator &i) = 0; // Return normalised density
  virtual BoutReal temperature(const DataIterator &i) = 0; // Return normalised temperature
  virtual BoutReal velocity(const DataIterator &i) = 0;  // normalised velocity

  /// Add a particle source. 
  /// 
  /// ddt(n) = ... + Sn
  ///
  /// @param[in] Sn  normalised particle source rate
  virtual void addParticles(const DataIterator &i, BoutReal Sn) = 0;
  
  /// Add an energy source (negative for energy loss)
  ///
  /// ddt(P) = ... + (2/3)*Se
  ///
  /// Note the factor of 2/3 to convert from power to change in pressure
  /// for monatomic gas
  ///
  /// @param[in] Se  Normalised power source [eV/m^-3/s]
  virtual void addEnergy(const DataIterator &i, BoutReal Se) = 0;

  /// Add a parallel momentum source
  ///
  /// ddt(NVi) = ... + Sm
  ///
  /// @param[in] Sm   Normalised momentum source [kg m/s /m^3]
  virtual void addMomentum(const DataIterator &i, BoutReal Sm) = 0;
  
protected:
  
};

#endif // __SD1D_FLUID_H__
  
