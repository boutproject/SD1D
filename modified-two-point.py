#!/usr/bin/env python

# Modified two-point model
# 

from boutdata import collect
from boutdata.data import BoutData 
import numpy as np

def momentum_loss_fraction(path, upstream_index=0, tind=-1):
    """
    The fraction of momentum lost between an upstream
    index and the end of the domain
    
    f_mom = \int F dl / p_u

    where F is the friction, dl the cell length, and 
    p_u the total (static+dynamic) pressure at upstream_index

    Inputs
    ------

    path  The directory containing the data
    upstream_index   The index to use as the upstream 
    tind   The time index. -1 is the final time
    
    Return
    ------
    
    Fraction of momentum lost
    
    """
    
    # Read the upstream static pressure 
    p = collect("P", path=path, yind=upstream_index, tind=tind)
    p_u_static = np.asscalar(p) # Convert to scalar
    
    # Calculate the upstream dynamic pressure
    n = collect("Ne", path=path, yind=upstream_index, tind=tind)
    nvi = collect("NVi", path=path, yind=upstream_index, tind=tind)
    p_u_dynamic = np.asscalar( nvi*nvi/n )

    # Total pressure p_u
    p_u = p_u_static + p_u_dynamic
    
    # Read the friction
    try:
        F = collect("F", path=path, yind=[upstream_index,-1], tind=tind).flatten()
    except:
        print("No friction force 'F' found")
        F = 0.0
    
    # Need to integrate along length, so get cell spacing
    dy = collect("dy", path=path, yind=[upstream_index,-1]).flatten()
    g_22 = collect("g_22", path=path, yind=[upstream_index,-1]).flatten()
    
    dl = dy * np.sqrt(g_22) # Cell length
    
    # F is already integrated over each cell, so just sum
    # contribution from each cell
    return np.sum( F * dl ) / p_u
    
    
def power_loss_fraction(path, upstream_index=0, tind=-1):
    """
    The fraction of internal energy lost between an upstream
    index and the end of the domain.
    
    This takes the input power from BOUT.inp as the upstream
    power flux.

    
    """
    
    # Read the energy exchange with neutrals and radiation
    try:
        E = collect("E", path=path, yind=[upstream_index,-1], tind=tind).flatten()
    except:
        print("Energy loss rate 'E' not found")
        E = 0.0
    try:
        R = collect("R", path=path, yind=[upstream_index,-1], tind=tind).flatten()
    except:
        print("Radiation rate 'R' not found")
        R = 0.0
    
    # Compression term P * Div_par(Vi)
    n = collect("Ne", path=path, yguards=True, tind=tind).flatten()
    nvi = collect("NVi", path=path, yguards=True, tind=tind).flatten()
    vi = nvi / n
    J = collect("J", path=path, yguards=True).flatten()
    dy = collect("dy", path=path, yguards=True).flatten()
    g_22 = collect("g_22", path=path, yguards=True).flatten()
    
    # Parallel divergence
    Jv = J*vi
    div_par_v = (Jv[2:] - Jv[:-2])/(2.*dy[1:-1]*np.sqrt(g_22[1:-1])*J[1:-1])
    div_par_v = div_par_v[1:-1] # Remove guard cells
    div_par_v = div_par_v[upstream_index:]
    
    p = collect("P", path=path, yind=[upstream_index,-1], tind=tind).flatten()
    compression = p*div_par_v

    # Need to integrate along length, so get cell spacing
    dy = collect("dy", path=path, yind=[upstream_index,-1]).flatten()
    g_22 = collect("g_22", path=path, yind=[upstream_index,-1]).flatten()
    # Jacobian to calculate cell volume
    J = collect("J", path=path, yind=[upstream_index,-1]).flatten()
    # Get the energy flux
    data = BoutData(path)
    power_flux = data["options"]["p"]["powerflux"] # in W/m^2
    
    # Normalisations
    nnorm = collect("Nnorm", path=path)
    tnorm = collect("Tnorm", path=path)
    pnorm = nnorm*tnorm*1.602e-19 # Converts p to Pascals
    cs0 = collect("Cs0", path=path)
    
    # Normalise, converting from SI to SD1D units
    power_flux_normalised = power_flux / (pnorm * cs0)
    
    return np.sum( (E+R+compression) * J * dy ) / (J[0] * power_flux_normalised / np.sqrt(g_22[0]))

if __name__ == "__main__":
    import argparse
    
    import sys
    import matplotlib
    matplotlib.rcParams.update({'font.size': 20})
    import matplotlib.pyplot as plt
    
    # Create argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="Directory containing the BOUT++ data")
    parser.add_argument("-t", "--tind", type=int, default=-1, help="Time index. Default is -1 (last)")
    parser.add_argument("-u", "--upstream", type=int, default=None, help="Y index of upstream point. Default is the end of PeSource")
    args = parser.parse_args()
    
    path = args.path  # Path to the data
    tind = args.tind  # Time index
    
    upstream_index = args.upstream # Upstream index, at the end of the source region
    if upstream_index is None:
        try:
            # Locate the end of the source
            pesource = collect("PeSource", path=path).flatten()
            upstream_index = np.where(pesource < 1e-10)[0][0]  # Return is a tuple, take first element
            print("Setting upstream index = %d" % upstream_index)
        except:
            print("Error: Could not locate end of the source region")
            print(" -> Setting upstream index to 0")
            upstream_index = 0
    
    # Calculate loss fractions
    f_mom = momentum_loss_fraction(path, upstream_index=upstream_index, tind=tind)
    f_pwr = power_loss_fraction(path, upstream_index=upstream_index, tind=tind)
    
    # Read the upstream static pressure 
    p = collect("P", path=path, yind=upstream_index, tind=tind)
    p_u_static = np.asscalar(p) # Convert to scalar
    
    # Calculate the upstream dynamic pressure
    n = collect("Ne", path=path, yind=upstream_index, tind=tind)
    nvi = collect("NVi", path=path, yind=upstream_index, tind=tind)
    p_u_dynamic = np.asscalar( nvi*nvi/n )

    # Total pressure p_u
    p_u = p_u_static + p_u_dynamic

    # Normalisations
    nnorm = collect("Nnorm", path=path)
    tnorm = collect("Tnorm", path=path)
    pnorm = nnorm*tnorm*1.602e-19 # Converts p to Pascals
    cs0 = collect("Cs0", path=path)
    
    # Get target density, temperature
    n = collect("Ne", path=path, yguards=True, tind=tind).flatten()
    n_target = 0.5*(n[-2] + n[-3]) * nnorm
    p = collect("P", path=path, yguards=True, tind=tind).flatten()
    te = 0.5*p/n * tnorm
    te_target = 0.5*(te[-2] + te[-3])

    # Jacobian to calculate cell areas
    J = collect("J", path=path, yind=[upstream_index,-1]).flatten()
    
    # Area expansion factor
    f_R = J[-1] / J[0]  # B_u / B_t = A_t / A_u
    
    print("\n\n------------------------------------\n")
    print("Integrated momentum loss fraction: %e" % f_mom)
    print("Integrated internal energy loss fraction: %e [WRONG]" % f_pwr)

    print("Area expansion factor = %e" % f_R)

    p_u *= pnorm # Convert to SI [Pascals]

    print("Upstream pressure = %e Pa" % p_u)

    # Get the energy flux from input options
    data = BoutData(path)
    power_flux = data["options"]["p"]["powerflux"] # in W/m^2
    
    print("Upstream power flux = %e W/m^2" % power_flux)
    
    # Sheath heat transmission coefficient
    sheath_gamma = data["options"]["sd1d"]["sheath_gamma"]
    
    print("Sheath heat transmission gamma = %e" % sheath_gamma)
    
    m_i = 2.*1.67e-27 # Ion mass [kg]

    # Power into the sheath
    sheath_power = sheath_gamma * 1.602e-19*n_target*te_target*np.sqrt(2.*1.602e-19*te_target/m_i)
    print("Sheath power flux = %e W/m^2" % sheath_power)

    f_pwr = 1. - sheath_power/power_flux
    print("Total power loss fraction = %e" % f_pwr)
    
    # A fraction 1/gamma of the power into the sheath is as kinetic energy
    
    f_pwr = 1. - sheath_power*(sheath_gamma-1.0)/sheath_gamma/power_flux
    print("Internal energy loss fraction = %e" % f_pwr)

    #################################################
    # Density

    print("\n------------------------------------\n\nSD1D density at target = %e" % n_target)
    
    m2pm_density_1 = sheath_gamma**2 / (32.*m_i)
    m2pm_density_2 = p_u**3 / power_flux**2
    m2pm_density_3 = (1. - f_mom)**3/(1. - f_pwr)**2
    m2pm_density_4 = 1. # 4/(1+Ti/Te)^2
    m2pm_density_5 = 1. # 8*M/(1+M^2)^3
    m2pm_density_6 = f_R**2

    m2pm_density = m2pm_density_1 * m2pm_density_2 * m2pm_density_3 * m2pm_density_4 * m2pm_density_5 * m2pm_density_6
    
    print("Modified 2-point model density = %e" % m2pm_density)
    print("        [gamma^2/(32m)] * [p_u^3/q_u^2] * [(1-f_mom)^3/(1-f_pwr)^2] * [4/(1+Ti/Te)^2] * [8M^2/(1+M^2)^3] * [(B_u/B_t)^2]")
    
    print("          %e  *  %e *        %e       *  %e   *   %e   *  %e" % (
        m2pm_density_1, m2pm_density_2, m2pm_density_3, m2pm_density_4, m2pm_density_5, m2pm_density_6))
    
    #################################################
    # Temperature
    print("\n------------------------------------\n\nSD1D temperature at target = %e" % te_target)
    
    m2pm_te_1 = 8.*m_i/(1.602e-19*sheath_gamma**2)
    m2pm_te_2 = (power_flux / p_u)**2
    m2pm_te_3 = ((1. - f_pwr) / (1. - f_mom))**2
    m2pm_te_4 = 1. # (1+Ti/Te)/2
    m2pm_te_5 = 1. # (1+M^2)^2/4M^2
    m2pm_te_6 = 1./f_R**2
    
    m2pm_te = m2pm_te_1 * m2pm_te_2 * m2pm_te_3 * m2pm_te_4 * m2pm_te_5 * m2pm_te_6
    print("Modified 2-point model temperature = %e" % m2pm_te)
    
    print("        [8m/e*gamma^2] * [q_u^2/p_u^2] * [(1-f_pwr)^2/(1-f_mom)^2] * [(1+Ti/Te)/2] * [(1+M^2)^2/(4M^2)] * [(B_t/B_u)^2]")
    print("         %e  *  %e *        %e       *  %e   *   %e   *  %e" % (
        m2pm_te_1, m2pm_te_2, m2pm_te_3, m2pm_te_4, m2pm_te_5, m2pm_te_6))
    
