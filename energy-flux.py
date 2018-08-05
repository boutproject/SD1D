#!/usr/bin/env python
#
# Evaluate energy sources, sinks, and fluxes

from boutdata import collect
from numpy import zeros, sum

import os

def heat_conduction(pos, Te, kappa0=2293.8117):
    """
    Calculate heat conduction in W/m^2 
    given input 1-D profiles
    
    pos[y] position in meters
    Te[y]  Electron temperature [eV]

    kappa0 Coefficient of heat conduction, so kappa = kappa0 * Te^5/2    

    Note: The value of kappa0 is what is used in SD1D
    """
    grad_Te = (Te[1:] - Te[:-1]) / (pos[1:] - pos[:-1]) # One-sided differencing

    Te_p = 0.5*(Te[1:] + Te[:-1])
    
    # Position of the result
    result_pos = 0.5*(pos[1:] + pos[:-1])
    
    return result_pos, -2293.8117*Te_p**(5./2)*grad_Te

def ke_convection(pos, n, vi, AA=2):
    """
    Calculate kinetic energy convection in W/m^2 
    given input 1-D profiles
    
    pos[y] position in meters
    n[y]   Density in m^-3
    vi[y]  Parallel flow velocity in m/s
    
    AA  Atomic mass number
    """

    # Interpolate onto cell boundaries
    vi_p = 0.5*(vi[1:] + vi[:-1])
    n_p = 0.5*(n[1:] + n[:-1])

    # Position of the result
    result_pos = 0.5*(pos[1:] + pos[:-1])

    return result_pos, 0.5*n_p*vi_p**3 * AA * 1.67e-27

def thermal_convection(pos, p, vi):
    """
    Calculate thermal energy convection in W/m^2 
    given input 1-D profiles
    
    pos[y] position in meters
    n[y]   Density in m^-3
    vi[y]  Parallel flow velocity in m/s
    """

    # Interpolate onto cell boundaries
    vi_p = 0.5*(vi[1:] + vi[:-1])
    p_p = 0.5*(p[1:] + p[:-1])

    # Position of the result
    result_pos = 0.5*(pos[1:] + pos[:-1])

    return result_pos, (5./2)*p_p*vi_p

def energy_flux(path, tind=-1):
    """
    Calculates the energy flux due to conduction and convection
    
    path  Path to the data files
    tind  Time index. By default the final time
    """
    
    # Evolving variables, remove extra guard cells so just one each side
    P = collect("P", path=path, tind=tind, yguards=True)[-1,0,1:-1,0]
    Ne = collect("Ne", path=path, tind=tind, yguards=True)[-1,0,1:-1,0]
    NVi = collect("NVi", path=path, tind=tind, yguards=True)[-1,0,1:-1,0]

    # Normalisations
    nnorm = collect("Nnorm", path=path, tind=tind)
    tnorm = collect("Tnorm", path=path, tind=tind)
    pnorm = nnorm*tnorm*1.602e-19 # Converts p to Pascals
    cs0 = collect("Cs0", path=path)

    try:
        kappa_epar = collect("kappa_epar", path=path, tind=tind)
    except:
        kappa_epar = None
    
    # electron temperature
    Te = (0.5*P/Ne) * tnorm

    # ion parallel velocity
    Vi = (NVi/Ne) * cs0

    NVi *= nnorm * cs0
    Ne *= nnorm
    P *= pnorm
    
    # Source
    pesource = collect("PeSource", path=path, yguards=True)
    
    dy = collect("dy", path=path, yguards=True)[0,1:-1]
    n = len(dy)
    pos = zeros(n)

    # position at the centre of the grid cell
    pos[0] = -0.5*dy[1]
    pos[1] = 0.5*dy[1]
    for i in range(2,n):
        pos[i] = pos[i-1] + 0.5*dy[i-1] + 0.5*dy[i]
    
    # Calculate energy transport
    flux_pos, conduction = heat_conduction(pos, Te)
    _, convect_ke = ke_convection(pos, Ne, Vi)
    _, convect_therm = thermal_convection(pos, P, Vi)

    if kappa_epar is None:
        conduction = zeros(len(flux_pos))
    
    # Return results as a dictionary
    return {"position":flux_pos,
            "conduction":conduction,
            "convection_ke":convect_ke,
            "convection_thermal":convect_therm}

    
if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="Directory containing the BOUT++ data")
    parser.add_argument("-s", "--fontsize", type=int, default=14, help="Font size, default 14")
    parser.add_argument('-x', "--nox", action='store_true', help="Don't use X display. Uses Agg backend.")
    
    args = parser.parse_args()
    
    import matplotlib
    matplotlib.rcParams.update({'font.size': args.fontsize})
    if args.nox:
        matplotlib.use('Agg')

    path = args.path
        
    import matplotlib.pyplot as plt
        
    result = energy_flux(path)

    pos = result["position"]

    plt.plot(pos, result["conduction"], label="Conduction")
    
    plt.plot(pos, result["convection_ke"], label="Kinetic convection")
    plt.plot(pos, result["convection_thermal"], label="Thermal convection")
    plt.plot(pos, result["convection_ke"] + result["convection_thermal"], label="Total convection")

    plt.plot(pos, result["convection_ke"] + result["convection_thermal"] + result["conduction"], label="Total energy flux")
    plt.xlabel("Position [m]")
    plt.ylabel(r"Heat flux [W/m$^2$]")
    plt.legend(loc="best")

    plt.savefig(os.path.join(path, "energy-flux.pdf"))
    plt.savefig(os.path.join(path, "energy-flux.png"))
    
    print("Saved figure as " + os.path.join(path, "energy-flux.pdf"))

    if args.nox:
        sys.exit(0)
    
    plt.show()
