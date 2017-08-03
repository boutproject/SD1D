#!/usr/bin/env python
#
#
# Reads the final time slice from the given path
# and plots transfer channels of particles, momentum and energy
#
# NOTE: This assumes that SD1D was run with diagnose=true
#       so that the different channels are saved

import sys


import matplotlib.pyplot as plt

from boutdata import collect
from numpy import zeros

# Check command-line arguments
if len(sys.argv) != 2:
    # Print usage information
    print("Usage: {0} path\n e.g. {0} data".format(sys.argv[0]))
    sys.exit(1)
    
# First (and only) argument is the path
path = sys.argv[1]

t = collect("t_array", path=path)
tind = len(t)-1 # Get the last time point

var_list = ["Ne", "P",       # Plasma profiles
            "Srec", "Siz",   # Particle sources / sinks
            "Frec", "Fiz", "Fcx", "Fel",   # Momentum source / sinks to neutrals
            "Rrec", "Riz", "Rzrad", "Rex", # Radiation, energy loss from system
            "Erec", "Eiz", "Ecx", "Eel",   # Energy transfer between neutrals and plasma
            "Nnorm", "Tnorm", "Omega_ci", "rho_s0", "Cs0"]  # Normalisations

########################################################
# Position
dy = collect("dy", path=path)[0,:]
n = len(dy)
pos = zeros(n)

# position at the centre of the grid cell
pos[0] = 0.5*dy[0]
for i in range(1,n):
    pos[i] = pos[i-1] + 0.5*dy[i-1] + 0.5*dy[i]
    
########################################################
# Read the data into a dictionary

data = {}
for var in var_list:
    try:
        data[var] = collect(var, tind=tind, path=path)
        
        if len(data[var].shape) == 4:
            # 4D variable
            data[var] = data[var][0,0,:,0] # Make 1D [y]
    except:
        print("Variable '%s' not found" % (var,))
        data[var] = None


        
########################################################
# Particle sources

Snorm = data["Nnorm"]*data["Omega_ci"] # Normalisation factor

plt.plot(pos, data["Srec"]*Snorm, label="Recombination (Srec)")
plt.plot(pos, data["Siz"]*Snorm, label="Ionisation (Siz)")
plt.xlabel("Parallel location [m]")
plt.ylabel(r"Plasma particle loss rate [m$^{-3}$s$^{-1}$]")
plt.legend(loc="upper left")

plt.savefig("processes_part.pdf")

plt.show()

########################################################
# Momentum sources

Fnorm = data["Nnorm"]*data["Cs0"]

plt.plot(pos, data["Frec"]*Fnorm, label="Recombination (Frec)")
plt.plot(pos, data["Fiz"]*Fnorm, label="Ionisation (Fiz)")
plt.plot(pos, data["Fcx"]*Fnorm, label="Charge exchange (Fcx)")
if data["Fel"] is not None:
    plt.plot(pos, data["Fel"]*Enorm, label="Elastic scattering (Fel)")
    
plt.xlabel("Parallel location [m]")
plt.ylabel(r"Plasma momentum transfer to neutrals [m$^{-2}$s$^{-1}$]")
plt.legend(loc="upper left")

plt.savefig("processes_mom.pdf")

plt.show()

########################################################
# Energy losses

Enorm = 1.602e-19*data["Tnorm"]*data["Nnorm"]*data["Omega_ci"]

plt.plot(pos, (data["Rrec"]+data["Erec"])*Enorm, label="Recombination (Rrec+Erec)")
plt.plot(pos, (data["Riz"]+data["Eiz"])*Enorm, label="Ionisation (Riz+Eiz)")
plt.plot(pos, data["Ecx"]*Enorm, label="Charge exchange (Ecx)")
plt.plot(pos, data["Rzrad"]*Enorm, label="Impurity radiation (Rzrad)")
if data["Rex"] is not None:
    plt.plot(pos, data["Rex"]*Enorm, label="Hydrogen excitation (Rex)")
if data["Eel"] is not None:
    plt.plot(pos, data["Eel"]*Enorm, label="Elastic scattering (Eel)")
    
plt.xlabel("Parallel location [m]")
plt.ylabel(r"Plasma energy loss rate [W m$^{-3}$]")
plt.legend(loc="upper left")

plt.xlim([29,30])

plt.savefig("processes_energy.pdf")

plt.show()


########################################################
# All plots
import matplotlib
matplotlib.rcParams.update({'font.size': 6})

f, axarr = plt.subplots(4, sharex=True)
axarr[0].plot(pos, data["Ne"]*data["Nnorm"], color='b')
axarr[0].set_ylabel(r"Electron density [m$^{-3}$]", color='b')
axarr[0].tick_params('y', colors='b')

ax2 = axarr[0].twinx()
Te = 0.5*data["P"]/data["Ne"]
ax2.plot(pos, Te*data["Tnorm"], color='r')
ax2.set_ylabel("Electron temperature [eV]", color='r')
ax2.tick_params('y', colors='r')

axarr[1].plot(pos, data["Srec"]*Snorm, label="Recombination (Srec)")
axarr[1].plot(pos, data["Siz"]*Snorm, label="Ionisation (Siz)")
axarr[1].set_title(r"Plasma particle loss rate [m$^{-3}$s$^{-1}$]")
axarr[1].legend(loc="upper left")

axarr[2].plot(pos, data["Frec"]*Fnorm, label="Recombination (Frec)")
axarr[2].plot(pos, data["Fiz"]*Fnorm, label="Ionisation (Fiz)")
axarr[2].plot(pos, data["Fcx"]*Fnorm, label="Charge exchange (Fcx)")
if data["Fel"] is not None:
    axarr[2].plot(pos, data["Fel"]*Enorm, label="Elastic scattering (Fel)")

axarr[2].set_title(r"Plasma momentum transfer to neutrals [m$^{-2}$s$^{-1}$]")
axarr[2].legend(loc="upper left")

axarr[3].plot(pos, (data["Rrec"]+data["Erec"])*Enorm, label="Recombination (Rrec+Erec)")
axarr[3].plot(pos, (data["Riz"]+data["Eiz"])*Enorm, label="Ionisation (Riz+Eiz)")
axarr[3].plot(pos, data["Ecx"]*Enorm, label="Charge exchange (Ecx)")
axarr[3].plot(pos, data["Rzrad"]*Enorm, label="Impurity radiation (Rzrad)")

if data["Rex"] is not None:
    axarr[3].plot(pos, data["Rex"]*Enorm, label="Hydrogen excitation (Rex)")
if data["Eel"] is not None:
    axarr[3].plot(pos, data["Eel"]*Enorm, label="Elastic scattering (Eel)")
    
axarr[3].set_title(r"Plasma energy loss rate [W m$^{-3}$]")
axarr[3].set_xlabel("Parallel location [m]")
axarr[3].legend(loc="upper left")

axarr[3].set_xlim([29,30])

plt.savefig("processes_all.pdf")

plt.show()
