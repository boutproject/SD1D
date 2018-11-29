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
import numpy as np

# Check command-line arguments
if len(sys.argv) < 2:
    # Print usage information
    print("Usage: {0} path\n e.g. {0} data [zoom]".format(sys.argv[0]))
    sys.exit(1)
    
# First argument is the path
path = sys.argv[1]

if len(sys.argv) > 2:
    # Optional second argument is length from target
    zoomlength = float(sys.argv[2])
else:
    zoomlength = 1.0 # 1m default

t = collect("t_array", path=path)
J = collect("J", path=path)[0,:]
dy = collect("dy", path=path)[0,:]
dV = J*dy # Volume element

tind = len(t)-1 # Get the last time point

var_list = ["PeSource",      # Input pressure source
            "Ne", "P",       # Plasma profiles
            "Srec", "Siz",   # Particle sources / sinks
            "Frec", "Fiz", "Fcx", "Fel",   # Momentum source / sinks to neutrals
            "Rrec", "Riz", "Rzrad", "Rex", # Radiation, energy loss from system
            "E", "Erec", "Eiz", "Ecx", "Eel",   # Energy transfer between neutrals and plasma
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

try:
    plt.plot(pos, data["Srec"]*Snorm, label="Recombination (Srec)")
    plt.plot(pos, data["Siz"]*Snorm, label="Ionisation (Siz)")
except:
    pass
plt.xlabel("Parallel location [m]")
plt.ylabel(r"Plasma particle loss rate [m$^{-3}$s$^{-1}$]")
plt.legend(loc="upper left")

plt.savefig(path+"/processes_part.pdf")

plt.show()

########################################################
# Momentum sources

Fnorm = data["Nnorm"]*data["Cs0"]

try:
    plt.plot(pos, data["Frec"]*Fnorm, label="Recombination (Frec)")
    plt.plot(pos, data["Fiz"]*Fnorm, label="Ionisation (Fiz)")
    plt.plot(pos, data["Fcx"]*Fnorm, label="Charge exchange (Fcx)")
    if data["Fel"] is not None:
        plt.plot(pos, data["Fel"]*Enorm, label="Elastic scattering (Fel)")
except:
    pass
plt.xlabel("Parallel location [m]")
plt.ylabel(r"Plasma momentum transfer to neutrals [m$^{-2}$s$^{-1}$]")
plt.legend(loc="upper left")

plt.savefig(path+"/processes_mom.pdf")

plt.show()

########################################################
# Energy losses

input_power = (3./2)*np.sum(data["PeSource"] * dV)
if data["Rzrad"] is not None:
    impurity_loss = np.sum(data["Rzrad"] * dV)
else:
    impurity_loss = 0.0

if data["Riz"] is not None:
    ionisation_loss = np.sum(data["Riz"] * dV)
else:
    ionisation_loss = 0.0

if data["Rrec"] is not None:
    recombination_loss = np.sum(data["Rrec"] * dV) # Note: Negative
else:
    recombination_loss = 0.0

if data["E"] is not None:
    neut_loss = np.sum(data["E"] * dV)  # Note: can be negative
else:
    neut_loss = 0.0

if data["Rex"] is not None:
    excitation_loss = np.sum(data["Rex"] * dV)
    ionisation_loss += excitation_loss

if data["Siz"] is not None:
    ion_source = np.sum(data["Siz"] * dV)
    print("Effective Eionise = {0}".format(-data["Tnorm"] * ionisation_loss / ion_source))
    
iz_neut_loss = ionisation_loss + neut_loss

print("Input power: {0}".format(input_power))
print("Impurity loss: {0} ({1} %)".format(impurity_loss, 100.*impurity_loss/input_power))
print("Ionisation loss: {0} ({1} %)".format(ionisation_loss, 100.*ionisation_loss/input_power))
print("Ionisation loss + neutral exchange: {0} ({1} %)".format(iz_neut_loss, 100.*iz_neut_loss/input_power))
print("Ionisation / (Input - Impurity - Recombination - Neutrals): {0} %".format(100.*ionisation_loss / (input_power - impurity_loss - recombination_loss - neut_loss)))



Enorm = 1.602e-19*data["Tnorm"]*data["Nnorm"]*data["Omega_ci"]

try:
    plt.plot(pos, (data["Rrec"]+data["Erec"])*Enorm, label="Recombination (Rrec+Erec)")
    plt.plot(pos, (data["Riz"]+data["Eiz"])*Enorm, label="Ionisation (Riz+Eiz)")
    plt.plot(pos, data["Ecx"]*Enorm, label="Charge exchange (Ecx)")
    plt.plot(pos, data["Rzrad"]*Enorm, label="Impurity radiation (Rzrad)")
    if data["Rex"] is not None:
        plt.plot(pos, data["Rex"]*Enorm, label="Hydrogen excitation (Rex)")
    if data["Eel"] is not None:
        plt.plot(pos, data["Eel"]*Enorm, label="Elastic scattering (Eel)")
except:
    pass

plt.xlabel("Parallel location [m]")
plt.ylabel(r"Plasma energy loss rate [W m$^{-3}$]")
plt.legend(loc="upper left")

plt.savefig(path+"/processes_energy.pdf")

plt.xlim([pos[-1]-zoomlength,pos[-1]]) # Zoom to last 1m

plt.savefig(path+"/processes_energy_zoom.pdf")

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

try:
    axarr[1].plot(pos, data["Srec"]*Snorm, label="Recombination (Srec)")
    axarr[1].plot(pos, data["Siz"]*Snorm, label="Ionisation (Siz)")
except:
    pass
axarr[1].set_title(r"Plasma particle loss rate [m$^{-3}$s$^{-1}$]")
axarr[1].legend(loc="upper left")

try:
    axarr[2].plot(pos, data["Frec"]*Fnorm, label="Recombination (Frec)")
    axarr[2].plot(pos, data["Fiz"]*Fnorm, label="Ionisation (Fiz)")
    axarr[2].plot(pos, data["Fcx"]*Fnorm, label="Charge exchange (Fcx)")
except:
    pass
if data["Fel"] is not None:
    axarr[2].plot(pos, data["Fel"]*Enorm, label="Elastic scattering (Fel)")

axarr[2].set_title(r"Plasma momentum transfer to neutrals [m$^{-2}$s$^{-1}$]")
axarr[2].legend(loc="upper left")

try:
    axarr[3].plot(pos, (data["Rrec"]+data["Erec"])*Enorm, label="Recombination (Rrec+Erec)")
    axarr[3].plot(pos, (data["Riz"]+data["Eiz"])*Enorm, label="Ionisation (Riz+Eiz)")
    axarr[3].plot(pos, data["Ecx"]*Enorm, label="Charge exchange (Ecx)")
    axarr[3].plot(pos, data["Rzrad"]*Enorm, label="Impurity radiation (Rzrad)")
except:
    pass
    
if data["Rex"] is not None:
    axarr[3].plot(pos, data["Rex"]*Enorm, label="Hydrogen excitation (Rex)")
if data["Eel"] is not None:
    axarr[3].plot(pos, data["Eel"]*Enorm, label="Elastic scattering (Eel)")
    
axarr[3].set_title(r"Plasma energy loss rate [W m$^{-3}$]")
axarr[3].set_xlabel("Parallel location [m]")
axarr[3].legend(loc="upper left")


plt.savefig(path+"/processes_all.pdf")

axarr[3].set_xlim([pos[-1]-zoomlength,pos[-1]])

plt.savefig(path+"/processes_all_zoom.pdf")

plt.show()
