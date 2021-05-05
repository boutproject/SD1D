#!/usr/bin/env python
#
# Plot profiles of density, temperature and Mach number

from boutdata import collect
from numpy import zeros, sum, sqrt

import sys
import matplotlib
matplotlib.rcParams.update({'font.size': 20})
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    # Print usage information
    print("Usage: {0} path\n e.g. {0} data [zoom]".format(sys.argv[0]))
    sys.exit(1)


path = sys.argv[1]

if len(sys.argv) > 2:
    # Optional second argument is length from target
    zoomlength = float(sys.argv[2])
else:
    zoomlength = 1.0 # 1m default

tind = -1

# Evolving variables, remove extra guard cells so just one each side
P = collect("P", path=path, tind=tind, yguards=True)[-1,0,1:-1,0]
Ne = collect("Ne", path=path, tind=tind, yguards=True)[-1,0,1:-1,0]
NVi = collect("NVi", path=path, tind=tind, yguards=True)[-1,0,1:-1,0]

try:
    Nn = collect("Nn", path=path, tind=tind, yguards=True)[-1,0,1:-1,0]
except IndexError:
    Nn = zeros(Ne.shape)

# Normalisations
nnorm = collect("Nnorm", path=path, tind=tind)
tnorm = collect("Tnorm", path=path, tind=tind)
pnorm = nnorm*tnorm*1.602e-19 # Converts p to Pascals
cs0 = collect("Cs0", path=path)

# electron temperature
Te = (0.5*P/Ne) * tnorm

NVi *= nnorm * cs0
Ne *= nnorm
Nn *= nnorm
P *= pnorm

dy = collect("dy", path=path, yguards=True)[0,1:-1]
n = len(dy)
pos = zeros(n)

# position at the centre of the grid cell
pos[0] = -0.5*dy[1]
pos[1] = 0.5*dy[1]
for i in range(2,n):
    pos[i] = pos[i-1] + 0.5*dy[i-1] + 0.5*dy[i]

def replace_guards(var):
    """
    This in-place replaces the points in the guard cells with the points on the boundary
    
    """
    var[0] = 0.5*(var[0] + var[1])
    var[-1] = 0.5*(var[-1] + var[-2])

replace_guards(pos)
replace_guards(Ne)
replace_guards(Nn)

replace_guards(Te)
Te[-1] = Te[-2] # Zero-gradient Te
#replace_guards(NVi) # Note: NVi in guard cells is Nout*Vout

# ion parallel velocity
Vi = NVi/Ne

# Sound speed
cs = cs0 * sqrt(2.*Te/tnorm)

# Mach number
M = Vi / cs

########################################
# Density, temperature

fig, ax1 = plt.subplots()
ax1.set_xlabel("Position [m]")
ax1.set_ylabel("Temperature [eV]", color="r")
ax1.tick_params('y', colors='r')

ax2 = ax1.twinx()
ax2.set_ylabel(r"Density $m^{-3}$", color='g')
ax2.tick_params('y', colors='g')

ax1.plot(pos, Te, '-r')
ax1.plot(pos[-1], Te[-1], 'or')
ax2.plot(pos, Ne, '-g')
ax2.plot(pos[-1], Ne[-1], 'og')

fig.savefig(path+"/profiles_ne_te.pdf")
fig.savefig(path+"/profiles_ne_te.png")

ax1.set_xlim([pos[-1] - zoomlength, pos[-1]])
ax2.set_xlim([pos[-1] - zoomlength, pos[-1]])

fig.savefig(path+"/profiles_ne_te_zoom.pdf")
fig.savefig(path+"/profiles_ne_te_zoom.png")

print("Saved figures in path '"+path+"'")

plt.show()

########################################
# Density, Mach number

fig, ax1 = plt.subplots()
ax1.set_xlabel("Position [m]")
ax1.set_ylabel("Mach number", color="k")
ax1.tick_params('y', colors='k')

ax2 = ax1.twinx()
ax2.set_ylabel(r"Density $m^{-3}$", color='g')
ax2.tick_params('y', colors='g')

ax1.plot(pos, M, '-k')
ax1.plot(pos[-1], M[-1], 'ok')
ax2.plot(pos, Ne, '-g')
ax2.plot(pos[-1], Ne[-1], 'og')

fig.savefig(path+"/profiles_ne_m.pdf")
fig.savefig(path+"/profiles_ne_m.png")

ax1.set_xlim([pos[-1] - zoomlength, pos[-1]])
ax2.set_xlim([pos[-1] - zoomlength, pos[-1]])

fig.savefig(path+"/profiles_ne_m_zoom.pdf")
fig.savefig(path+"/profiles_ne_m_zoom.png")

########################################
# Plasma, neutral density

fig, ax1 = plt.subplots()

ax1.set_ylabel(r"Plasma density $m^{-3}$", color='g')
ax1.tick_params('y', colors='g')
ax1.set_xlabel("Position [m]")

ax2 = ax1.twinx()
ax2.set_ylabel("Neutral density $m^{-3}$", color="b")
ax2.tick_params('y', colors='b')

ax1.plot(pos, Ne, '-g')
ax1.plot(pos[-1], Ne[-1], 'og')
ax2.plot(pos, Nn, '-b')
ax2.plot(pos[-1], Nn[-1], 'ob')

fig.savefig(path+"/profiles_ne_nn.pdf")
fig.savefig(path+"/profiles_ne_nn.png")

ax1.set_xlim([pos[-1] - zoomlength, pos[-1]])
ax2.set_xlim([pos[-1] - zoomlength, pos[-1]])

fig.savefig(path+"/profiles_ne_nn_zoom.pdf")
fig.savefig(path+"/profiles_ne_nn_zoom.png")

print("Saved figures in path '"+path+"'")

plt.show()



