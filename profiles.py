#!/usr/bin/env python
#
# Plot profiles of density and temperature

from boutdata import collect
from numpy import zeros, sum

import sys
import matplotlib
matplotlib.rcParams.update({'font.size': 20})
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    # Print usage information
    print("Usage: {0} path\n e.g. {0} data".format(sys.argv[0]))
    sys.exit(1)

path = sys.argv[1]

tind = -1

# Evolving variables, remove extra guard cells so just one each side
P = collect("P", path=path, tind=tind, yguards=True)[-1,0,1:-1,0]
Ne = collect("Ne", path=path, tind=tind, yguards=True)[-1,0,1:-1,0]
NVi = collect("NVi", path=path, tind=tind, yguards=True)[-1,0,1:-1,0]

# Normalisations
nnorm = collect("Nnorm", path=path, tind=tind)
tnorm = collect("Tnorm", path=path, tind=tind)
pnorm = nnorm*tnorm*1.602e-19 # Converts p to Pascals
cs0 = collect("Cs0", path=path)

# electron temperature
Te = (0.5*P/Ne) * tnorm

# ion parallel velocity
Vi = (NVi/Ne) * cs0

NVi *= nnorm * cs0
Ne *= nnorm
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
replace_guards(Te)
replace_guards(Vi)

fig, ax1 = plt.subplots()
ax1.set_xlabel("Position [m]")
ax1.set_ylabel("Temperature [eV]", color="r")
ax1.tick_params('y', colors='r')

ax2 = ax1.twinx()
ax2.set_ylabel(r"Density $m^{-3}$", color='b')
ax2.tick_params('y', colors='b')

ax1.plot(pos, Te, '-r')
ax2.plot(pos, Ne, '-b')

fig.savefig("profiles.pdf")
fig.savefig("profiles.png")

ax1.set_xlim([pos[-1] - 1.0, pos[-1]])
ax2.set_xlim([pos[-1] - 1.0, pos[-1]])

fig.savefig("profiles_zoom.pdf")
fig.savefig("profiles_zoom.png")

plt.show()
