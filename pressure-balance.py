from boutdata import collect
from numpy import zeros, sum, sqrt, sign

import sys
import matplotlib
matplotlib.rcParams.update({'font.size': 18})
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    # Print usage information
    print("Usage: {0} path\n e.g. {0} data".format(sys.argv[0]))
    sys.exit(1)

    
path = sys.argv[1]
tind = -1

# Evolving variables, remove extra guard cells so just one each side
p = collect("P", path=path, tind=tind, yguards=True)[-1,0,1:-1,0]
pn = collect("Pn", path=path, tind=tind, yguards=True)[-1,0,1:-1,0]

nvi = collect("NVi", path=path, tind=tind, yguards=True)[-1,0,1:-1,0]
ne = collect("Ne", path=path, tind=tind, yguards=True)[-1,0,1:-1,0]
nvn = collect("NVn", path=path, tind=tind, yguards=True)[-1,0,1:-1,0]
nn = collect("Nn", path=path, tind=tind, yguards=True)[-1,0,1:-1,0]


# Normalisations
nnorm = collect("Nnorm", path=path, tind=tind)
tnorm = collect("Tnorm", path=path, tind=tind)
pnorm = nnorm*tnorm*1.602e-19 # Converts p to Pascals
cs0 = collect("Cs0", path=path)


p *= pnorm
pn *= pnorm

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

replace_guards(p)
replace_guards(pn)
replace_guards(ne)
replace_guards(nn)

dynamic_p = (nvi**2/ne) * pnorm #* sign(nvi)
dynamic_n = (nvn**2/nn) * pnorm #* sign(nvn)

########################################
# Static pressures

fig, ax1 = plt.subplots()
ax1.set_xlabel("Position [m]")
ax1.set_ylabel("Plasma pressure [Pa]", color="b")
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
ax2.set_ylabel(r"Neutral pressure [Pa]", color='g')
ax2.tick_params('y', colors='g')

ax1.plot(pos, p, '-r')
ax1.plot(pos[0], p[0], 'ob')
#ax1.annotate(r"$p = %.2e$Pa" % (p[0],), xy=(pos[0], p[0]), textcoords='data')

ax2.plot(pos, pn, '-g')
ax2.plot(pos[-1], pn[-1], 'og')

#ax2.annotate(r"$p_n = %.2e$Pa" % (pn[-1],), xy=(pos[-1], pn[-1]), textcoords='data')

fig.savefig(path+"/pressure_balance.pdf")
fig.savefig(path+"/pressure_balance.png")

ax1.set_xlim([pos[-1] - 1.0, pos[-1]])
ax2.set_xlim([pos[-1] - 1.0, pos[-1]])

fig.savefig(path+"/pressure_balance_zoom.pdf")
fig.savefig(path+"/pressure_balance_zoom.png")

plt.show()

########################################
# Static, dynamic


fig, ax1 = plt.subplots()
ax1.set_xlabel("Position [m]")
ax1.set_ylabel("Pressure [Pa]")

ax1.plot(pos, p, '-b', label="Static plasma")
ax1.plot(pos, pn, '--b', label="Static neutral")

ax1.plot(pos, dynamic_p, '-r', label="Dynamic plasma")
ax1.plot(pos, dynamic_n, '--r', label="Dynamic neutral")

ax1.plot(pos, p + pn + dynamic_p + dynamic_n, '-k', label="Total")

l1 = ax1.legend()
l1.set_zorder(0)   # Put legend at the back

print("Upstream pressure: {0}".format(p[0]))
print("""
Target pressures
---
Static plasma: {0}
Dynamic plasma: {1}
Static neutral: {2}
Dynamic neutral: {3}
---
Total plasma: {4} (fmomloss = {5})
Total: {6}
---\n""".format(p[-1],
                dynamic_p[-1],
                pn[-1],
                dynamic_n[-1],
                p[-1] + dynamic_p[-1], 1.0 - (p[-1] + dynamic_p[-1])/p[0],
                p[-1] + dynamic_p[-1] + pn[-1] + dynamic_n[-1]))

fig.savefig(path+"/pressure_plasma.pdf")
fig.savefig(path+"/pressure_plasma.png")

ax1.set_xlim([pos[-1] - 2.0, pos[-1]])

fig.savefig(path+"/pressure_plasma_zoom.pdf")
fig.savefig(path+"/pressure_plasma_zoom.png")

plt.show()
