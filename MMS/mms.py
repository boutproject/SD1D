from __future__ import print_function
from __future__ import division
#
# Generate test case using SymPy
# 
#

from boutdata.mms import *

# Define metric tensor
metric = Metric() # Identity

# Define manufactured solution

Ne = 1.0
NVi = 0.0
P = 1.0

Vi = NVi / Ne

# Density
ddt_Ne = -Div_par(Ne*Vi)

# Parallel momentum
ddt_NVi = (
    - Div_par(NVi*Vi)
    - Grad_par(P)
    )
