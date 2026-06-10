# This Python module is inteded to be used on EOS tables downloaded from
# stellarcollapse.org. It adjusts the baryon density rho so that the
# specific internal energy is always positive. It also cleans up the
# spurious sound speed table entries.
#
# Author: Leonardo Werneck
#
# This module is based on the stellarcollapse.py module from David Radice.
# Source: https://bitbucket.org/FreeTHC/thcextra/src/master/tools/eos_tables/import/

# Step 0a: Initialize core Python modules
import os
import sys
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import Image

# Step 0b: Get user input
# Step 0b.i: Check correct usage
if len(sys.argv) != 2:
    print("(ERROR) Correct usage is: python [this script] [input eos table path]")
    sys.exit(1)

# Step 0b.ii: Input table path
original_tablepath = sys.argv[1]
compute_relativistic_cs2 = False
if len(sys.argv) > 2:
    if sys.argv[2].lower() == "true" or sys.argv[2].lower() == "yes":
        compute_relativistic_cs2 = True

# Step 0b.iii: Output table path
output_tablepath   = original_tablepath.split(".h5")[0]+"_adj_mineps0.h5"

# Step 0b.iv: Make a copy of the original table. This will be used as the output file.
os.system("cp %s %s"%(original_tablepath,output_tablepath))

# Step 0b.v: Print basic info to the user
print("Input  table file: ",original_tablepath)
print("Output table file: ",output_tablepath)

# Step 0c: Define physical constants
# Change mass definition for APR4
# Step 0c.i: Free Neutron Mass (MeV)
MN_MEV      = 939.5654133
# Step 0c.ii: Vacuum speed of light
C_CGS         = 29979245800.0
# Step 0c.iii: Conversion factor from eV to cgs
EV_CGS        = 1.6021765649999999e-12
# Step 0c.iv: Conversion factor from MeV to cgs
MEV_CGS       = 1e6 * EV_CGS
# Step 0c.v: Unified atomic mass unit (cgs)
MN_CGS      = MN_MEV * MEV_CGS / (C_CGS**2)

# Step 1: Read in the output table
table = h5.File(output_tablepath,"r+")

# Step 2: Set hydro quantities
# Step 2a: Set rho: table entry is log10(rho)
rho = 10**np.array(table["logrho"])[:]

# Step 2b: Set the pressure: table entry is log10(press)
press = 10**np.array(table["logpress"])[:]

# Step 2c: Table entry is log10(eps + eps_{0}), set eps s.t. min(eps + eps_{0}) = 0
eps0 = table["energy_shift"][0]
eps_tab  = 10**np.array(table["logenergy"])[:]
eps  = eps_tab - eps0
eps_min = np.min(eps)

# Step 2d: Set cs2
cs2 = np.array(table["cs2"])[:]

# Step 3: Set the new mass factor
# Step 3a: Set the atomic mass unit (in MeV)
mu = MN_MEV

# Step 3b: Set the atomic mass unit (in cgs)
mu_cgs = MN_CGS

# Step 3c: Compute m_f = (1 - (eps_{0} - eps_min))*m_u
mf = MN_MEV*(1.0 + eps_min/C_CGS**2)

# Step 3d: Convert mf from MeV to cgs units
mf_cgs = (mf*MEV_CGS)/(C_CGS**2)
print("THE ADJUSTED BARYON MASS IS {:e} GRAMS.".format(mf_cgs))

# Step 4: Compute the total energy density: e = rho(1+eps)
energy = rho*(1+eps/(C_CGS**2))

# Step 5: Adjust rho
#rho *= (mf_cgs/mu_cgs)
rho *= (1.0 + eps_min/C_CGS**2)

# Step 6: Adjust eps
eps = (energy/rho - 1)*C_CGS**2

# Step 7: Adjust cs2
# Step 7a: Compute h = c^{2} + eps + P/rho
h    = C_CGS**2 + eps + press/rho

# Step 7b: Compute the relativistic sound speed
cs2 /= h

# Step 7c: Fix NaN entries
cs2[np.isnan(cs2)] = 1.0-1e-15

# Step 7d: Fix infinite entries
cs2[np.isinf(cs2)] = 1.0-1e-15

# Step 7e: Fix superluminal entries
cs2[cs2 > 1]   = 1.0

# Step 7f: Fix negative c_{s}^{2} entries
cs2[cs2 < 0] = 0.0

# Step 8: Output the adjusted table
# Step 8a: Density
table["logrho"      ][...] = np.log10(rho)
# Step 8b: Specific internal energy
table["logenergy"   ][...] = np.log10(eps + 100) # To avoid log(0)
# Step 8c: Sound speed
table["cs2"         ][...] = cs2*h
# Step 8d: Energy shift
table["energy_shift"][...] = 100
# Step 8e: Close the table file, which should now have been adjusted.
table.close()

# Step 9: Validation
# Step 9a: Read in table quantities
with h5.File(output_tablepath,"r") as table:
    nrho  = table["pointsrho"][0]
    nye  = table["pointsye"][0]
    ntemp  = table["pointstemp"][0]
    rho   = 10**np.array(table["logrho"      ])[:]
    press = 10**np.array(table["logpress"    ])[:]
    eps0  =              table["energy_shift"][0]
    eps   = 10**np.array(table["logenergy"   ])[:] - eps0
    cs2   =     np.array(table["cs2"         ])[:]

# Step 9b: Compute h
h    = C_CGS**2 + eps + press/rho

# Step 9c: Compute the relativistic sound speed
cs2 /= h

# Step 9d: Convert h to dimensionless units
h   /= C_CGS**2

# Step 9e: Print validation notes
print("\nValidating adjusted table:")
print(" h   > 1:",np.all(h>1))
print("eps >= 0:",np.all(eps>=0))
print("eps0 = 0:",eps0==0)
print("cs2 >= 0:",np.all(cs2>=0))
print("cs2 <= 1:",np.all(cs2/h<=1))

# All done!
print("\nAll done!")
