# This code was was developed by Leo Werneck

# Load necessary Python modules
from bisect import bisect_left
from numpy import array, empty
from h5py import File as h5open
#from scipy.interpolate import interp1d, interp2d
from scipy.interpolate import interp1d, RectBivariateSpline

class EOS_beq:
    def __init__(self,
                 rho,
                 temperature,
                 rho_of_P_interpolator,
                 P_of_rho_interpolator,
                 eps_of_rho_interpolator,
                 Ye_of_rho_interpolator,
                 units_of_density,
                 units_of_pressure,
                 units_of_spec_int_energy
                 ):
        self.rho                      = rho
        self.temperature              = temperature
        self.rho_of_P_interpolator    = rho_of_P_interpolator
        self.P_of_rho_interpolator    = P_of_rho_interpolator
        self.eps_of_rho_interpolator  = eps_of_rho_interpolator
        self.Ye_of_rho_interpolator   = Ye_of_rho_interpolator
        self.units_of_density         = units_of_density
        self.units_of_pressure        = units_of_pressure
        self.units_of_spec_int_energy = units_of_spec_int_energy

def read_and_slice_EOS_table_at_given_temperature(eos_file_path, T_in_MeV):
    eos_file = h5open(eos_file_path, "r")

    # Read in the ye, P, eps, T, rho, and the neutrino chemical potential
    tab_Ye   =     array(eos_file['ye'       ])[:]
    tab_P    = 10**array(eos_file['logpress' ])[:]
    tab_eps  = 10**array(eos_file['logenergy'])[:] - eos_file['energy_shift'][0]
    tab_temp = 10**array(eos_file['logtemp'  ])[:]
    tab_rho  = 10**array(eos_file['logrho'   ])[:]
    munu     =     array(eos_file['munu'     ])[:]

    # Find temperature index
    T_idx = bisect_left(tab_temp, T_in_MeV)
    print("Slicing table at index %d. T[%d] = %.3lf MeV | Input temperature: %.3lf MeV"%(T_idx, T_idx, tab_temp[T_idx], T_in_MeV))

    # Set 2d interpolators for P and eps
    # as functions of rho and Ye
    P_interp   = RectBivariateSpline(tab_rho,tab_Ye,tab_P[  :,T_idx,:].T, kx = 1, ky = 1)
    eps_interp = RectBivariateSpline(tab_rho,tab_Ye,tab_eps[:,T_idx,:].T, kx = 1, ky = 1)

    # Now, loop over the density and set the array
    Ye_be  = empty(len(tab_rho))
    P_be   = empty(len(tab_rho))
    eps_be = empty(len(tab_rho))

    ye_min = tab_Ye[0]

    for rho_idx in range(len(tab_rho)):
        # Set the interpolator
        Ye_ito_munu_of_T_and_rho = interp1d(munu[:,T_idx,rho_idx],tab_Ye, kind="cubic")
        # Then find the value of Ye which sets munu to zero
        try:
            Ye_be[rho_idx]  = Ye_ito_munu_of_T_and_rho(0.0)
        except:
						# Some high density points do not contain a beta equilibrium solution. Assume that the table minimum is a suitable fix.
            Ye_be[rho_idx] = ye_min
        # Finally, determine P and eps in beta equilibrium
        P_be[  rho_idx] = P_interp(  tab_rho[rho_idx],Ye_be[rho_idx])
        eps_be[rho_idx] = eps_interp(tab_rho[rho_idx],Ye_be[rho_idx])

    # Convert table to code units
    # Set physical constants in cgs units
    G    = 6.6738480e-8   # cm^3/g/s^2
    c    = 2.99792458e10  # cm/s
    Msun = 1.9884e33      # g

    # Set the length factor:
    # .----------------------.
    # | [Msun * G / c^2] = L |
    # .----------------------.
    units_of_length = Msun * G / c**2

    # Set the time factor:
    # .----------------------.
    # | [Msun * G / c^3] = T |
    # .----------------------.
    units_of_time = units_of_length/c

    # Pressure has units of M * L^(-1) * T^(-2), so
    units_of_pressure = Msun * units_of_length**(-1) * units_of_time**(-2)

    # Density has units of M * L^(-3), so
    units_of_density = Msun * units_of_length**(-3)

    # Specific internal energy has units of L^(2) * T^(-2),so
    units_of_spec_int_energy = c**2

    # Convert the table to dimensionless units
    P_be    /= units_of_pressure
    tab_rho /= units_of_density
    eps_be  /= units_of_spec_int_energy

    # Functions to interpolate P and eps
    P_of_rho   = interp1d(P_be   ,tab_rho, kind="cubic")
    rho_of_P   = interp1d(tab_rho,P_be, kind="cubic")
    eps_of_rho = interp1d(tab_rho,eps_be, kind="cubic")
    Ye_of_rho  = interp1d(tab_rho,Ye_be, kind="cubic")

    # Set the EOS tuple
    return EOS_beq(tab_rho, tab_temp[T_idx], P_of_rho, rho_of_P, eps_of_rho, Ye_of_rho,
                   units_of_density, units_of_pressure, units_of_spec_int_energy)

if __name__ == '__main__':
    from sys import argv

    if len(argv) != 3:
        print(f"Correct usage: python3 {argv[0]} <Table Path> <Temperature (in MeV)>")
        exit(1)

    eos = read_and_slice_EOS_table_at_given_temperature(argv[1], float(argv[2]))

    from numpy import savetxt
    savetxt("Ye_beq.txt", list(zip(eos.Ye_of_rho_interpolator(eos.rho))))
    print("Wrote output to file Ye_beq.txt")
