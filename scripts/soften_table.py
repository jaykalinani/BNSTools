import os
import sys
import bisect
import h5py as h5
import numpy as np
from scipy.ndimage import median_filter

# UNIT CONVERSIONS
## APR stores entropy in units of k_B / baryon, while helmholtz uses erg / K / g. Convert between the two using m_n / k_B.
ENTFAC = 1.2131450484808232e-08
K_TO_MEV = 8.61733326214518e-11 

def stitch_table(extend_tab, soft_tab, irp, irm, nrho_off, chi_rho):
    # In dens merge region
    extend_tab[:, :, irm:irp+1] = chi_rho * extend_tab[:, :, irm:irp+1] + (1 - chi_rho) * soft_tab[:, :, irm-nrho_off:irp-nrho_off+1]

    # In high dens region
    extend_tab[:, :, :irm+1] = soft_tab[:, :, irp-nrho_off+1:]

    return extend_tab

def main(args):

    # Create New Table
    original_tablepath = args.table
    soft_tablepath = args.soft_table
    output_tablepath   = original_tablepath.split(".h5")[0]+"_collapse.h5"

    os.system("cp %s %s"%(original_tablepath,output_tablepath))

    print("Master table file: ",original_tablepath)
    print("Soft   table file: ",soft_tablepath)
    print("Output table file: ",output_tablepath)
    os.remove(output_tablepath)
    eos_file = h5.File(original_tablepath,"r")
    soft_file = h5.File(soft_tablepath,"r")
    out_file = h5.File(output_tablepath, "a")

    # Open Master EOS table, read in values
    tab_Ye   = np.array(eos_file['ye'       ])[:]
    tab_logtemp = np.array(eos_file['logtemp'  ])[:]
    tab_logrho  = np.array(eos_file['logrho'   ])[:]

    dYe = np.abs(tab_Ye[1] - tab_Ye[0])
    dlogT = np.abs(tab_logtemp[1] - tab_logtemp[0])
    dlogrho = np.abs(tab_logrho[1] - tab_logrho[0])
    
    tab_P    = 10**np.array(eos_file['logpress' ])[:]
    tab_shift = eos_file['energy_shift'][0]
    tab_eps  = 10**np.array(eos_file['logenergy'])[:] - tab_shift 
    tab_S = np.array(eos_file['entropy'])[:]
    tab_cs2 = np.array(eos_file['cs2'])[:]
    tab_muhat = np.array(eos_file['muhat'])[:]
    tab_mu_n = np.array(eos_file['mu_n'])[:]
    tab_mu_p = np.array(eos_file['mu_p'])[:]
    tab_mu_e = np.array(eos_file['mu_e'])[:]
    tab_munu = np.array(eos_file['munu'])[:]
    tab_Xn = np.array(eos_file['Xn'])[:]
    tab_Xp = np.array(eos_file['Xp'])[:]
    tab_Xa = np.array(eos_file['Xa'])[:]
    tab_Xh = np.array(eos_file['Xh'])[:]
    tab_Abar = np.array(eos_file['Abar'])[:]
    tab_Zbar = np.array(eos_file['Zbar'])[:]
    tab_gamma = np.array(eos_file['gamma'])[:]
    tab_dedt = np.array(eos_file['dedt'])[:]
    tab_dpdrhoe = np.array(eos_file['dpdrhoe'])[:]
    tab_dpderho = np.array(eos_file['dpderho'])[:]

    # Open Soft EOS table, read in values
    soft_tab_Ye   = np.array(soft_file['ye'       ])[:]
    soft_tab_logtemp = np.array(soft_file['logtemp'  ])[:]
    soft_tab_logrho  = np.array(soft_file['logrho'   ])[:]

    soft_dYe = np.abs(soft_tab_Ye[1] - soft_tab_Ye[0])
    soft_dlogT = np.abs(soft_tab_logtemp[1] - soft_tab_logtemp[0])
    soft_dlogrho = np.abs(soft_tab_logrho[1] - soft_tab_logrho[0])
    
    soft_tab_P    = 10**np.array(soft_file['logpress' ])[:]
    soft_tab_shift = soft_file['energy_shift'][0]
    soft_tab_eps  = 10**np.array(soft_file['logenergy'])[:] - soft_tab_shift 
    soft_tab_S = np.array(soft_file['entropy'])[:]
    soft_tab_cs2 = np.array(soft_file['cs2'])[:]
    soft_tab_muhat = np.array(soft_file['muhat'])[:]
    soft_tab_mu_n = np.array(soft_file['mu_n'])[:]
    soft_tab_mu_p = np.array(soft_file['mu_p'])[:]
    soft_tab_mu_e = np.array(soft_file['mu_e'])[:]
    soft_tab_munu = np.array(soft_file['munu'])[:]
    soft_tab_Xn = np.array(soft_file['Xn'])[:]
    soft_tab_Xp = np.array(soft_file['Xp'])[:]
    soft_tab_Xa = np.array(soft_file['Xa'])[:]
    soft_tab_Xh = np.array(soft_file['Xh'])[:]
    soft_tab_Abar = np.array(soft_file['Abar'])[:]
    soft_tab_Zbar = np.array(soft_file['Zbar'])[:]
    soft_tab_gamma = np.array(soft_file['gamma'])[:]
    soft_tab_dedt = np.array(soft_file['dedt'])[:]
    soft_tab_dpdrhoe = np.array(soft_file['dpdrhoe'])[:]
    soft_tab_dpderho = np.array(soft_file['dpderho'])[:]

    # Check table compatibility
    og_shape = (len(tab_Ye), len(tab_logtemp), len(tab_logrho))
    assert(len(tab_logtemp) == len(soft_tab_logtemp))
    assert(len(tab_Ye) == len(soft_tab_Ye))
    nrho_offset = len(tab_logrho) - len(soft_tab_logrho)

    # STITCHING PARAMETERS
    rho_stitch = args.rho_stitch # log10(rho in g cm^-3)
    width = args.width
    tol = args.tol

    ## Transition region bounds. See MERGE in SROEOS User Guide
    rho_t_plus = rho_stitch - width * np.arctanh(2 * tol - 1)
    rho_t_minus =  rho_stitch + width * np.arctanh(2 * tol - 1)
    irho_plus = bisect.bisect_left(tab_logrho, rho_t_plus)
    irho_minus = bisect.bisect_left(tab_logrho, rho_t_minus)

    ## Stitch param check
    print("rho_t_+ = {:.6f}; rho_t_- = {:.6f}".format(rho_t_plus, rho_t_minus))
    print("irho_+ = {:d}; irho_- = {:d}".format(irho_plus, irho_minus))

    # INITIALIZE NEW TABLES TO MASTER EOS
    P = tab_P
    eps = tab_eps  
    S = tab_S 
    cs2 = tab_cs2
    muhat = tab_muhat
    mu_n = tab_mu_n
    mu_p = tab_mu_p
    mu_e = tab_mu_e
    munu = tab_munu
    Xn = tab_Xn
    Xp = tab_Xp
    Xa = tab_Xa
    Xh = tab_Xh
    Abar = tab_Abar
    Zbar = tab_Zbar
    gamma = tab_gamma
    dedt = tab_dedt
    dpdrhoe = tab_dpdrhoe
    dpderho = tab_dpderho

    # ADJUST HIGH DENS REGION

    ## Set stitching functions
    chi_rho = 0.5 * (1 + np.tanh((tab_logrho[irho_minus:irho_plus+1] - rho_stitch) / width))
    chi_rho = np.tile(chi_rho, (og_shape[0], og_shape[1], 1))

    # Stitch tables
    P = stitch_table(P, soft_tab_P, irho_plus, irho_minus, nrho_offset, chi_rho)
    eps = stitch_table(eps, soft_tab_eps, irho_plus, irho_minus, nrho_offset, chi_rho)
    S = stitch_table(S, soft_tab_S, irho_plus, irho_minus, nrho_offset, chi_rho)
    cs2 = stitch_table(cs2, soft_tab_cs2, irho_plus, irho_minus, nrho_offset, chi_rho)
    muhat = stitch_table(muhat, soft_tab_muhat, irho_plus, irho_minus, nrho_offset, chi_rho)
    mu_n = stitch_table(mu_n, soft_tab_mu_n, irho_plus, irho_minus, nrho_offset, chi_rho)
    mu_p = stitch_table(mu_p, soft_tab_mu_p, irho_plus, irho_minus, nrho_offset, chi_rho)
    mu_e = stitch_table(mu_e, soft_tab_mu_e, irho_plus, irho_minus, nrho_offset, chi_rho)
    munu = stitch_table(munu, soft_tab_munu, irho_plus, irho_minus, nrho_offset, chi_rho)
    Xn = stitch_table(Xn, soft_tab_Xn, irho_plus, irho_minus, nrho_offset, chi_rho)
    Xp = stitch_table(Xp, soft_tab_Xp, irho_plus, irho_minus, nrho_offset, chi_rho)
    Xa = stitch_table(Xa, soft_tab_Xa, irho_plus, irho_minus, nrho_offset, chi_rho)
    Xh = stitch_table(Xh, soft_tab_Xh, irho_plus, irho_minus, nrho_offset, chi_rho)
    Abar = stitch_table(Abar, soft_tab_Abar, irho_plus, irho_minus, nrho_offset, chi_rho)
    Zbar = stitch_table(Zbar, soft_tab_Zbar, irho_plus, irho_minus, nrho_offset, chi_rho)
    gamma =stitch_table(gamma, soft_tab_gamma, irho_plus, irho_minus, nrho_offset, chi_rho)
    dedt = stitch_table(dedt, soft_tab_dedt, irho_plus, irho_minus, nrho_offset, chi_rho)
    dpdrhoe = stitch_table(dpdrhoe, soft_tab_dpdrhoe, irho_plus, irho_minus, nrho_offset, chi_rho)
    dpderho = stitch_table(dpderho, soft_tab_dpderho, irho_plus, irho_minus, nrho_offset, chi_rho)

    # Update new HDF5
    print("Outputting to file")
    out_file.create_dataset("pointsrho", data=np.array([len(tab_logrho)]), shape=(1,))
    out_file.create_dataset("pointstemp", data=np.array([len(tab_logtemp)]), shape=(1,))
    out_file.create_dataset("pointsye", data=np.array([len(tab_Ye)]), shape=(1,))
    out_file.create_dataset("logrho", data=tab_logrho)
    out_file.create_dataset("logtemp", data=tab_logtemp)
    out_file.create_dataset("ye", data=tab_Ye)
    out_file.create_dataset("logpress", data=np.log10(P))
    out_file.create_dataset("entropy", data=S)
    out_file.create_dataset("logenergy", data=np.log10(eps + tab_shift))
    out_file.create_dataset("energy_shift", data=np.array([tab_shift]), shape=(1,))
    out_file.create_dataset("muhat", data=muhat)
    out_file.create_dataset("mu_n", data=mu_n)
    out_file.create_dataset("mu_p", data=mu_p)
    out_file.create_dataset("mu_e", data=mu_e)
    out_file.create_dataset("munu", data=munu)
    out_file.create_dataset("dedt", data=dedT)
    out_file.create_dataset("dpdrhoe", data=dPdrhoe)
    out_file.create_dataset("dpderho", data=dPde)
    out_file.create_dataset("cs2", data=cs2)
    out_file.create_dataset("gamma", data=gamma)
    out_file.create_dataset("Xn", data=Xn)
    out_file.create_dataset("Xp", data=Xp)
    out_file.create_dataset("Xa", data=Xa)
    out_file.create_dataset("Xh", data=Xh)
    out_file.create_dataset("Abar", data=Abar)
    out_file.create_dataset("Zbar", data=Zbar)

    # Close files
    out_file.close()
    eos_file.close()
    soft_file.close()

    return
 
if __name__ == '__main__':
    #from sys import argv
    import argparse
    # -----------------------------------------------------------------------------
    parser = argparse.ArgumentParser()
    # -----------------------------------------------------------------------------
    parser.add_argument("-t", "--table", dest="table", required=True,
        help="EOS table")
    parser.add_argument("-s", "--soft-table", dest="soft_table", required=True,
        help="Soft EOS table")
    parser.add_argument("--stitch-dens", dest="rho_stitch", required=True, type=float,
        help="Transition density")
    parser.add_argument("--width", dest="width", required=True, type=float,
        help="Width of transition region")
    parser.add_argument("--tol", dest="tol", required=True, type=float,
        help="Tolerance of transition region")

    args = parser.parse_args()
    main(args)
    print("Done :)")
