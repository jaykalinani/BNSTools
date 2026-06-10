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

def stitch(chi, tab, helm):
    return chi * tab + (1 - chi) * helm

def extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, orig_tab):
    # Assume composition does not change in low density, temperature region (see Hayashi et al 2023)
    new_tab = np.full(new_shape, np.nan)
    new_tab[:, num_new_pts_T:, num_new_pts_rho:] = orig_tab
    for iYe in range(new_shape[0]):
        for iT in range(new_shape[1]):
            for irho in range(new_shape[2]):
                if np.isnan(new_tab[iYe, iT, irho]):
                    irho_star = max(0, irho - num_new_pts_rho)
                    iT_star = max(0, iT - num_new_pts_T)
                    new_tab[iYe, iT, irho] = orig_tab[iYe, iT_star, irho_star]
    return new_tab

def extend_std_table_at_Ye(iYe, nT_new, nrho_new, irp, irm, iTp, iTm, chi_rho, chi_temp, extend_tab, orig_tab, h_tab, do_temp_ext):
    # In dens merge region
    extend_tab[iYe, nT_new:, irm:irp+1] = chi_rho * orig_tab[iYe, :, irm-nrho_new:irp-nrho_new+1] + (1 - chi_rho) * h_tab[nT_new:, irm:irp+1]

    # In low dens region
    extend_tab[iYe,nT_new:, :irm+1] = h_tab[nT_new:, :irm+1]

    if do_temp_ext:
        # In temp merge region
        extend_tab[iYe, iTm:iTp+1, :] = chi_temp * extend_tab[iYe, iTm:iTp+1, :] + (1 - chi_temp) * h_tab[iTm:iTp+1, :]

        # low temp region
        extend_tab[iYe,:iTm+1,:] = h_tab[:iTm+1, :]

    return extend_tab

def deriv_median_filter(tab, width, thresh):
    tab_cpy = tab
    # Get median value of neighbors at each point
    filt_tab = median_filter(tab, size=width, mode="nearest")

    # If value differs from neighbors too much, set value to median
    bad_pts = (np.abs(tab - filt_tab) / np.abs(filt_tab)) > thresh

    arr_iter = np.nditer(tab, flags=['multi_index'])
    for val in arr_iter:
        if bad_pts[arr_iter.multi_index]:
            tab_cpy[arr_iter.multi_index] = filt_tab[arr_iter.multi_index]

    return tab_cpy

def main(args):
    try:
        import helmholtz
    except:
        print("Did you install python-helmholtz?")
        return
 
    # Set New Table Parameters

    do_temp_extend = args.do_temp_extend
    print("TEMP EXTEND???")
    print(do_temp_extend)

    # Final Params
    logrho_min = args.rhomin
    logT_min = args.Tmin

    # Create New Table
    original_tablepath = args.table
    output_tablepath   = original_tablepath.split(".h5")[0]+"_ext.h5"

    os.system("cp %s %s"%(original_tablepath,output_tablepath))

    print("Input  table file: ",original_tablepath)
    print("Output table file: ",output_tablepath)
    os.remove(output_tablepath)
    eos_file = h5.File(original_tablepath,"r")
    out_file = h5.File(output_tablepath, "a")

    # Open EOS file, read in values
    tab_Ye   = np.array(eos_file['ye'       ])[:]
    tab_logtemp = np.array(eos_file['logtemp'  ])[:]
    tab_logrho  = np.array(eos_file['logrho'   ])[:]

    dYe = np.abs(tab_Ye[1] - tab_Ye[0])
    dlogT = np.abs(tab_logtemp[1] - tab_logtemp[0])
    dlogrho = np.abs(tab_logrho[1] - tab_logrho[0])
    
    if args.og_rhomin != 999:
        irho_ogmin = bisect.bisect_left(tab_logrho, args.og_rhomin)
        tab_logrho = tab_logrho[irho_ogmin:]
    else:
        irho_ogmin = 0

    if args.og_Tmin != 999:
        iT_ogmin = bisect.bisect_left(tab_logtemp, args.og_Tmin)
        tab_logtemp = tab_logtemp[iT_ogmin:]
    else:
        iT_ogmin = 0

    tab_P    = 10**np.array(eos_file['logpress' ])[:, iT_ogmin:, irho_ogmin:]
    tab_shift = eos_file['energy_shift'][0]
    tab_eps  = 10**np.array(eos_file['logenergy'])[:, iT_ogmin:, irho_ogmin:] - tab_shift 
    tab_S = np.array(eos_file['entropy'])[:, iT_ogmin:, irho_ogmin:]
    tab_cs2 = np.array(eos_file['cs2'])[:, iT_ogmin:, irho_ogmin:]
    tab_muhat = np.array(eos_file['muhat'])[:, iT_ogmin:, irho_ogmin:]
    tab_mu_n = np.array(eos_file['mu_n'])[:, iT_ogmin:, irho_ogmin:]
    tab_mu_p = np.array(eos_file['mu_p'])[:, iT_ogmin:, irho_ogmin:]
    tab_mu_e = np.array(eos_file['mu_e'])[:, iT_ogmin:, irho_ogmin:]
    tab_munu = np.array(eos_file['munu'])[:, iT_ogmin:, irho_ogmin:]
    tab_Xn = np.array(eos_file['Xn'])[:, iT_ogmin:, irho_ogmin:]
    tab_Xp = np.array(eos_file['Xp'])[:, iT_ogmin:, irho_ogmin:]
    tab_Xa = np.array(eos_file['Xa'])[:, iT_ogmin:, irho_ogmin:]
    tab_Xh = np.array(eos_file['Xh'])[:, iT_ogmin:, irho_ogmin:]
    tab_Abar = np.array(eos_file['Abar'])[:, iT_ogmin:, irho_ogmin:]
    tab_Zbar = np.array(eos_file['Zbar'])[:, iT_ogmin:, irho_ogmin:]
    tab_gamma = np.array(eos_file['gamma'])[:, iT_ogmin:, irho_ogmin:]

    # We want the new table to line up with the old one where they overlap
    num_new_pts_rho = int((tab_logrho[0] - logrho_min) / dlogrho)
    rho_min_adj = tab_logrho[0] - num_new_pts_rho * dlogrho
    total_pts_rho = num_new_pts_rho + len(tab_logrho)
    rho_space_final = np.linspace(rho_min_adj, tab_logrho[-1], total_pts_rho)

    if do_temp_extend:
        num_new_pts_T = int((tab_logtemp[0] - logT_min) / dlogT)
        T_min_adj = tab_logtemp[0] - num_new_pts_T * dlogT
        total_pts_T = num_new_pts_T + len(tab_logtemp)
        T_space_final = np.linspace(T_min_adj, tab_logtemp[-1], total_pts_T)

    else:
        num_new_pts_T = 0
        T_min_adj = tab_logtemp[0]
        total_pts_T = len(tab_logtemp)
        T_space_final = np.linspace(T_min_adj, tab_logtemp[-1], total_pts_T)

    new_shape = (len(tab_Ye), len(T_space_final), len(rho_space_final))

    # STITCHING PARAMETERS
    rho_stitch = args.rho_stitch # log10(rho in g cm^-3)
    T_stitch = args.T_stitch  # log10(T in MeV)
    width = args.width
    tol = args.tol

    ## Transition region bounds. See MERGE in SROEOS User Guide
    rho_t_plus = rho_stitch - width * np.arctanh(2 * tol - 1)
    rho_t_minus =  rho_stitch + width * np.arctanh(2 * tol - 1)
    irho_plus = bisect.bisect_left(rho_space_final, rho_t_plus)
    irho_minus = bisect.bisect_left(rho_space_final, rho_t_minus)

    if do_temp_extend:
        T_t_plus = T_stitch - width * np.arctanh(2 * tol - 1)
        T_t_minus =  T_stitch + width * np.arctanh(2 * tol - 1)
        iT_plus = bisect.bisect_left(T_space_final, T_t_plus)
        iT_minus = bisect.bisect_left(T_space_final, T_t_minus)

    else:
        # Dummy values
        T_t_plus = T_stitch
        T_t_minus =  T_stitch
        iT_plus = 0
        iT_minus = 0

    rho_space_orig = np.arange(tab_logrho[0], tab_logrho[-1] + dlogrho, dlogrho)
    T_space_orig = np.arange(tab_logtemp[0], tab_logtemp[-1] + dlogT, dlogT)

    ## Stitch param check
    print("new logrho_min = {:.6f}; new logT_min = {:.6f}".format(rho_min_adj, T_min_adj))
    print("total pts rho = {:d}; total pts T = {:d}".format(len(rho_space_final), len(T_space_final)))
    print("rho_t_+ = {:.6f}; rho_t_- = {:.6f}; T_t_+ = {:.6f}; T_t_- = {:.6f}".format(rho_t_plus, rho_t_minus, T_t_plus, T_t_minus))
    print("irho_+ = {:d}; irho_- = {:d}; iT_+ = {:d}; iT_- = {:d}; ".format(irho_plus, irho_minus, iT_plus, iT_minus))
    print("new_pts_rho = {:d}; new_pts_T = {:d}".format(num_new_pts_rho, num_new_pts_T))
    print("Consistency check: rho[irho_+] = {:.4e}; rho[irho_-] = {:.4e}; T[iT_+] = {:.4e}; T[iT_-] = {:.4e}; ".format(rho_space_final[irho_plus], rho_space_final[irho_minus], T_space_final[iT_plus], T_space_final[iT_minus]))
    if (iT_minus < num_new_pts_T) or (irho_minus < num_new_pts_rho) or (iT_plus > len(T_space_final)) or (irho_plus > len(rho_space_final)):
        print("Error: transition region out of bounds.")
        return 
    

    # CREATE NEW COMPOSITION TABLES
    muhat = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_muhat)
    mu_n = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_mu_n)
    mu_p = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_mu_p)
    mu_e = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_mu_e)
    munu = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_munu)
    Xp = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_Xp)
    Xn = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_Xn)
    Xa = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_Xa)
    Xh = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_Xh)
    Abar = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_Abar)
    Zbar = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_Zbar)

    # INITIALIZE EXTENDED TABLES
    P = np.empty(new_shape)
    S = np.empty(new_shape)
    cs2 = np.empty(new_shape)
    gamma = np.empty(new_shape)
    # eps needs special treatment due to eps_nuc (see Hayashi et al 2023)
    eps = np.full(new_shape, np.nan)

    # FILL ORIGINAL REGION
    P[:, num_new_pts_T:, num_new_pts_rho:] = tab_P
    eps[:, num_new_pts_T:, num_new_pts_rho:] = tab_eps
    S[:, num_new_pts_T:, num_new_pts_rho:] = tab_S
    cs2[:, num_new_pts_T:, num_new_pts_rho:] = tab_cs2
    gamma[:, num_new_pts_T:, num_new_pts_rho:] = tab_gamma

    # EXTEND REMAINING TABLES

    ## Set stitching functions
    chi_rho = 0.5 * (1 + np.tanh((rho_space_final[irho_minus:irho_plus+1] - rho_stitch) / width))
    chi_rho = np.tile(chi_rho, (new_shape[1] - num_new_pts_T, 1))

    chi_temp = 0.5 * (1 + np.tanh((T_space_final[iT_minus:iT_plus+1] - T_stitch) / width))
    chi_temp = np.tile(chi_temp, (new_shape[2], 1)).T

    ## Loop over Y_e
    for iYe, valYe in enumerate(tab_Ye):
        print("Calculating Ye = {:.4f}".format(valYe))

        ## Prepare input arrays for Helmholtz EOS
        helm_rhospace = 10**rho_space_final
        ## Helmholtz EOS has maximum rho * Y_e = 1e15. This will only be violated in the APREOS region or the forbidden high density, low temp region
        for irho in range(new_shape[2]):
            if helm_rhospace[irho] * valYe > 1e15:
                helm_rhospace[irho] = (1 - 1e-3) * 1e15 / valYe
        helm_rhospace = np.tile(helm_rhospace, (new_shape[1], 1))

        helm_Tspace = np.tile(10**T_space_final / K_TO_MEV, (new_shape[2], 1)).T

        Abar_tot = Xp[iYe,:,:] + Xn[iYe,:,:] + Xa[iYe,:,:] * 4 + Xh[iYe,:,:] * Abar[iYe,:,:]

        # Create Helmholtz 2D table
        h = helmholtz.helmeos(dens=helm_rhospace, temp=helm_Tspace, abar=Abar_tot, zbar=Abar_tot*valYe) # Y_e = zbar / abar 

        # Extend tables
        P = extend_std_table_at_Ye(iYe, num_new_pts_T, num_new_pts_rho, irho_plus, irho_minus, iT_plus, iT_minus, chi_rho, chi_temp, P, tab_P, h.ptot, do_temp_extend)
        S = extend_std_table_at_Ye(iYe, num_new_pts_T, num_new_pts_rho, irho_plus, irho_minus, iT_plus, iT_minus, chi_rho, chi_temp, S, tab_S, ENTFAC * h.stot, do_temp_extend)
        cs2 = extend_std_table_at_Ye(iYe, num_new_pts_T, num_new_pts_rho, irho_plus, irho_minus, iT_plus, iT_minus, chi_rho, chi_temp, cs2, tab_cs2, h.cs**2, do_temp_extend)
        gamma = extend_std_table_at_Ye(iYe, num_new_pts_T, num_new_pts_rho, irho_plus, irho_minus, iT_plus, iT_minus, chi_rho, chi_temp, gamma, tab_gamma, h.gam1, do_temp_extend)

        # eps needs special treatment due to eps_nuc (see Hayashi et al 2023)
        for iT in range(new_shape[1]):
            iT_star = max(0, iT - num_new_pts_T)
            for irho in range(new_shape[2]):
                irho_star = max(0, irho - num_new_pts_rho)
                if np.isnan(eps[iYe, iT, irho]):
                    eps_nuc = tab_eps[iYe, iT_star, irho_star] - h.etot[max(iT, num_new_pts_T), max(irho, num_new_pts_rho)]
                    eps[iYe, iT, irho] = h.etot[iT, irho] + eps_nuc

    ## End of Ye Loop

    # CALCULATE DERIVATIVES
    print("Calculating Derivatives")
    dT = 10**T_space_final
    drho = 10**rho_space_final

    dedT = np.gradient(eps, dT, axis=1)
    dPdT = np.gradient(P, dT, axis=1)
    dPde = dPdT / dedT
    dPdrho = np.gradient(P, drho, axis=2)
    dedrho = np.gradient(eps, drho, axis=2)
    dPdrhoe = dPdrho - dPdT * dedrho / dedT

    # Default to original values in unmodified region

    ## Load in original data
    tab_dedt = np.array(eos_file["dedt"])[:, iT_ogmin:, irho_ogmin:]
    tab_dPde = np.array(eos_file["dpderho"])[:, iT_ogmin:, irho_ogmin:]
    tab_dPdrhoe = np.array(eos_file["dpdrhoe"])[:, iT_ogmin:, irho_ogmin:]

    dedT[:,iT_plus:,irho_plus:] = tab_dedt[:,iT_plus-num_new_pts_T:,irho_plus-num_new_pts_rho:] 
    dPde[:,iT_plus:,irho_plus:] = tab_dPde[:,iT_plus-num_new_pts_T:,irho_plus-num_new_pts_rho:] 
    dPdrhoe[:,iT_plus:,irho_plus:] = tab_dPdrhoe[:,iT_plus-num_new_pts_T:,irho_plus-num_new_pts_rho:] 

    # Perform median filter as is done in Jonah Miller's Singularity EOS code
    MF_W = 3 ## vals from Singularity
    MF_THRESH = 10
    dedT = deriv_median_filter(dedT, MF_W, MF_THRESH)
    dPde = deriv_median_filter(dPde, MF_W, MF_THRESH)
    dPdrhoe = deriv_median_filter(dPdrhoe, MF_W, MF_THRESH)

    # Update new HDF5
    print("Outputting to file")
    out_file.create_dataset("pointsrho", data=np.array([len(rho_space_final)]), shape=(1,))
    out_file.create_dataset("pointstemp", data=np.array([len(T_space_final)]), shape=(1,))
    out_file.create_dataset("pointsye", data=np.array([len(tab_Ye)]), shape=(1,))
    out_file.create_dataset("logrho", data=rho_space_final)
    out_file.create_dataset("logtemp", data=T_space_final)
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

    return
 
if __name__ == '__main__':
    #from sys import argv
    import argparse
    # -----------------------------------------------------------------------------
    parser = argparse.ArgumentParser()
    # -----------------------------------------------------------------------------
    parser.add_argument("-t", "--table", dest="table", required=True,
        help="EOS table")
    parser.add_argument("--ex-temp", dest="do_temp_extend", action="store_true",
        help="Extend temperature?")
    parser.add_argument("--orig-min-logdens", dest="og_rhomin", type=float, default=999,
        help="Set minimum density of original table in log10(g cm^-3)")
    parser.add_argument("--orig-min-logtemp", dest="og_Tmin", type=float, default=999,
        help="Set minimum temperature of original table in log10(MeV)")
    parser.add_argument("--min-logtemp", dest="Tmin", required=True, type=float,
        help="New minimum temperature")
    parser.add_argument("--min-dens", dest="rhomin", required=True, type=float,
        help="New minimum density")
    parser.add_argument("--stitch-dens", dest="rho_stitch", required=True, type=float,
        help="Transition density")
    parser.add_argument("--stitch-temp", dest="T_stitch", required=True, type=float,
        help="Transition temperature")
    parser.add_argument("--width", dest="width", required=True, type=float,
        help="Width of transition region")
    parser.add_argument("--tol", dest="tol", required=True, type=float,
        help="Tolerance of transition region")

    args = parser.parse_args()
    main(args)
    print("Done :)")
