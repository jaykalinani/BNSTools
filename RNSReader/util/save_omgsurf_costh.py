import numpy as np

# Set options
rhoatmo = 1e-11
thermopath = "/home/acwen/Software/HydroRNS_DiffRot/Hydro_RNSID_DiffRot/Code/Run/IGL/solutions/solution_thermoTOT.dat"
omgpath = "/home/acwen/Software/HydroRNS_DiffRot/Hydro_RNSID_DiffRot/Code/Run/IGL/solutions/solution_OmegaTOT.dat"
savename = "omgsurf_atmo1e-11_IGL_jc.txt"

def main():
    thermofile = np.loadtxt(thermopath)
    omgfile = np.loadtxt(omgpath)

    cosax = np.unique(omgfile[:,1])
    out_costh_omg = np.empty((len(cosax), 2))

    for ii, costh in enumerate(cosax):
        rho_cos = thermofile[thermofile[:,1] == costh][:,5]
        r_cos = thermofile[thermofile[:,1] == costh][:,0]
        rsurf_idx = np.argmin(np.abs(rho_cos - rhoatmo))
        r_surf = r_cos[rsurf_idx]
        select_pt = np.logical_and(omgfile[:,0] == r_surf, omgfile[:,1] == costh)
        omgsurf = omgfile[select_pt][0, 2] # The grid is printed twice for some reason, if there are errors in the future check here
        out_costh_omg[ii, 0] = costh
        out_costh_omg[ii, 1] = omgsurf

    np.savetxt(savename, out_costh_omg)
    return

main()
