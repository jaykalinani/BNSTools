
/* Load the TOV_interpolate_1D() function */
#include "tov_interp.h"

// Declare initial data input struct:
//          stores data from initial data solver,
//          so they can be put on the numerical grid.
typedef struct __ID_inputs {
    CCTK_REAL Rbar;
    int Rbar_idx;
    int interp_stencil_size;
    int numlines_in_file;
    CCTK_REAL *r_Schw_arr,*rho_arr,*rho_baryon_arr,*P_arr,*M_arr,*expnu_arr,*exp4phi_arr,*rbar_arr;
} ID_inputs;

/* This function returns the TOV quantities at point rr
 * by interpolating the data in the TOV initial data file.
 */
void interpolate_TOV_solution_to_point(const CCTK_REAL rr, ID_inputs other_inputs,
                           CCTK_REAL *exp_4phi, CCTK_REAL *expnu,
                           CCTK_REAL *Pressure, CCTK_REAL *rho_baryon, CCTK_REAL *rho__total_energy_density) {

  /* The mass valus is not used, but we have to
   * store it in this dummy variable because the
   * initial data file contains it.
   */
  CCTK_REAL M;

  /* Perform the interpolation, returning:
   *  - rho__total_energy_density
   *  - rho_baryon
   *  - Pressure
   *  - Mass (dummy variable, unused)
   *  - exp(nu)
   *  - exp(4phi)
   */
  TOV_interpolate_1D(rr,other_inputs.Rbar,other_inputs.Rbar_idx,other_inputs.interp_stencil_size,
                     other_inputs.numlines_in_file,
                     other_inputs.r_Schw_arr,other_inputs.rho_arr,other_inputs.rho_baryon_arr,other_inputs.P_arr,other_inputs.M_arr,
                     other_inputs.expnu_arr,other_inputs.exp4phi_arr,other_inputs.rbar_arr,
                     rho__total_energy_density,rho_baryon,Pressure,&M,expnu,exp_4phi);

}
