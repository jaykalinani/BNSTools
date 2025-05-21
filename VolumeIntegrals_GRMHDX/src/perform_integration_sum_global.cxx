#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "util_Table.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "number_of_reductions.cxx"
#include "../../../CarpetX/CarpetX/src/driver.hxx"
#include "../../../CarpetX/CarpetX/src/reduction.hxx"

extern "C" void VI_GRMHDX_DoSum(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_VI_GRMHDX_DoSum;
  DECLARE_CCTK_PARAMETERS;

  int which_integral = NumIntegrals - *IntegralCounter + 1;

  /* FIXME: Add this symmetry stuff... Should be straightforward. */
  CCTK_REAL sym_factor1, sym_factor2, sym_factor3;

  /*
    if (CCTK_EQUALS(domain,"bitant")){
      sym_factor1 = 2.0e0;
      sym_factor2 = 2.0e0;
      sym_factor3 = 0.0e0;
    } else if (CCTK_EQUALS(domain,"octant")){
      sym_factor1 = 8.0e0;
      sym_factor2 = 0.0e0;
      sym_factor3 = 0.0e0;
    } else {
      sym_factor1 = 1.0e0;
      sym_factor2 = 1.0e0;
      sym_factor3 = 1.0e0;
    }
  */

  sym_factor1 = 1.0e0;
  sym_factor2 = 1.0e0;
  sym_factor3 = 1.0e0;

  // CarpetX reductions currently include the cell volume in the sum. Keeping
  // this code here in case that ever changes. const amrex::Geometry &geom =
  // CarpetX::ghext->patchdata.at(0).amrcore->Geom(0); const CCTK_REAL *restrict
  // const dx = geom.CellSize(); double d3x = dx[0] * dx[1] * dx[2];

  /* Note: Must edit VI_GRMHDX_number_of_reductions() when adding new
     integrands! This function is defined in VI_GRMHDX_number_of_reductions.C */
  int num_reductions = VI_GRMHDX_number_of_reductions(which_integral);

  if (verbose >= 1)
    printf("VolumeIntegrals_GRMHDX: Iter %d, num_reductions=%d, Integ. "
           "quantity=%s, sphere moves/tracks AMR centre=%d/%d | INSIDE center "
           "x,y,z=%e,%e,%e ; r=%e | OUTSIDE center x,y,z=%e,%e,%e ; r=%e\n",
           which_integral, num_reductions,
           Integration_quantity_keyword[which_integral],
           amr_centre__tracks__volintegral_inside_sphere[which_integral],
           volintegral_sphere__tracks__amr_centre[which_integral],
           volintegral_inside_sphere__center_x[which_integral],
           volintegral_inside_sphere__center_y[which_integral],
           volintegral_inside_sphere__center_z[which_integral],
           volintegral_inside_sphere__radius[which_integral],
           volintegral_outside_sphere__center_x[which_integral],
           volintegral_outside_sphere__center_y[which_integral],
           volintegral_outside_sphere__center_z[which_integral],
           volintegral_outside_sphere__radius[which_integral]);

  /* Perform the reduction sums across all MPI processes */

  for (int i = 0; i < num_reductions; i++) {
    char integralname[100];
    sprintf(integralname, "VolumeIntegrals_GRMHDX::VolIntegrand%d", i + 1);

    // Get group and var idx for reduction
    int varindex = CCTK_VarIndex(integralname);
    const int gi = CCTK_GroupIndexFromVarI(varindex);
    assert(gi >= 0);
    const int v0 = CCTK_FirstVarIndexI(gi);
    assert(v0 >= 0);
    const int vi = varindex - v0;

    // Perform reduction. We only want the summed integrands.
    const CarpetX::reduction<CCTK_REAL, 3> red = CarpetX::reduce(gi, vi, 0);

    VolIntegral[4 * (which_integral) + i] =
        red.sum; // * d3x; // <- Multiply the integrand by d3x

    if (verbose == 2)
      printf("VolumeIntegrals_GRMHDX: Iteration %d, reduction %d of %d. "
             "Reduction value=%e\n",
             which_integral, i + 1, num_reductions,
             VolIntegral[4 * (which_integral) + i]);
  }

  /* AMR box centre tracks volume integral output */
  if (num_reductions == 4 &&
      amr_centre__tracks__volintegral_inside_sphere[which_integral] != -1) {
    double norm = sym_factor1 * VolIntegral[4 * (which_integral) + 3];
    volintegral_inside_sphere__center_x[which_integral] =
        sym_factor2 * VolIntegral[4 * (which_integral) + 0] / norm;
    volintegral_inside_sphere__center_y[which_integral] =
        sym_factor2 * VolIntegral[4 * (which_integral) + 1] / norm;
    volintegral_inside_sphere__center_z[which_integral] =
        sym_factor3 * VolIntegral[4 * (which_integral) + 2] / norm;

    int which_centre =
        amr_centre__tracks__volintegral_inside_sphere[which_integral];

    if (verbose >= 1)
      printf("VolumeIntegrals_GRMHDX: AMR centre #%d tracks Integral %d: "
             "(x,y,z)=(%e,%e,%e) [norm=%e]. Prev centre @ (%e,%e,%e).\n",
             amr_centre__tracks__volintegral_inside_sphere[which_integral],
             which_integral,
             volintegral_inside_sphere__center_x[which_integral],
             volintegral_inside_sphere__center_y[which_integral],
             volintegral_inside_sphere__center_z[which_integral], norm,
             position_x[which_centre], position_y[which_centre],
             position_z[which_centre]);

    /* Activate AMR box tracking for this centre */
    active[which_centre] = 1;
    /* Update AMR box centre position.
       Note that this will have no effect until cctk_iteration%regrid_every==0
     */
    position_x[which_centre] =
        volintegral_inside_sphere__center_x[which_integral];
    position_y[which_centre] =
        volintegral_inside_sphere__center_y[which_integral];
    position_z[which_centre] =
        volintegral_inside_sphere__center_z[which_integral];
  } else {
    for (int i = 0; i < num_reductions; i++)
      VolIntegral[4 * (which_integral) + i] *= sym_factor1;
  }

  /* Set global CoM tracker */
  if (which_integral == 1 && CCTK_EQUALS(Integration_quantity_keyword[which_integral], "centerofmass")
        && set_origin_with_VIX ) {
    double norm = sym_factor1 * VolIntegral[4 * (which_integral) + 3];
    *comx = sym_factor2 * VolIntegral[4 * (which_integral) + 0] / norm;
    *comy = sym_factor2 * VolIntegral[4 * (which_integral) + 1] / norm;
    *comz = sym_factor3 * VolIntegral[4 * (which_integral) + 2] / norm;
  }
    
}
