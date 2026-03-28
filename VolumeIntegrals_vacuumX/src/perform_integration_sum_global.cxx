#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include "number_of_reductions.cxx"
#include "../../../CarpetX/CarpetX/src/driver.hxx"
#include "../../../CarpetX/CarpetX/src/reduction.hxx"

extern "C" void VI_vacuumX_DoSum(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_VI_vacuumX_DoSum;
  DECLARE_CCTK_PARAMETERS;
  const int myproc = CCTK_MyProc(cctkGH);

  const int which_integral =
      NumIntegrals - static_cast<int>(*IntegralCounter) + 1;
  if (which_integral < 1 || which_integral > NumIntegrals ||
      which_integral > 100) {
    CCTK_VERROR("Invalid integral index: which_integral=%d NumIntegrals=%d IntegralCounter=%d",
                which_integral, NumIntegrals,
                static_cast<int>(*IntegralCounter));
  }

  CCTK_REAL sym_factor1, sym_factor2, sym_factor3;

  if (CCTK_EQUALS(domain, "bitant")) {
    sym_factor1 = 2.0e0;
    sym_factor2 = 2.0e0;
    sym_factor3 = 0.0e0;
  } else if (CCTK_EQUALS(domain, "octant")) {
    sym_factor1 = 8.0e0;
    sym_factor2 = 0.0e0;
    sym_factor3 = 0.0e0;
  } else {
    sym_factor1 = 1.0e0;
    sym_factor2 = 1.0e0;
    sym_factor3 = 1.0e0;
  }

  const int num_reductions = VI_vacuumX_number_of_reductions(which_integral);

  if (verbose >= 1)
    printf("VolumeIntegrals_vacuumX: Iter %d, num_reductions=%d, Integ. quantity=%s, sphere moves/tracks AMR centre=%d/%d | INSIDE center x,y,z=%e,%e,%e ; r=%e | OUTSIDE center x,y,z=%e,%e,%e ; r=%e\n",
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

  for (int i = 0; i < num_reductions; i++) {
    char integralname[100];
    sprintf(integralname, "VolumeIntegrals_vacuumX::VolIntegrand%d", i + 1);

    const int varindex = CCTK_VarIndex(integralname);
    const int gi = CCTK_GroupIndexFromVarI(varindex);
    assert(gi >= 0);
    const int v0 = CCTK_FirstVarIndexI(gi);
    assert(v0 >= 0);
    const int vi = varindex - v0;

    const CarpetX::reduction<CCTK_REAL, 3> red = CarpetX::reduce(gi, vi, 0);
    const CCTK_REAL redsum = red.sum;
    if (std::isfinite(redsum)) {
      VolIntegral[4 * (which_integral) + i] = redsum;
    } else {
      VolIntegral[4 * (which_integral) + i] = 0.0;
      if (myproc == 0 && verbose >= 1) {
        printf("VolumeIntegrals_vacuumX: Replaced non-finite reduction with 0: integral=%d reduction=%d/%d keyword=%s sum=%e min=%e max=%e vol=%e\n",
               which_integral, i + 1, num_reductions,
               Integration_quantity_keyword[which_integral], redsum, red.min,
               red.max, red.vol);
      }
    }

    if (verbose == 2)
      printf("VolumeIntegrals_vacuumX: Iteration %d, reduction %d of %d. Reduction value=%e\n",
             which_integral, i + 1, num_reductions,
             VolIntegral[4 * (which_integral) + i]);
  }

  if (num_reductions == 4 &&
      amr_centre__tracks__volintegral_inside_sphere[which_integral] != -1) {
    const double norm = sym_factor1 * VolIntegral[4 * (which_integral) + 3];
    if (std::isfinite(norm) && std::abs(norm) > 0.0) {
      volintegral_inside_sphere__center_x[which_integral] =
          sym_factor2 * VolIntegral[4 * (which_integral) + 0] / norm;
      volintegral_inside_sphere__center_y[which_integral] =
          sym_factor2 * VolIntegral[4 * (which_integral) + 1] / norm;
      volintegral_inside_sphere__center_z[which_integral] =
          sym_factor3 * VolIntegral[4 * (which_integral) + 2] / norm;

      const int which_centre =
          amr_centre__tracks__volintegral_inside_sphere[which_integral];
      if (which_centre < 0 || which_centre > 2) {
        CCTK_VERROR("Invalid BoxInBox centre index %d for integral %d; valid range is [0,2]",
                    which_centre, which_integral);
      }

      if (verbose >= 1)
        printf("VolumeIntegrals_vacuumX: AMR centre #%d tracks Integral %d: (x,y,z)=(%e,%e,%e) [norm=%e]. Prev centre @ (%e,%e,%e).\n",
               which_centre, which_integral,
               volintegral_inside_sphere__center_x[which_integral],
               volintegral_inside_sphere__center_y[which_integral],
               volintegral_inside_sphere__center_z[which_integral], norm,
               position_x[which_centre], position_y[which_centre],
               position_z[which_centre]);

      active[which_centre] = 1;
      position_x[which_centre] =
          volintegral_inside_sphere__center_x[which_integral];
      position_y[which_centre] =
          volintegral_inside_sphere__center_y[which_integral];
      position_z[which_centre] =
          volintegral_inside_sphere__center_z[which_integral];
    }
  } else {
    for (int i = 0; i < num_reductions; i++)
      VolIntegral[4 * (which_integral) + i] *= sym_factor1;
  }
}
