#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "util_Table.h"
#include <cassert>
#include <cmath>
#include <cstdlib>

extern "C" void VI_GRMHDX_InitializeIntegralCounterToZero(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_VI_GRMHDX_InitializeIntegralCounterToZero;
  DECLARE_CCTK_PARAMETERS;

  // Init Counter
  *IntegralCounter = 0;

  // Init CoM Coord Holders
  *comx = 0;
  *comy = 0;
  *comz = 0;

  // Init Array that Holds Results
  for (int ii = 0; ii < 101; ii++) {
    for (int jj = 0; jj < 4; jj++) {
      const int arridx = 4 * ii + jj;
      VolIntegral[arridx] = 0.0;
    }

    // Seed moving-sphere centers so analysis/file-output inputs are valid at
    // iteration 0 and after recovery.
    volintegral_inside_sphere__center_x[ii] =
        volintegral_sphere__center_x_initial[ii];
    volintegral_inside_sphere__center_y[ii] =
        volintegral_sphere__center_y_initial[ii];
    volintegral_inside_sphere__center_z[ii] =
        volintegral_sphere__center_z_initial[ii];

    volintegral_outside_sphere__center_x[ii] =
        volintegral_sphere__center_x_initial[ii];
    volintegral_outside_sphere__center_y[ii] =
        volintegral_sphere__center_y_initial[ii];
    volintegral_outside_sphere__center_z[ii] =
        volintegral_sphere__center_z_initial[ii];
  }

  if (verbose == 2)
    printf("VolumeIntegrals_GRMHDX: Just set IntegralCounter to %d\n",
           static_cast<int>(*IntegralCounter));
}

extern "C" void VI_GRMHDX_InitializeIntegralCounter(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_VI_GRMHDX_InitializeIntegralCounter;
  DECLARE_CCTK_PARAMETERS;

  if (VolIntegral_out_every <= 0 || NumIntegrals <= 0) {
    *IntegralCounter = 0;
    return;
  }

  if (cctk_iteration % VolIntegral_out_every == 0) {
    *IntegralCounter = NumIntegrals;
    if (verbose == 2)
      printf("VolumeIntegrals_GRMHDX: Just set IntegralCounter to %d == "
             "NumIntegrals\n",
             static_cast<int>(*IntegralCounter));
  } else {
    *IntegralCounter = 0;
  }
}

extern "C" void VI_GRMHDX_DecrementIntegralCounter(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_VI_GRMHDX_DecrementIntegralCounter;
  DECLARE_CCTK_PARAMETERS;

  (*IntegralCounter)--;
  if (verbose == 2)
    printf("VolumeIntegrals_GRMHDX: Just decremented IntegralCounter to %d\n",
           static_cast<int>(*IntegralCounter));
}
