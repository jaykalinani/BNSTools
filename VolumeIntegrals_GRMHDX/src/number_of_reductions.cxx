#ifndef VI_GRMHDX_REDUCTIONS_CXX
#define VI_GRMHDX_REDUCTIONS_CXX

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "util_Table.h"
#include <cassert>
#include <cmath>
#include <cstdlib>

#include <loop_device.hxx>

/* ADD TO THIS LIST IF YOU HAVE NEW INTEGRAND */
CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE static inline int
VI_GRMHDX_number_of_reductions(int which_integral) {
  DECLARE_CCTK_PARAMETERS;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                  "usepreviousintegrands"))
    return volintegral_usepreviousintegrands_num_integrands[which_integral];
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral], "centerofmass"))
    return 4;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral], "restmass"))
    return 1;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral], "one"))
    return 1;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                  "density_weighted_norm_B_field"))
    return 2;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                  "volume_average_norm_B_ab"))
    return 4;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                  "volume_average_norm_B_cd"))
    return 4;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                  "kinetic_energy"))
    return 3;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                  "kinetic_energy_palenzuela"))
    return 3;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                  "kinetic_energy_shibata"))
    return 3;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                  "thermal_energy"))
    return 3;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                  "magnetic_energy_total"))
    return 2;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                  "magnetic_energy_comov"))
    return 1;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                  "kinetic_energy_total"))
    return 1;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral], "em_energy_ab"))
    return 2;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral], "em_energy_cd"))
    return 2;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral], "coordvolume"))
    return 1;

  printf("VolumeIntegrals ERROR: You forgot to specify the number of "
         "reductions for Integration_quantity_keyword=%s!\n",
         Integration_quantity_keyword[which_integral]);
  exit(1);
  return -1000;
}
#endif
