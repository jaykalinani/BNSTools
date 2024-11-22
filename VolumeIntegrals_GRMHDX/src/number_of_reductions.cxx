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
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE static inline int
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

  printf("VolumeIntegrals ERROR: You forgot to specify the number of "
         "reductions for Integration_quantity_keyword=%s!\n",
         Integration_quantity_keyword[which_integral]);
  exit(1);
  return -1000;
}
#endif
