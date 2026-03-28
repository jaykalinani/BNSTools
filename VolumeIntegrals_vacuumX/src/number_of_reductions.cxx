#ifndef VI_VACUUMX_REDUCTIONS_CXX
#define VI_VACUUMX_REDUCTIONS_CXX

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

#include <cassert>
#include <cmath>
#include <cstdlib>

/* ADD TO THIS LIST IF YOU HAVE NEW INTEGRAND */
CCTK_ATTRIBUTE_ALWAYS_INLINE static inline int
VI_vacuumX_number_of_reductions(int which_integral) {
  DECLARE_CCTK_PARAMETERS;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                  "H_M_CnstraintsL2"))
    return 4;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                  "H_M2_CnstraintsL2"))
    return 2;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                  "usepreviousintegrands"))
    return volintegral_usepreviousintegrands_num_integrands[which_integral];
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                  "centeroflapse"))
    return 4;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral], "one"))
    return 1;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral], "ADM_Mass"))
    return 1;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                  "ADM_Momentum"))
    return 1;
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                  "ADM_Angular_Momentum"))
    return 1;

  CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
              "VolumeIntegrals_vacuumX ERROR: You forgot to specify the number of reductions for Integration_quantity_keyword=%s!\n",
              Integration_quantity_keyword[which_integral]);

  return -1000;
}

#endif
