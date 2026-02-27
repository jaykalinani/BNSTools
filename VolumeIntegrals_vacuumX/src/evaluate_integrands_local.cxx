#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

#include <loop_device.hxx>

#include <array>
#include <cassert>
#include <cmath>
#include <cstring>

#include "integrands.hxx"

using namespace Loop;

namespace {

enum class VI_vacuumX_centering_kind : int {
  vertex = 0,
  cell = 1,
};

struct VI_vacuumX_dynamic_gf {
  GF3D2<const CCTK_REAL> gf;
  VI_vacuumX_centering_kind centering;
};


inline int VI_vacuumX_varindex_or_error(const char *varname) {
  const int vi = CCTK_VarIndex(varname);
  if (vi < 0) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Couldn't find variable '%s'", varname);
  }
  return vi;
}

inline std::array<int, Loop::dim> VI_vacuumX_get_group_indextype(const int gi) {
  assert(gi >= 0);

  const int tags = CCTK_GroupTagsTableI(gi);
  assert(tags >= 0);

  std::array<CCTK_INT, Loop::dim> index;
  int iret = Util_TableGetIntArray(tags, Loop::dim, index.data(), "index");

  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    const int centering = CCTK_GroupCenteringTableI(gi);
    assert(centering >= 0);
    iret = Util_TableGetIntArray(centering, Loop::dim, index.data(),
                                 "centering");
  }

  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    index = {0, 0, 0}; // default vertex-centered
  } else if (iret >= 0) {
    assert(iret == Loop::dim);
  } else {
    assert(0);
  }

  std::array<int, Loop::dim> indextype;
  for (int d = 0; d < Loop::dim; ++d)
    indextype[d] = index[d];

  return indextype;
}

inline VI_vacuumX_centering_kind
VI_vacuumX_get_centering_or_error(const std::array<int, Loop::dim> &indextype,
                                  const char *name_for_error) {
  bool is_vertex = true;
  bool is_cell = true;
  for (int d = 0; d < Loop::dim; ++d) {
    is_vertex = is_vertex && (indextype[d] == 0);
    is_cell = is_cell && (indextype[d] == 1);
  }

  if (is_vertex)
    return VI_vacuumX_centering_kind::vertex;
  if (is_cell)
    return VI_vacuumX_centering_kind::cell;

  CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
              "Variable '%s' has unsupported centering index={%d,%d,%d}. "
              "Supported centerings are {0,0,0} (vertex) or {1,1,1} (cell).",
              name_for_error, indextype[0], indextype[1], indextype[2]);
}

inline VI_vacuumX_dynamic_gf
VI_vacuumX_get_dynamic_gf(const cGH *cctkGH, const int tl, const int vi,
                          const char *name_for_error) {
  const auto *ptr =
      static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vi));
  if (!ptr) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Couldn't get data pointer of input variable '%s'",
                name_for_error);
  }

  const int gi = CCTK_GroupIndexFromVarI(vi);
  assert(gi >= 0);

  const auto indextype = VI_vacuumX_get_group_indextype(gi);
  const GF3D2layout layout(cctkGH, indextype);
  return VI_vacuumX_dynamic_gf{
      GF3D2<const CCTK_REAL>(layout, ptr),
      VI_vacuumX_get_centering_or_error(indextype, name_for_error)};
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
VI_vacuumX_sample_dynamic_gf(const VI_vacuumX_dynamic_gf &f,
                             const PointDesc &p) {
  if (f.centering == VI_vacuumX_centering_kind::vertex) {
    return VI_vacuumX_avg_v2c_at(f.gf, p, p.I);
  }
  return f.gf(p.I);
}

} // namespace

extern "C" void VI_vacuumX_ComputeIntegrand(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_VI_vacuumX_ComputeIntegrand;
  DECLARE_CCTK_PARAMETERS;

  const int which_integral =
      NumIntegrals - static_cast<int>(*IntegralCounter) + 1;

  if (which_integral < 1 || which_integral > NumIntegrals ||
      which_integral > 100) {
    CCTK_VERROR("Invalid integral index: which_integral=%d NumIntegrals=%d IntegralCounter=%d",
                which_integral, NumIntegrals,
                static_cast<int>(*IntegralCounter));
  }

  const CCTK_INT timelevel = 0;

  /* Note: Must extend this if/else statement if adding a new integrand! */
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                  "H_M_CnstraintsL2")) {

    const int H_vi = VI_vacuumX_varindex_or_error(HamiltonianVarString);
    const int MU0_vi = VI_vacuumX_varindex_or_error(Momentum0VarString);
    const int MU1_vi = VI_vacuumX_varindex_or_error(Momentum1VarString);
    const int MU2_vi = VI_vacuumX_varindex_or_error(Momentum2VarString);

    const auto H_gf =
        VI_vacuumX_get_dynamic_gf(cctkGH, timelevel, H_vi, HamiltonianVarString);
    const auto MU0_gf =
        VI_vacuumX_get_dynamic_gf(cctkGH, timelevel, MU0_vi, Momentum0VarString);
    const auto MU1_gf =
        VI_vacuumX_get_dynamic_gf(cctkGH, timelevel, MU1_vi, Momentum1VarString);
    const auto MU2_gf =
        VI_vacuumX_get_dynamic_gf(cctkGH, timelevel, MU2_vi, Momentum2VarString);

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const CCTK_REAL H = VI_vacuumX_sample_dynamic_gf(H_gf, p);
          const CCTK_REAL M0 = VI_vacuumX_sample_dynamic_gf(MU0_gf, p);
          const CCTK_REAL M1 = VI_vacuumX_sample_dynamic_gf(MU1_gf, p);
          const CCTK_REAL M2 = VI_vacuumX_sample_dynamic_gf(MU2_gf, p);
          VolIntegrand1(p.I) = H * H;
          VolIntegrand2(p.I) = M0 * M0;
          VolIntegrand3(p.I) = M1 * M1;
          VolIntegrand4(p.I) = M2 * M2;
        });

  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "H_M2_CnstraintsL2")) {

    const int H_vi = VI_vacuumX_varindex_or_error(HamiltonianVarString);
    const auto H_gf =
        VI_vacuumX_get_dynamic_gf(cctkGH, timelevel, H_vi, HamiltonianVarString);

    // Z4c and CanudaX_BSSNMoL do not provide a default M^2 gridfunction.
    // If MomentumSquaredVarString is "UNSET", derive M^2 from momentum components.
    if (std::strcmp(MomentumSquaredVarString, "UNSET") == 0) {
      const int MU0_vi = VI_vacuumX_varindex_or_error(Momentum0VarString);
      const int MU1_vi = VI_vacuumX_varindex_or_error(Momentum1VarString);
      const int MU2_vi = VI_vacuumX_varindex_or_error(Momentum2VarString);
      const auto MU0_gf = VI_vacuumX_get_dynamic_gf(cctkGH, timelevel, MU0_vi,
                                                    Momentum0VarString);
      const auto MU1_gf = VI_vacuumX_get_dynamic_gf(cctkGH, timelevel, MU1_vi,
                                                    Momentum1VarString);
      const auto MU2_gf = VI_vacuumX_get_dynamic_gf(cctkGH, timelevel, MU2_vi,
                                                    Momentum2VarString);

      grid.loop_all_device<1, 1, 1>(
          grid.nghostzones,
          [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            const CCTK_REAL H = VI_vacuumX_sample_dynamic_gf(H_gf, p);
            const CCTK_REAL M0 = VI_vacuumX_sample_dynamic_gf(MU0_gf, p);
            const CCTK_REAL M1 = VI_vacuumX_sample_dynamic_gf(MU1_gf, p);
            const CCTK_REAL M2 = VI_vacuumX_sample_dynamic_gf(MU2_gf, p);
            const CCTK_REAL M2sq = M0 * M0 + M1 * M1 + M2 * M2;
            VolIntegrand1(p.I) = H * H;
            VolIntegrand2(p.I) = M2sq * M2sq;
          });
    } else {
      const int M2_vi = VI_vacuumX_varindex_or_error(MomentumSquaredVarString);
      const auto M2_gf = VI_vacuumX_get_dynamic_gf(cctkGH, timelevel, M2_vi,
                                                   MomentumSquaredVarString);

      grid.loop_all_device<1, 1, 1>(
          grid.nghostzones,
          [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            const CCTK_REAL H = VI_vacuumX_sample_dynamic_gf(H_gf, p);
            const CCTK_REAL M2 = VI_vacuumX_sample_dynamic_gf(M2_gf, p);
            VolIntegrand1(p.I) = H * H;
            VolIntegrand2(p.I) = M2 * M2;
          });
    }
  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "centeroflapse")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          VI_vacuumX_CoL_integrand(VolIntegrand1, VolIntegrand2, VolIntegrand3,
                                   VolIntegrand4, p, alp);
        });

  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "usepreviousintegrands")) {

    /* Do Nothing; the action for this is below. */

  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral], "one")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { VolIntegrand1(p.I) = 1.0; });

  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "ADM_Mass")) {

    const CCTK_REAL idx = 1.0 / CCTK_DELTA_SPACE(0);
    const CCTK_REAL idy = 1.0 / CCTK_DELTA_SPACE(1);
    const CCTK_REAL idz = 1.0 / CCTK_DELTA_SPACE(2);

    grid.loop_allmn_device<1, 1, 1>(
        grid.nghostzones, 1,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          VI_vacuumX_ADM_Mass_integrand_eval_derivs(
              VolIntegrand2, VolIntegrand3, VolIntegrand4, p, idx, idy, idz,
              alp, gxx, gxy, gxz, gyy, gyz, gzz);
        });

    grid.loop_allmn_device<1, 1, 1>(
        grid.nghostzones, 2,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          VI_vacuumX_ADM_Mass_integrand(VolIntegrand1, p, idx, idy, idz,
                                        VolIntegrand2, VolIntegrand3,
                                        VolIntegrand4);
        });

  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "ADM_Momentum")) {

    const CCTK_REAL idx = 1.0 / CCTK_DELTA_SPACE(0);
    const CCTK_REAL idy = 1.0 / CCTK_DELTA_SPACE(1);
    const CCTK_REAL idz = 1.0 / CCTK_DELTA_SPACE(2);

    grid.loop_allmn_device<1, 1, 1>(
        grid.nghostzones, 1,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          VI_vacuumX_ADM_Momentum_integrand_eval_derivs(
              VolIntegrand2, VolIntegrand3, VolIntegrand4, p, idx, idy, idz,
              alp, gxx, gxy, gxz, gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz,
              kzz);
        });

    grid.loop_allmn_device<1, 1, 1>(
        grid.nghostzones, 2,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          VI_vacuumX_ADM_Momentum_integrand(VolIntegrand1, p, idx, idy, idz,
                                            VolIntegrand2, VolIntegrand3,
                                            VolIntegrand4);
        });

  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "ADM_Angular_Momentum")) {

    const CCTK_REAL idx = 1.0 / CCTK_DELTA_SPACE(0);
    const CCTK_REAL idy = 1.0 / CCTK_DELTA_SPACE(1);
    const CCTK_REAL idz = 1.0 / CCTK_DELTA_SPACE(2);

    grid.loop_allmn_device<1, 1, 1>(
        grid.nghostzones, 1,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          VI_vacuumX_ADM_Angular_Momentum_integrand_eval_derivs(
              VolIntegrand2, VolIntegrand3, VolIntegrand4, p, idx, idy, idz,
              alp, gxx, gxy, gxz, gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz,
              kzz);
        });

    grid.loop_allmn_device<1, 1, 1>(
        grid.nghostzones, 2,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          VI_vacuumX_ADM_Angular_Momentum_integrand(
              VolIntegrand1, p, idx, idy, idz, VolIntegrand2, VolIntegrand3,
              VolIntegrand4);
        });

  } else {

    /* Print a warning if no integrand is computed because
     * Integration_quantity_keyword unrecognized. */
    printf("VolumeIntegrals_vacuumX: WARNING: Integrand not computed. Did not "
           "understand Integration_quantity_keyword[%d] = %s\n",
           which_integral, Integration_quantity_keyword[which_integral]);
    // Unknown keyword: clear integrands so stale values cannot leak into sums.
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          VolIntegrand1(p.I) = 0.0;
          VolIntegrand2(p.I) = 0.0;
          VolIntegrand3(p.I) = 0.0;
          VolIntegrand4(p.I) = 0.0;
        });
  }

  if (cctk_iteration == 0) {
    volintegral_inside_sphere__center_x[which_integral] =
        volintegral_sphere__center_x_initial[which_integral];
    volintegral_inside_sphere__center_y[which_integral] =
        volintegral_sphere__center_y_initial[which_integral];
    volintegral_inside_sphere__center_z[which_integral] =
        volintegral_sphere__center_z_initial[which_integral];

    volintegral_outside_sphere__center_x[which_integral] =
        volintegral_sphere__center_x_initial[which_integral];
    volintegral_outside_sphere__center_y[which_integral] =
        volintegral_sphere__center_y_initial[which_integral];
    volintegral_outside_sphere__center_z[which_integral] =
        volintegral_sphere__center_z_initial[which_integral];
  }

  if (volintegral_sphere__tracks__amr_centre[which_integral] != -1) {
    const int which_centre =
        volintegral_sphere__tracks__amr_centre[which_integral];
    if (which_centre < 0 || which_centre > 2) {
      CCTK_VERROR("Invalid BoxInBox centre index %d for integral %d; valid range is [0,2]",
                  which_centre, which_integral);
    }

    volintegral_inside_sphere__center_x[which_integral] = position_x[which_centre];
    volintegral_inside_sphere__center_y[which_integral] = position_y[which_centre];
    volintegral_inside_sphere__center_z[which_integral] = position_z[which_centre];

    volintegral_outside_sphere__center_x[which_integral] =
        position_x[which_centre];
    volintegral_outside_sphere__center_y[which_integral] =
        position_y[which_centre];
    volintegral_outside_sphere__center_z[which_integral] =
        position_z[which_centre];
  }

  /* ZERO OUT INTEGRATION REGIONS */

  /* Set integrands to zero outside a sphere centered at x,y,z.
     I.e., this results in the integral being restricted INSIDE sphere */
  if (volintegral_inside_sphere__radius[which_integral] > 0.0) {
    const double radius = volintegral_inside_sphere__radius[which_integral];
    double xprime = volintegral_sphere__center_x_initial[which_integral];
    double yprime = volintegral_sphere__center_y_initial[which_integral];
    double zprime = volintegral_sphere__center_z_initial[which_integral];

    if (cctk_iteration > 0 &&
        (amr_centre__tracks__volintegral_inside_sphere[which_integral] != -1 ||
         volintegral_sphere__tracks__amr_centre[which_integral] != -1)) {
      xprime = volintegral_inside_sphere__center_x[which_integral];
      yprime = volintegral_inside_sphere__center_y[which_integral];
      zprime = volintegral_inside_sphere__center_z[which_integral];
    }

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const double dx = p.x - xprime;
          const double dy = p.y - yprime;
          const double dz = p.z - zprime;
          if (sqrt(dx * dx + dy * dy + dz * dz) > radius) {
            VolIntegrand1(p.I) = 0.0;
            VolIntegrand2(p.I) = 0.0;
            VolIntegrand3(p.I) = 0.0;
            VolIntegrand4(p.I) = 0.0;
          }
        });
  }

  /* Set integrands to zero inside a sphere centered at x,y,z.
     I.e., this results in the integral being restricted OUTSIDE sphere. */
  if (volintegral_outside_sphere__radius[which_integral] > 0.0) {
    const double radius = volintegral_outside_sphere__radius[which_integral];
    const double xprime = volintegral_outside_sphere__center_x[which_integral];
    const double yprime = volintegral_outside_sphere__center_y[which_integral];
    const double zprime = volintegral_outside_sphere__center_z[which_integral];

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const double dx = p.x - xprime;
          const double dy = p.y - yprime;
          const double dz = p.z - zprime;
          if (sqrt(dx * dx + dy * dy + dz * dz) <= radius) {
            VolIntegrand1(p.I) = 0.0;
            VolIntegrand2(p.I) = 0.0;
            VolIntegrand3(p.I) = 0.0;
            VolIntegrand4(p.I) = 0.0;
          }
        });
  }
}
