#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <array>
#include <vector>

#ifdef CCTK_MPI
#include <mpi.h>
#endif

#include "number_of_reductions.cxx"
#include "../../../CarpetX/CarpetX/src/driver.hxx"
#include "../../../CarpetX/CarpetX/src/reduction.hxx"

namespace {

bool VI_vacuumX_is_surface_integral(const char *keyword) {
  return CCTK_Equals(keyword, "ADM_Mass_Surface") ||
         CCTK_Equals(keyword, "ADM_Momentum_Surface") ||
         CCTK_Equals(keyword, "ADM_Angular_Momentum_Surface");
}

CCTK_REAL VI_vacuumX_simpson_weight(const int i, const int n) {
  if (i == 0 || i == n)
    return 1.0;
  return i % 2 == 0 ? 2.0 : 4.0;
}

void VI_vacuumX_direct_surface_integral(
    const cGH *cctkGH, const int num_reductions, const CCTK_REAL center_x,
    const CCTK_REAL center_y, const CCTK_REAL center_z,
    const CCTK_REAL radius, const int ntheta, const int nphi,
    const char *interpolator_name, const char *interpolator_pars,
    CCTK_REAL result[4]) {
  if (radius <= 0.0)
    CCTK_VERROR("Direct spherical surface integration requires radius > 0");
  if (ntheta < 3 || nphi < 2 || ntheta % 2 == 0)
    CCTK_VERROR("Direct spherical Simpson integration requires odd "
                "surface_ntheta >= 3 and surface_nphi >= 2; found %d and %d",
                ntheta, nphi);
  if (num_reductions != 1 && num_reductions != 3)
    CCTK_VERROR("Direct spherical surface integration supports 1 or 3 "
                "components; found %d",
                num_reductions);

  const bool is_root = CCTK_MyProc(cctkGH) == 0;
  const int npoints = is_root ? ntheta * nphi : 0;
  std::vector<CCTK_REAL> x(npoints), y(npoints), z(npoints);
  std::array<std::vector<CCTK_REAL>, 4> values;
  for (int n = 0; n < num_reductions; ++n)
    values[n].resize(npoints);

  const CCTK_REAL dtheta = M_PI / (ntheta - 1);
  const CCTK_REAL dphi = 2.0 * M_PI / nphi;
  if (is_root) {
    for (int iphi = 0; iphi < nphi; ++iphi) {
      const CCTK_REAL phi = iphi * dphi;
      for (int itheta = 0; itheta < ntheta; ++itheta) {
        const CCTK_REAL theta = itheta * dtheta;
        const int p = itheta + ntheta * iphi;
        const CCTK_REAL sintheta = sin(theta);
        x[p] = center_x + radius * sintheta * cos(phi);
        y[p] = center_y + radius * sintheta * sin(phi);
        z[p] = center_z + radius * cos(theta);
      }
    }
  }

  const void *interp_coords[3] = {x.data(), y.data(), z.data()};
  const CCTK_INT input_indices[4] = {
      CCTK_VarIndex("VolumeIntegrals_vacuumX::VolIntegrand1"),
      CCTK_VarIndex("VolumeIntegrals_vacuumX::VolIntegrand2"),
      CCTK_VarIndex("VolumeIntegrals_vacuumX::VolIntegrand3"),
      CCTK_VarIndex("VolumeIntegrals_vacuumX::VolIntegrand4")};
  const CCTK_INT output_types[4] = {0, 0, 0, 0};
  CCTK_POINTER output_arrays[4] = {
      values[0].data(), values[1].data(), values[2].data(), values[3].data()};

  const int interp_handle = CCTK_InterpHandle(interpolator_name);
  if (interp_handle < 0)
    CCTK_VERROR("Could not obtain interpolator handle for '%s': %d",
                interpolator_name, interp_handle);

  const int table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  if (table_handle < 0)
    CCTK_VERROR("Could not create interpolation parameter table: %d",
                table_handle);
  const int table_status =
      Util_TableSetFromString(table_handle, interpolator_pars);
  if (table_status < 0)
    CCTK_VERROR("Could not parse surface_interpolator_pars '%s': %d",
                interpolator_pars, table_status);

  const int interp_status =
      DriverInterpolate(cctkGH, 3, interp_handle, table_handle, 0, npoints, 0,
                        interp_coords, num_reductions, input_indices,
                        num_reductions, output_types, output_arrays);
  const int destroy_status = Util_TableDestroy(table_handle);
  if (destroy_status < 0)
    CCTK_VERROR("Could not destroy interpolation parameter table: %d",
                destroy_status);
  if (interp_status < 0)
    CCTK_VERROR("Direct spherical interpolation failed with status %d",
                interp_status);

  CCTK_REAL local_result[4] = {0.0, 0.0, 0.0, 0.0};
  if (is_root) {
    const CCTK_REAL area_scale = radius * radius * dtheta * dphi / 3.0;
    for (int n = 0; n < num_reductions; ++n) {
      CCTK_REAL sum = 0.0;
      for (int iphi = 0; iphi < nphi; ++iphi) {
        for (int itheta = 0; itheta < ntheta; ++itheta) {
          const CCTK_REAL theta = itheta * dtheta;
          const CCTK_REAL wtheta =
              VI_vacuumX_simpson_weight(itheta, ntheta - 1);
          const int p = itheta + ntheta * iphi;
          if (!std::isfinite(values[n][p]))
            CCTK_VERROR("Non-finite direct-sphere interpolation result in "
                        "component %d at theta index %d, phi index %d",
                        n, itheta, iphi);
          sum += wtheta * values[n][p] * sin(theta);
        }
      }
      local_result[n] = area_scale * sum;
    }
  }

  for (int n = 0; n < num_reductions; ++n)
    result[n] = local_result[n];

#ifdef CCTK_MPI
  MPI_Comm comm = MPI_COMM_WORLD;
  if (CCTK_IsFunctionAliased("GetMPICommWorld"))
    comm = *static_cast<const MPI_Comm *>(GetMPICommWorld(cctkGH));

  MPI_Datatype mpi_real;
  if (sizeof(CCTK_REAL) == sizeof(float))
    mpi_real = MPI_FLOAT;
  else if (sizeof(CCTK_REAL) == sizeof(double))
    mpi_real = MPI_DOUBLE;
  else if (sizeof(CCTK_REAL) == sizeof(long double))
    mpi_real = MPI_LONG_DOUBLE;
  else
    CCTK_ERROR("Unsupported CCTK_REAL type for MPI broadcast");

  const int broadcast_status =
      MPI_Bcast(result, num_reductions, mpi_real, 0, comm);
  if (broadcast_status != MPI_SUCCESS)
    CCTK_VERROR("Could not broadcast direct spherical integral: %d",
                broadcast_status);
#endif
}

} // namespace

extern "C" void VI_vacuumX_DoSum(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_VI_vacuumX_DoSum;
  DECLARE_CCTK_PARAMETERS;
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
    printf("VolumeIntegrals_vacuumX: Iter %d, num_reductions=%d, Integ. quantity=%s, sphere moves/tracks AMR centre=%d/%d | SURFACE r=%e ntheta=%d nphi=%d interpolator=%s | INSIDE center x,y,z=%e,%e,%e ; r=%e | OUTSIDE center x,y,z=%e,%e,%e ; r=%e\n",
           which_integral, num_reductions,
           Integration_quantity_keyword[which_integral],
           amr_centre__tracks__volintegral_inside_sphere[which_integral],
           volintegral_sphere__tracks__amr_centre[which_integral],
           volintegral_surface_sphere__radius[which_integral],
           surface_ntheta, surface_nphi, surface_interpolator_name,
           volintegral_inside_sphere__center_x[which_integral],
           volintegral_inside_sphere__center_y[which_integral],
           volintegral_inside_sphere__center_z[which_integral],
           volintegral_inside_sphere__radius[which_integral],
           volintegral_outside_sphere__center_x[which_integral],
           volintegral_outside_sphere__center_y[which_integral],
           volintegral_outside_sphere__center_z[which_integral],
           volintegral_outside_sphere__radius[which_integral]);

  if (VI_vacuumX_is_surface_integral(
          Integration_quantity_keyword[which_integral])) {
    if (!CCTK_EQUALS(domain, "full"))
      CCTK_VERROR("Direct spherical surface integration currently requires "
                  "domain='full'");

    CCTK_REAL result[4] = {0.0, 0.0, 0.0, 0.0};
    VI_vacuumX_direct_surface_integral(
        cctkGH, num_reductions,
        volintegral_inside_sphere__center_x[which_integral],
        volintegral_inside_sphere__center_y[which_integral],
        volintegral_inside_sphere__center_z[which_integral],
        volintegral_surface_sphere__radius[which_integral], surface_ntheta,
        surface_nphi, surface_interpolator_name, surface_interpolator_pars,
        result);
    for (int i = 0; i < num_reductions; ++i)
      VolIntegral[4 * which_integral + i] = result[i];
    return;
  }

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
    if (!std::isfinite(redsum)) {
      CCTK_VERROR("Non-finite reduction for integral %d component %d/%d "
                  "keyword='%s': sum=%e min=%e max=%e volume=%e",
                  which_integral, i + 1, num_reductions,
                  Integration_quantity_keyword[which_integral], redsum,
                  red.min, red.max, red.vol);
    }

    const bool is_constraint_norm =
        CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                    "H_M_CnstraintsL2") ||
        CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                    "H_M2_CnstraintsL2");
    if (is_constraint_norm) {
      if (redsum < 0.0)
        CCTK_VERROR("Negative squared-constraint integral %e for integral "
                    "%d component %d",
                    redsum, which_integral, i);
      VolIntegral[4 * which_integral + i] = sqrt(redsum);
    } else {
      VolIntegral[4 * which_integral + i] = redsum;
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
      if (which_centre < 0 || which_centre >= 100) {
        CCTK_VERROR("Invalid BoxInBox centre index %d for integral %d; valid range is [0,99]",
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
