#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "rnsreader_utils.hxx"
#include "consts.h"
#include "aster_utils.hxx"

#include <loop.hxx>
#include <loop_device.hxx>

namespace RNSReader {

using namespace Loop;
using namespace amrex;
using namespace std;
using namespace AsterUtils;

Omega_th_reader* vel_th_reader; // storage for vel_atmo(theta)

extern "C" void RNSReader_Init_VelAtmo(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_RNSReader_Init_VelAtmo;
  DECLARE_CCTK_PARAMETERS;

  // Open the Y_e file, which should countain Y_e(rho) for the EOS table slice
  FILE *vel_th_file = fopen(vel_th_filename, "r");
    
  // Check if everything is OK with the file
  if ((vel_th_file = fopen(vel_th_filename, "r")) == NULL) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
      "File \"%s\" does not exist. ABORTING", vel_th_filename);
  } else {
    // Set interpolation stencil size
    const int interp_stencil_size = 5;

    vel_th_reader =
      (Omega_th_reader *)The_Managed_Arena()->alloc(sizeof *vel_th_reader);
    vel_th_reader->init(vel_th_file, interp_stencil_size);

    // Close the file
    fclose(vel_th_file);
  }
}


extern "C" void RNSReader_Set_VelAtmo(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_RNSReader_Set_VelAtmo;
  DECLARE_CCTK_PARAMETERS;

  /* Gridfunctions */
  const vec<GF3D2<const CCTK_REAL>, dim> gf_beta{betax, betay, betaz};

  grid.loop_all<1, 1, 1>(grid.nghostzones,
                         [=] CCTK_HOST(const Loop::PointDesc &p)
                             CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const CCTK_REAL velxL = velx(p.I);
    const CCTK_REAL velyL = vely(p.I);
    const CCTK_REAL radial_distance = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
    
    // Grading rho
    CCTK_REAL rho_atm = (radial_distance > r_atmo)
      ? (rho_abs_min * pow((r_atmo / radial_distance), n_rho_atmo))
      : rho_abs_min;
    const CCTK_REAL rho_atmo_cut = rho_atm * (1 + atmo_tol);

    if (rho(p.I) <= rho_atmo_cut) {
      CCTK_REAL costh = std::abs(p.z) / radial_distance; // cos(theta) on RNS grid runs from 0 to 1
      CCTK_REAL omg_surf;
      vel_th_reader->interpolate_1d_quantity_as_function_of_th(
          MDIV - 1, costh, &omg_surf); // Interp omega to cos(theta) from file

      // Grading Omega
      CCTK_REAL omg_atm = (radial_distance > r_omg_atmo)
        ? (omg_surf * pow((r_omg_atmo / radial_distance), n_omg_atmo))
        : omg_surf;
      
      double alpL = calc_avg_v2c(alp, p);

      const vec<CCTK_REAL, 3> betas_avg(
          [&](int i) ARITH_INLINE { return calc_avg_v2c(gf_beta(i), p); });

      velx(p.I) = (-p.y * omg_atm + betas_avg(0)) / alpL;
      vely(p.I) = (p.x * omg_atm + betas_avg(1)) / alpL;

    } else {
      velx(p.I) = velxL;
      vely(p.I) = velyL;
    }
                             });

}

} // namespace RNSReader
