#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

#include "seedavec_utils.hxx"

namespace SeedAvecInspiral {
using namespace std;
using namespace Loop;
using namespace AsterUtils;

extern "C" void SeedAvecInspiral_Initialize_Seeding_Flags(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_SeedAvecInspiral_Initialize_Seeding_Flags;
  DECLARE_CCTK_PARAMETERS;

  *SeedNow = 0;
  *DoneSeeding = 0;
  for (int ii=0; ii<3; ii++) {
    vel_NS1[ii] = 0.0;
    vel_NS2[ii] = 0.0;
  }
  return;
}

extern "C" void SeedAvecInspiral_Set_Seeding_Flags(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_SeedAvecInspiral_Set_Seeding_Flags;
  DECLARE_CCTK_PARAMETERS;

  bool seed_trigger = false;
  if (use_time && (cctk_time > seeding_time))
    seed_trigger = true;
  if (use_separation) {
    const CCTK_REAL ns_Dx = abs(*comx1 - *comx2);
    const CCTK_REAL ns_Dy = abs(*comy1 - *comy2);
    const CCTK_REAL ns_sep = sqrt(ns_Dx * ns_Dx + ns_Dy * ns_Dy);
    if (ns_sep < seeding_separation)
      seed_trigger = true;
  }

  if ((!*DoneSeeding) && (cctk_iteration % seed_every == 0) && seed_trigger) {
    CCTK_VINFO("Seeding Initial Avec...");
    *SeedNow = 1;
  } else {
    *SeedNow = 0;
  }
  
  return;
}

extern "C" void SeedAvecInspiral_InitializeCenteredAvec_BNS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_SeedAvecInspiral_InitializeCenteredAvec_BNS;
  DECLARE_CCTK_PARAMETERS;

  // Set origin according to parameter or NS CoMs
  CCTK_REAL x01, y01, z01, x02, y02, z02;
  if (seed_every != 0) {
    x01 = *comx1;
    y01 = *comy1;
    z01 = *comz1;
    x02 = *comx2;
    y02 = *comy2;
    z02 = *comz2;
  }
  else {
    x01 = dipole_x[0]; 
    y01 = dipole_y[0]; 
    z01 = dipole_z[0]; 
    x02 = dipole_x[1]; 
    y02 = dipole_y[1]; 
    z02 = dipole_z[1]; 
  }

  if (CCTK_EQUALS(Afield_config, "internal dipole")) {

    /* computing cell centered vector potential components */
    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          // dummy initialization
          CCTK_REAL x_local = 0.0;
          CCTK_REAL y_local = 0.0;

          // For star 1 at minus side
          if (p.x < 0) {
            x_local = p.x - x01;
            y_local = p.y - y01;
          }
          // For star 2 at plus side
          else {
            x_local = p.x - x02;
            y_local = p.y - y02;
          }

          CCTK_REAL Pcut = press_max * press_cut;
          CCTK_REAL Pdiff = std::max(press(p.I) - Pcut, 0.0);
          CCTK_REAL Aphi_local = Ab * pow(Pdiff, Avec_kappa);
          Avec_x_cent(p.I) = -y_local * Aphi_local;
          Avec_y_cent(p.I) = x_local * Aphi_local;
          Avec_z_cent(p.I) = 0.0;
        });

  } else if (CCTK_EQUALS(Afield_config, "external dipole")) {

    /* computing cell centered vector potential components */
    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          // For star 1 at minus side
          CCTK_REAL x_local_s1 = p.x - x01;
          CCTK_REAL y_local_s1 = p.y - y01;
          CCTK_REAL z_local_s1 = p.z - z01;
          CCTK_REAL cylrad2_s1 =
              x_local_s1 * x_local_s1 + y_local_s1 * y_local_s1;
          CCTK_REAL rsph_s1 =
              sqrt(x_local_s1 * x_local_s1 + y_local_s1 * y_local_s1 +
                   z_local_s1 * z_local_s1);

          CCTK_REAL Aphi_local_s1 =
              B0 * (pow(r0, 3.0) / (pow(r0, 3.0) + pow(rsph_s1, 3.0))) /
              sqrt(cylrad2_s1 + 1.0e-16);

          // For star 2 at minus side
          CCTK_REAL x_local_s2 = p.x - x02;
          CCTK_REAL y_local_s2 = p.y - y02;
          CCTK_REAL z_local_s2 = p.z - z02;
          CCTK_REAL cylrad2_s2 =
              x_local_s2 * x_local_s2 + y_local_s2 * y_local_s2;
          CCTK_REAL rsph_s2 =
              sqrt(x_local_s2 * x_local_s2 + y_local_s2 * y_local_s2 +
                   z_local_s2 * z_local_s2);

          CCTK_REAL Aphi_local_s2 =
              B0 * (pow(r0, 3.0) / (pow(r0, 3.0) + pow(rsph_s2, 3.0))) /
              sqrt(cylrad2_s2 + 1.0e-16);

          Avec_x_cent(p.I) =
              -(y_local_s1 * Aphi_local_s1 + y_local_s2 * Aphi_local_s2);
          Avec_y_cent(p.I) =
              x_local_s1 * Aphi_local_s1 + x_local_s2 * Aphi_local_s2;
          Avec_z_cent(p.I) = 0.0;
        });

  } else if (CCTK_EQUALS(Afield_config, "external dipole UIUC")) {

    /* computing cell centered vector potential components */
    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

          const CCTK_REAL pi = 2 * acos(0.0);

          // For star 1 at minus side
          CCTK_REAL x_local_s1 = p.x - x01;
          CCTK_REAL y_local_s1 = p.y - y01;
          CCTK_REAL z_local_s1 = p.z - z01;
          CCTK_REAL cylrad2_s1 =
              x_local_s1 * x_local_s1 + y_local_s1 * y_local_s1;
          CCTK_REAL sphrad2_s1 =
              x_local_s1 * x_local_s1 + y_local_s1 * y_local_s1 +
                   z_local_s1 * z_local_s1;
          CCTK_REAL r02 = r0 * r0;

          // See e.g. Ruiz+ 2020, Eq. 1. Here, B0 is I0 from Eq. 1, and we have preemptively canceled out the factor
          // of cylrad2.
          CCTK_REAL Aphi_local_s1 =
              pi * B0 * r02 / pow(r02 + sphrad2_s1, 1.5) 
              * (1.0 + (15.0 * r02 * (r02 + cylrad2_s1) / (8.0 * pow(r02 + sphrad2_s1, 2.0))));

          // For star 2 at minus side
          CCTK_REAL x_local_s2 = p.x - x02;
          CCTK_REAL y_local_s2 = p.y - y02;
          CCTK_REAL z_local_s2 = p.z - z02;
          CCTK_REAL cylrad2_s2 =
              x_local_s2 * x_local_s2 + y_local_s2 * y_local_s2;
          CCTK_REAL sphrad2_s2 =
              x_local_s2 * x_local_s2 + y_local_s2 * y_local_s2 +
                   z_local_s2 * z_local_s2;

          CCTK_REAL Aphi_local_s2 =
              pi * B0 * r02 / pow(r02 + sphrad2_s2, 1.5) 
              * (1.0 + (15.0 * r02 * (r02 + cylrad2_s2) / (8.0 * pow(r02 + sphrad2_s2, 2.0))));

          Avec_x_cent(p.I) =
              -(y_local_s1 * Aphi_local_s1 + y_local_s2 * Aphi_local_s2);
          Avec_y_cent(p.I) =
              x_local_s1 * Aphi_local_s1 + x_local_s2 * Aphi_local_s2;
          Avec_z_cent(p.I) = 0.0;
        });
  } else {
    CCTK_ERROR("Vector potential configuration not defined");
  }

  *DoneSeeding = 1;
}

extern "C" void SeedAvecInspiral_InitializeStagAvec_BNS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_SeedAvecInspiral_InitializeStagAvec_BNS;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int<1, 0, 0>(grid.nghostzones,
                         [=] CCTK_HOST(const Loop::PointDesc &p)
                             CCTK_ATTRIBUTE_ALWAYS_INLINE {
                               Avec_x(p.I) = calc_avg_c2e<0>(Avec_x_cent, p);
                             });

  grid.loop_int<0, 1, 0>(grid.nghostzones,
                         [=] CCTK_HOST(const Loop::PointDesc &p)
                             CCTK_ATTRIBUTE_ALWAYS_INLINE {
                               Avec_y(p.I) = calc_avg_c2e<1>(Avec_y_cent, p);
                             });

  grid.loop_int<0, 0, 1>(grid.nghostzones,
                         [=] CCTK_HOST(const Loop::PointDesc &p)
                             CCTK_ATTRIBUTE_ALWAYS_INLINE {
                               Avec_z(p.I) = calc_avg_c2e<2>(Avec_z_cent, p);
                             });
}

extern "C" void SeedAvecInspiral_InitializeAvectoZero(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_SeedAvecInspiral_InitializeAvectoZero;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all<1, 0, 0>(grid.nghostzones,
                         [=] CCTK_HOST(const Loop::PointDesc &p)
                             CCTK_ATTRIBUTE_ALWAYS_INLINE {
                               Avec_x(p.I) = 0.0;
                             });

  grid.loop_all<0, 1, 0>(grid.nghostzones,
                         [=] CCTK_HOST(const Loop::PointDesc &p)
                             CCTK_ATTRIBUTE_ALWAYS_INLINE {
                               Avec_y(p.I) = 0.0;
                             });

  grid.loop_all<0, 0, 1>(grid.nghostzones,
                         [=] CCTK_HOST(const Loop::PointDesc &p)
                             CCTK_ATTRIBUTE_ALWAYS_INLINE {
                               Avec_z(p.I) = 0.0;
                             });
}

} // namespace SeedAvecInspiral
