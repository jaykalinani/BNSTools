#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "util_Table.h"
#include <cassert>
#include <cmath>
#include <cstring>

#include "integrands.cxx"

using namespace Loop;
extern "C" void VI_GRMHDX_ComputeIntegrand(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_VI_GRMHDX_ComputeIntegrand;
  DECLARE_CCTK_PARAMETERS;

  // rhs: m_neutron_MeV * MeV_to_erg * PRESSGF * LENGTHGF^3
  // TODO: Put in param.ccl
  constexpr double my_baryon_mass = 939.565379 * 1.60217733e-6 * 1.80123683248503e-39 * (6.77269222552442e-06*6.77269222552442e-06*6.77269222552442e-06);

  int which_integral = NumIntegrals - *IntegralCounter + 1;
  if (CCTK_MyProc == 0) {
    CCTK_VINFO("Computing Integrand %d", which_integral);
  }

  /* Note: Must extend this if/else statement if adding a new integrand! */
  if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                  "centerofmass")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CoM_integrand(VolIntegrand1, VolIntegrand2, VolIntegrand3,
                        VolIntegrand4, p, w_lorentz, rho, gxx, gxy, gxz,
                        gyy, gyz, gzz);
        });
  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "coordvolume")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CoordVol_integrand(VolIntegrand1, p, rho, dens_atmo);
        });
  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "restmass")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          M0_integrand(VolIntegrand1, p, w_lorentz, rho, gxx, gxy, gxz,
                       gyy, gyz, gzz);
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
                         "density_weighted_norm_B_field")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          mean_density_weighted_B(VolIntegrand1, VolIntegrand2, p, w_lorentz, rho, Bvecx, Bvecy, Bvecz, gxx, gxy, gxz,
                       gyy, gyz, gzz);
        });
  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "kinetic_energy")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

          // TODO: Get BNS COM from --BNSTrackerGen-- somewhere
          // Make sure we have updated values 
          // before this integral is called!

          double cms_x = 0.0;
          double cms_y = 0.0;

          kinetic(VolIntegrand1, VolIntegrand2, VolIntegrand3, p, velx, vely, velz, w_lorentz, rho, eps, press, gxx, gxy, gxz,
                       gyy, gyz, gzz, cms_x, cms_y);
        });
  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "kinetic_energy_palenzuela")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

          // TODO: Get BNS COM from --BNSTrackerGen-- somewhere
          // Make sure we have updated values 
          // before this integral is called!

          double cms_x = 0.0;
          double cms_y = 0.0;

          kinetic_palenzuela(VolIntegrand1, VolIntegrand2, VolIntegrand3, p, velx, vely, velz, w_lorentz, rho, eps, press,
                       alp, betax, betay, betaz, gxx, gxy, gxz,
                       gyy, gyz, gzz, cms_x, cms_y);
        });
  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "kinetic_energy_shibata")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

          // TODO: Get BNS COM from --BNSTrackerGen-- somewhere
          // Make sure we have updated values 
          // before this integral is called!

          double cms_x = 0.0;
          double cms_y = 0.0;

          kinetic_shibata(VolIntegrand1, VolIntegrand2, VolIntegrand3, p, velx, vely, velz, w_lorentz, rho, eps, press,
                       alp, betax, betay, betaz, gxx, gxy, gxz,
                       gyy, gyz, gzz, cms_x, cms_y);
        });
  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "kinetic_energy_total")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

          kinetic_tot(VolIntegrand1, p, velx, vely, velz, w_lorentz, rho, eps, press,
                       gxx, gxy, gxz, gyy, gyz, gzz);
        });
  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "thermal_energy")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

          thermal(VolIntegrand1, VolIntegrand2, VolIntegrand3, p, w_lorentz, rho, eps, entropy,
                       gxx, gxy, gxz, gyy, gyz, gzz, my_baryon_mass);
        });
  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "magnetic_energy_total")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

          // TODO: Get BNS COM from --BNSTrackerGen-- somewhere
          // Make sure we have updated values 
          // before this integral is called!

          double cms_x = 0.0;
          double cms_y = 0.0;

          magnetic_tot(VolIntegrand1, VolIntegrand2, p, velx, vely, velz, Bvecx, Bvecy, Bvecz,
                       gxx, gxy, gxz, gyy, gyz, gzz, cms_x, cms_y);
        });
  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "em_energy_ab")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

          // TODO: Get BNS COM from --BNSTrackerGen-- somewhere
          // Make sure we have updated values 
          // before this integral is called!

          double cms_x = 0.0;
          double cms_y = 0.0;

          magnetic_tot_12(VolIntegrand1, VolIntegrand2, p, rho, dens_a, dens_b, velx, vely, velz, Bvecx, Bvecy, Bvecz,
                       gxx, gxy, gxz, gyy, gyz, gzz, cms_x, cms_y);
        });
  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "em_energy_cd")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

          // TODO: Get BNS COM from --BNSTrackerGen-- somewhere
          // Make sure we have updated values 
          // before this integral is called!

          double cms_x = 0.0;
          double cms_y = 0.0;

          magnetic_tot_12(VolIntegrand1, VolIntegrand2, p, rho, dens_c, dens_d, velx, vely, velz, Bvecx, Bvecy, Bvecz,
                       gxx, gxy, gxz, gyy, gyz, gzz, cms_x, cms_y);
        });
  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "volume_average_norm_B_ab")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

          // TODO: Get BNS COM from --BNSTrackerGen-- somewhere
          // Make sure we have updated values 
          // before this integral is called!

          double cms_x = 0.0;
          double cms_y = 0.0;

          volume_norm_B_12(VolIntegrand1, VolIntegrand2, VolIntegrand3, VolIntegrand4, p, rho, dens_a, dens_b, Bvecx, Bvecy, Bvecz,
                       gxx, gxy, gxz, gyy, gyz, gzz, cms_x, cms_y);
        });
  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "volume_average_norm_B_cd")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

          // TODO: Get BNS COM from --BNSTrackerGen-- somewhere
          // Make sure we have updated values 
          // before this integral is called!

          double cms_x = 0.0;
          double cms_y = 0.0;

          volume_norm_B_12(VolIntegrand1, VolIntegrand2, VolIntegrand3, VolIntegrand4, p, rho, dens_c, dens_d, Bvecx, Bvecy, Bvecz,
                       gxx, gxy, gxz, gyy, gyz, gzz, cms_x, cms_y);
        });
  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "magnetic_energy_comov")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

          magnetic_co(VolIntegrand1, p, smallb2, w_lorentz,
                       gxx, gxy, gxz, gyy, gyz, gzz);
        });
  } else if (CCTK_EQUALS(Integration_quantity_keyword[which_integral],
                         "magnetic_energy_comov")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

          magnetic_co(VolIntegrand1, p, smallb2, w_lorentz,
                       gxx, gxy, gxz, gyy, gyz, gzz);
        });
  } else {
    /* Print a warning if no integrand is computed because
     * Integration_quantity_keyword unrecognized. */
    printf("VolumeIntegrals: WARNING: Integrand not computed. Did not "
           "understand Integration_quantity_keyword[%d] = %s\n",
           which_integral, Integration_quantity_keyword[which_integral]);
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
    int which_centre = volintegral_sphere__tracks__amr_centre[which_integral];

    volintegral_inside_sphere__center_x[which_integral] =
        position_x[which_centre];
    volintegral_inside_sphere__center_y[which_integral] =
        position_y[which_centre];
    volintegral_inside_sphere__center_z[which_integral] =
        position_z[which_centre];

    volintegral_outside_sphere__center_x[which_integral] =
        position_x[which_centre];
    volintegral_outside_sphere__center_y[which_integral] =
        position_y[which_centre];
    volintegral_outside_sphere__center_z[which_integral] =
        position_z[which_centre];
  }

  /* ZERO OUT INTEGRATION REGIONS */

  /* The below code supports zeroing out of arbitrary spherical shells.
     In the case of integration INSIDE a full sphere, this code also supports
     moving spheres if AMR centre tracking is enabled. This can be used to
     track compact objects.

     Here's one way:
     Generally the lapse at the center of these objects is minimized.
     If the object has a length scale of R and is known to be centered at x,y,z,
     compute X^i = Integral [ (1-lapse) x^i dV] / Integral [(1-lapse) dV] over
     the spherical volume centered at x,y,z with radius R. This should yield
     X,Y,Z = x,y,z to good approximation. If you enable AMR centre tracking, it
     should work just as well
     as any other method for tracking compact objects, if not better. */

  /* Set integrands to zero outside a sphere centered at x,y,z.
     I.e., this results in the integral being restricted INSIDE sphere */
  if (volintegral_inside_sphere__radius[which_integral] > 0.0) {
    double radius = volintegral_inside_sphere__radius[which_integral];
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
          double x_minus_xprime = p.x - xprime;
          double y_minus_xprime = p.y - yprime;
          double z_minus_xprime = p.z - zprime;
          if (sqrt(x_minus_xprime * x_minus_xprime +
                   y_minus_xprime * y_minus_xprime +
                   z_minus_xprime * z_minus_xprime) > radius)
            VolIntegrand1(p.I) = VolIntegrand2(p.I) = VolIntegrand3(p.I) =
                VolIntegrand4(p.I) = 0.0;
        });
  }

  /* Set integrands to zero inside a sphere centered at x,y,z.
     I.e., this results in the integral being restricted OUTSIDE sphere.
     Combine this with above to get spherical shell.

     Note that volume integrals outside a sphere are fixed at the
     original sphere center position for all time, unlike volume
     integrals inside a sphere. We do this because the latter are
     generally used for tracking moving compact objects or other things. */
  if (volintegral_outside_sphere__radius[which_integral] > 0.0) {
    double radius = volintegral_outside_sphere__radius[which_integral];
    double xprime = volintegral_outside_sphere__center_x[which_integral];
    double yprime = volintegral_outside_sphere__center_y[which_integral];
    double zprime = volintegral_outside_sphere__center_z[which_integral];
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          double x_minus_xprime = p.x - xprime;
          double y_minus_xprime = p.y - yprime;
          double z_minus_xprime = p.z - zprime;
          if (sqrt(x_minus_xprime * x_minus_xprime +
                   y_minus_xprime * y_minus_xprime +
                   z_minus_xprime * z_minus_xprime) <= radius)
            VolIntegrand1(p.I) = VolIntegrand2(p.I) = VolIntegrand3(p.I) =
                VolIntegrand4(p.I) = 0.0;
        });
  }
}
