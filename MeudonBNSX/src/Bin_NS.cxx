
#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>
#include <array>
#include <vec.hxx>
#include <vector>
#include <ios>
#include <iostream>
#include <stdlib.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <bin_ns.h>
#include <unites.h>

#include <loop.hxx>
#include <loop_device.hxx>
#include <meudonbnsx_utils.hxx>

namespace MeudonBNSX {
using namespace std;
using namespace Arith;
using namespace Loop;

// define namespace here for old versions of Lorene that don't do so
namespace Lorene {}
using namespace Lorene;

extern "C"
void MeudonBNSX_initialise(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTSX_MeudonBNSX_initialise;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO("Setting up LORENE Bin_NS initial data");

  // Meudon data are distributed in SI units (MKSA).  Here are some
  // conversion factors.
  // Be aware: these are the constants Lorene uses. They do differ from other
  // conventions, but they gave the best results in some tests.

  CCTK_REAL const c_light  = Unites::c_si;      // speed of light [m/s]
  CCTK_REAL const nuc_dens = Unites::rhonuc_si; // Nuclear density as used in Lorene units [kg/m^3]
  CCTK_REAL const G_grav   = Unites::g_si;      // gravitational constant [m^3/kg/s^2]
  CCTK_REAL const M_sun    = Unites::msol_si;   // solar mass [kg]

  // Cactus units in terms of SI units:
  // (These are derived from M = M_sun, c = G = 1, and using 1/M_sun
  // for the magnetic field)
  CCTK_REAL const cactusM = M_sun;
  CCTK_REAL const cactusL = cactusM * G_grav / pow(c_light,2);
  CCTK_REAL const cactusT = cactusL / c_light;

  // Other quantities in terms of Cactus units
  CCTK_REAL const coord_unit = cactusL / 1.0e+3;         // from km (~1.477)
  CCTK_REAL const rho_unit   = cactusM / pow(cactusL,3); // from kg/m^3

  CCTK_INFO ("Setting up coordinates");

  int const npoints = (cctk_lsh[0]-1) * (cctk_lsh[1]-1) * (cctk_lsh[2]-1);
  vector<double> xx(npoints), yy(npoints), zz(npoints);
  
  //TODO: currently works only with polytropic EOS 
  CCTK_INFO("MeudonBNSX will use the polytropic equation of state.");
    
  grid.loop_all<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const array<CCTK_INT, dim> indextype = {1, 1, 1};
    const GF3D2layout layout(cctkGH, indextype);
    CCTK_INT idx = layout.linear(p.I);
    xx[idx] = p.X[0] * coord_unit;
    yy[idx] = p.X[1] * coord_unit;
    zz[idx] = p.X[2] * coord_unit;
  });
  
  // --------------------------------------------------------------
  //   CHECKING FILE NAME EXISTENCE
  // --------------------------------------------------------------
  FILE *file;
  if ((file = fopen(filename, "r")) != NULL) 
     fclose(file);
  else {
     CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                 "File \"%s\" does not exist. ABORTING", filename);
  }
  // Handle potentially different EOS table directory. LORENE recieves that via
  // environment variable
  if (strlen(eos_table_filepath) > 0) {
    if (setenv("LORENE_TABULATED_EOS_PATH", eos_table_filepath, 1)) {
      CCTK_ERROR("Unable to set environment variable LORENE_TABULATED_EOS_PATH");

    }
  }


  CCTK_VInfo (CCTK_THORNSTRING, "Reading from file \"%s\"", filename);

  //try {
  Bin_NS bin_ns (npoints, &xx[0], &yy[0], &zz[0], filename);

  CCTK_VInfo (CCTK_THORNSTRING, "omega [rad/s]:       %g", bin_ns.omega);
  CCTK_VInfo (CCTK_THORNSTRING, "dist [km]:           %g", bin_ns.dist);
  CCTK_VInfo (CCTK_THORNSTRING, "dist_mass [km]:      %g", bin_ns.dist_mass);
  CCTK_VInfo (CCTK_THORNSTRING, "mass1_b [M_sun]:     %g", bin_ns.mass1_b);
  CCTK_VInfo (CCTK_THORNSTRING, "mass2_b [M_sun]:     %g", bin_ns.mass2_b);
  CCTK_VInfo (CCTK_THORNSTRING, "mass_ADM [M_sun]:    %g", bin_ns.mass_adm);
  CCTK_VInfo (CCTK_THORNSTRING, "L_tot [G M_sun^2/c]: %g", bin_ns.angu_mom);
  CCTK_VInfo (CCTK_THORNSTRING, "rad1_x_comp [km]:    %g", bin_ns.rad1_x_comp);
  CCTK_VInfo (CCTK_THORNSTRING, "rad1_y [km]:         %g", bin_ns.rad1_y);
  CCTK_VInfo (CCTK_THORNSTRING, "rad1_z [km]:         %g", bin_ns.rad1_z);
  CCTK_VInfo (CCTK_THORNSTRING, "rad1_x_opp [km]:     %g", bin_ns.rad1_x_opp);
  CCTK_VInfo (CCTK_THORNSTRING, "rad2_x_comp [km]:    %g", bin_ns.rad2_x_comp);
  CCTK_VInfo (CCTK_THORNSTRING, "rad2_y [km]:         %g", bin_ns.rad2_y);
  CCTK_VInfo (CCTK_THORNSTRING, "rad2_z [km]:         %g", bin_ns.rad2_z);
  CCTK_VInfo (CCTK_THORNSTRING, "rad2_x_opp [km]:     %g", bin_ns.rad2_x_opp);
  // LORENE's EOS is in terms on number density n = rho/m_nucleon:
  // P = K n^Gamma
  // to convert to SI units:
  // K_SI(n) = K_LORENE rho_nuc c^2 / n_nuc^gamma
  // Converting this to be in terms of the mass density rho = n m_nucleon gets
  // changes n_nuc to rho_nuc:
  // K_SI(rho) = K_LORENE c^2 / rho_nuc^(gamma-1)
  // In SI units P has units of M / (L T^2) and rho has units of M/L^3 thus
  // K_SI has units of (L^3/M)^Gamma M/(L T^2).
  // In Cactus units P and rho have the same units thus K_Cactus is unitless.
  // Conversion between K_SI and K_Cactus thus amounts to dividing out the
  // units of the SI quantity.
  double K = bin_ns.kappa_poly1 * pow((pow(c_light, 6.0) /
             ( pow(G_grav, 3.0) * M_sun * M_sun *
               nuc_dens )),bin_ns.gamma_poly1-1.);
  CCTK_VInfo (CCTK_THORNSTRING, "K [ET unit]:         %.15g", K);

  assert (bin_ns.np == npoints);

  CCTK_INFO ("Filling in Cactus grid points");

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const array<CCTK_INT, dim> indextype = {1, 1, 1};
    const GF3D2layout layout(cctkGH, indextype);
    CCTK_INT idx = layout.linear(p.I);

    if (CCTK_EQUALS(initial_lapse, "MeudonBNSX")) { 
      alp_cc(p.I) = bin_ns.nnn[idx];
    }

    if (CCTK_EQUALS(initial_shift, "MeudonBNSX")) { 
      betax_cc(p.I) = -bin_ns.beta_x[idx];
      betay_cc(p.I) = -bin_ns.beta_y[idx];
      betaz_cc(p.I) = -bin_ns.beta_z[idx];
    }

    if (CCTK_EQUALS(initial_data, "MeudonBNSX")) {
      gxx_cc(p.I) = bin_ns.g_xx[idx];
      gxy_cc(p.I) = bin_ns.g_xy[idx];
      gxz_cc(p.I) = bin_ns.g_xz[idx];
      gyy_cc(p.I) = bin_ns.g_yy[idx];
      gyz_cc(p.I) = bin_ns.g_yz[idx];
      gzz_cc(p.I) = bin_ns.g_zz[idx];

      kxx_cc(p.I) = bin_ns.k_xx[idx] * coord_unit;
      kxy_cc(p.I) = bin_ns.k_xy[idx] * coord_unit;
      kxz_cc(p.I) = bin_ns.k_xz[idx] * coord_unit;
      kyy_cc(p.I) = bin_ns.k_yy[idx] * coord_unit;
      kyz_cc(p.I) = bin_ns.k_yz[idx] * coord_unit;
      kzz_cc(p.I) = bin_ns.k_zz[idx] * coord_unit;
    }

    if (CCTK_EQUALS(initial_data, "MeudonBNSX")) {
      rho(p.I) = bin_ns.nbar[idx] / rho_unit;
      if (!recalculate_eps)
        eps(p.I) = bin_ns.ener_spec[idx];

      //TODO: currently works only with polytropic EOS 
      press(p.I) = K * pow(rho(p.I), bin_ns.gamma_poly1);

      velx(p.I) = bin_ns.u_euler_x[idx];
      vely(p.I) = bin_ns.u_euler_y[idx];
      velz(p.I) = bin_ns.u_euler_z[idx];

      // Especially the velocity is set to strange values outside of the
      // matter region, so take care of this in the following way
      if (rho(p.I) < 1.e-30) {
        rho(p.I) = 1.e-30;
        velx(p.I) = 0.0;
        vely(p.I) = 0.0;
        velz(p.I) = 0.0;
        eps(p.I) = K * pow(rho(p.I), bin_ns.gamma_poly1-1.) / (bin_ns.gamma_poly1-1.);
        press(p.I) = K * pow(rho(p.I), bin_ns.gamma_poly1);
      }
    }

//  });

//  {

    if (CCTK_EQUALS(initial_lapse, "MeudonBNSX")) { 
      if (CCTK_EQUALS (initial_dtlapse, "MeudonBNSX")) {
        CCTK_ERROR("Code for computing time derivatives of lapse is not yet implemented");
      } else if (CCTK_EQUALS (initial_dtlapse, "none") or CCTK_EQUALS(initial_dtlapse,"zero")) {
        dtalp_cc(p.I) = 0.0;
      } else {
        CCTK_WARN (CCTK_WARN_ABORT, "internal error");
      }
    }

    if (CCTK_EQUALS(initial_shift, "MeudonBNSX")) { 
      if (CCTK_EQUALS (initial_dtshift, "MeudonBNSX")) {
        CCTK_ERROR("Code for calculating time derivatives of shift is not yet implemented");
      } else if (CCTK_EQUALS (initial_dtshift, "none") or CCTK_EQUALS(initial_dtshift,"zero")) {
        dtbetax_cc(p.I) = 0.0;
        dtbetay_cc(p.I) = 0.0;
        dtbetaz_cc(p.I) = 0.0;
      } else {
        CCTK_WARN (CCTK_WARN_ABORT, "internal error");
      }
    }
  });

  CCTK_INFO ("Done.");
//  } catch (ios::failure e) {
//    CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
//                "Could not read initial data from file '%s': %s", filename, e.what());
//  }
}

extern "C" void MeudonBNSX_Interpolation_C2V(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MeudonBNSX_Interpolation_C2V;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int<0, 0, 0>(grid.nghostzones,
                         [=] CCTK_HOST(const Loop::PointDesc &p)
                             CCTK_ATTRIBUTE_ALWAYS_INLINE {
                               gxx(p.I) = calc_avg_c2v(gxx_cc, p);
                               gxy(p.I) = calc_avg_c2v(gxy_cc, p);
                               gxz(p.I) = calc_avg_c2v(gxz_cc, p);
                               gyy(p.I) = calc_avg_c2v(gyy_cc, p);
                               gyz(p.I) = calc_avg_c2v(gyz_cc, p);
                               gzz(p.I) = calc_avg_c2v(gzz_cc, p);

                               kxx(p.I) = calc_avg_c2v(kxx_cc, p);
                               kxy(p.I) = calc_avg_c2v(kxy_cc, p);
                               kxz(p.I) = calc_avg_c2v(kxz_cc, p);
                               kyy(p.I) = calc_avg_c2v(kyy_cc, p);
                               kyz(p.I) = calc_avg_c2v(kyz_cc, p);
                               kzz(p.I) = calc_avg_c2v(kzz_cc, p);

                               alp(p.I) = calc_avg_c2v(alp_cc, p);
                               betax(p.I) = calc_avg_c2v(betax_cc, p);
                               betay(p.I) = calc_avg_c2v(betay_cc, p);
                               betaz(p.I) = calc_avg_c2v(betaz_cc, p);

                               dtalp(p.I) = calc_avg_c2v(dtalp_cc, p);
                               dtbetax(p.I) = calc_avg_c2v(dtbetax_cc, p);
                               dtbetay(p.I) = calc_avg_c2v(dtbetay_cc, p);
                               dtbetaz(p.I) = calc_avg_c2v(dtbetaz_cc, p);
                             });

  CCTK_INFO("Done interpolation for TOV ADM variables.");

} 

} //namespace
