#ifndef INTEGRANDS_CXX
#define INTEGRANDS_CXX

#include <loop_device.hxx>
#include <VI_utils.hxx>

/*
  NRPy+ Python code that generates below C code:
  (Note: used commit 6a86d39aa3e4c863caa51c2fc7a04f3d06727167 in
         https://github.com/zachetienne/nrpytutorial)

import sympy as sp
import grid as gri
import indexedexp as ixp
from outputC import *

DIM=3

# Need to read from memory & write to memory manually due to speed limiter if()
statement. with open("compute_rhostar.h","w") as file: file.write("""
// This function computes
//  rho_* = alpha*u^0*sqrt(gamma)*rho_0
//        = w_lorentz*sqrt(gamma)*rho_0 ,
// where gamma = determinant of 3-metric,
// rho_0 is the rest-mass density, and
// w_lorentz is the Lorentz factor\n\n""")

    for i in range(DIM):
        for j in range(i,DIM):
            file.write("const CCTK_REAL gammaDD"+str(i)+str(j)+" =
g"+chr(ord('x')+i)+chr(ord('x')+j)+"[index];\n") file.write(""" const CCTK_REAL
w_lorentz = w_lorentz[index]; const CCTK_REAL rho0      = rhozero[  index];

CCTK_REAL w_lorentz_limited = w_lorentz;
if(w_lorentz > CoM_integrand_GAMMA_SPEED_LIMIT) w_lorentz_limited =
CoM_integrand_GAMMA_SPEED_LIMIT;
""")

    rho0,w_lorentz_limited = sp.symbols("rho0 w_lorentz_limited")
    gammaDD                = ixp.declarerank2("gammaDD", "sym01",DIM=3)

    dummy, detgamma = ixp.symm_matrix_inverter3x3(gammaDD)

    rhostar = w_lorentz_limited*sp.sqrt(detgamma)*rho0

    rho_star_str_ugly = outputC(rhostar, "const CCTK_REAL rhostar",
filename="returnstring",params="includebraces=False") # Beautify
rho_star_str_ugly rho_star_str = rho_star_str_ugly.replace( "pow(gammaDD02,
2)","(gammaDD02*gammaDD02)").replace( "pow(gammaDD12,
2)","(gammaDD12*gammaDD12)").replace( "pow(gammaDD01,
2)","(gammaDD01*gammaDD01)")

    file.write(rho_star_str)
*/

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
compute_rho_star(
    const Loop::PointDesc &p, const Loop::GF3D2<const double> w_lorentzGF,
    const Loop::GF3D2<const double> rho0GF, const Loop::GF3D2<const double> gxx,
    const Loop::GF3D2<const double> gxy, const Loop::GF3D2<const double> gxz,
    const Loop::GF3D2<const double> gyy, const Loop::GF3D2<const double> gyz,
    const Loop::GF3D2<const double> gzz, const CCTK_REAL gamma_lim) {

  // This function computes
  //  rho_* = alpha*u^0*sqrt(gamma)*rho_0
  //        = w_lorentz*sqrt(gamma)*rho_0 ,
  // where gamma = determinant of 3-metric,
  // rho_0 is the rest-mass density, and
  // w_lorentz is the Lorentz factor

  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};
  /* Get covariant metric */
  const smat<CCTK_REAL, 3> glo(
      [&](int i, int j) ARITH_INLINE { return calc_avg_v2c(gf_g(i, j), p); });

  const CCTK_REAL gammaDD00 = glo(0, 0);
  const CCTK_REAL gammaDD01 = glo(0, 1);
  const CCTK_REAL gammaDD02 = glo(0, 2);
  const CCTK_REAL gammaDD11 = glo(1, 1);
  const CCTK_REAL gammaDD12 = glo(1, 2);
  const CCTK_REAL gammaDD22 = glo(2, 2);

  const CCTK_REAL w_lorentz = w_lorentzGF(p.I);
  const CCTK_REAL rho0 = rho0GF(p.I);

  CCTK_REAL w_lorentz_limited = w_lorentz;
  if (w_lorentz > gamma_lim)
    w_lorentz_limited = gamma_lim;
  /*
   *  Original SymPy expression:
   *  "rhostar = rho0*w_lorentz_limited*sqrt(gammaDD00*gammaDD11*gammaDD22 -
   * gammaDD00*gammaDD12**2 - gammaDD01**2*gammaDD22 +
   * 2*gammaDD01*gammaDD02*gammaDD12 - gammaDD02**2*gammaDD11)"
   */
  const CCTK_REAL rhostar = rho0 * w_lorentz_limited *
                            sqrt(gammaDD00 * gammaDD11 * gammaDD22 -
                                 gammaDD00 * (gammaDD12 * gammaDD12) -
                                 (gammaDD01 * gammaDD01) * gammaDD22 +
                                 2 * gammaDD01 * gammaDD02 * gammaDD12 -
                                 (gammaDD02 * gammaDD02) * gammaDD11);

  return rhostar;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
compute_sqrtgamma(const Loop::PointDesc &p, const Loop::GF3D2<const double> gxx,
                  const Loop::GF3D2<const double> gxy,
                  const Loop::GF3D2<const double> gxz,
                  const Loop::GF3D2<const double> gyy,
                  const Loop::GF3D2<const double> gyz,
                  const Loop::GF3D2<const double> gzz) {

  const CCTK_REAL gammaDD00 = gxx(p.I);
  const CCTK_REAL gammaDD01 = gxy(p.I);
  const CCTK_REAL gammaDD02 = gxz(p.I);
  const CCTK_REAL gammaDD11 = gyy(p.I);
  const CCTK_REAL gammaDD12 = gyz(p.I);
  const CCTK_REAL gammaDD22 = gzz(p.I);

  const CCTK_REAL sqrtgamma = sqrt(gammaDD00 * gammaDD11 * gammaDD22 -
                                   gammaDD00 * (gammaDD12 * gammaDD12) -
                                   (gammaDD01 * gammaDD01) * gammaDD22 +
                                   2 * gammaDD01 * gammaDD02 * gammaDD12 -
                                   (gammaDD02 * gammaDD02) * gammaDD11);

  return sqrtgamma;
}

/* Coord Volume: */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
CoordVol_integrand(const Loop::GF3D2<double> VolIntegrand1,
                   const Loop::PointDesc &p,
                   const Loop::GF3D2<const double> rho0,
                   const CCTK_REAL dens_1) {

  const CCTK_REAL my_rho = rho0(p.I);
  if (my_rho > dens_1) {
    VolIntegrand1(p.I) = 1.0;
  } else {
    VolIntegrand1(p.I) = 0.0;
  }
}

/* Center of Mass: */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void CoM_integrand(
    const Loop::GF3D2<double> VolIntegrand1,
    const Loop::GF3D2<double> VolIntegrand2,
    const Loop::GF3D2<double> VolIntegrand3,
    const Loop::GF3D2<double> VolIntegrand4, const Loop::PointDesc &p,
    const Loop::GF3D2<const double> w_lorentz,
    const Loop::GF3D2<const double> rho0, const Loop::GF3D2<const double> gxx,
    const Loop::GF3D2<const double> gxy, const Loop::GF3D2<const double> gxz,
    const Loop::GF3D2<const double> gyy, const Loop::GF3D2<const double> gyz,
    const Loop::GF3D2<const double> gzz, const CCTK_REAL gamma_lim) {
  double rho_starL = compute_rho_star(p, w_lorentz, rho0, gxx, gxy, gxz, gyy,
                                      gyz, gzz, gamma_lim);
  VolIntegrand1(p.I) = rho_starL * p.x;
  VolIntegrand2(p.I) = rho_starL * p.y;
  VolIntegrand3(p.I) = rho_starL * p.z;
  VolIntegrand4(p.I) = rho_starL;
}

/* Rest Mass: */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void M0_integrand(
    const Loop::GF3D2<double> VolIntegrand1, const Loop::PointDesc &p,
    const Loop::GF3D2<const double> w_lorentz,
    const Loop::GF3D2<const double> rho0, const Loop::GF3D2<const double> gxx,
    const Loop::GF3D2<const double> gxy, const Loop::GF3D2<const double> gxz,
    const Loop::GF3D2<const double> gyy, const Loop::GF3D2<const double> gyz,
    const Loop::GF3D2<const double> gzz, const CCTK_REAL gamma_lim) {
  double rho_starL = compute_rho_star(p, w_lorentz, rho0, gxx, gxy, gxz, gyy,
                                      gyz, gzz, gamma_lim);
  VolIntegrand1(p.I) = rho_starL;
}

/* Density weighted norm of magnetic field strength: */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
mean_density_weighted_B(
    const Loop::GF3D2<double> VolIntegrand1,
    const Loop::GF3D2<double> VolIntegrand2, const Loop::PointDesc &p,
    const Loop::GF3D2<const double> Bvecx,
    const Loop::GF3D2<const double> Bvecy,
    const Loop::GF3D2<const double> Bvecz,
    const Loop::GF3D2<const double> w_lorentz,
    const Loop::GF3D2<const double> rho0, const Loop::GF3D2<const double> gxx,
    const Loop::GF3D2<const double> gxy, const Loop::GF3D2<const double> gxz,
    const Loop::GF3D2<const double> gyy, const Loop::GF3D2<const double> gyz,
    const Loop::GF3D2<const double> gzz, const CCTK_REAL gamma_lim) {
  const CCTK_REAL gammaDD00 = gxx(p.I);
  const CCTK_REAL gammaDD01 = gxy(p.I);
  const CCTK_REAL gammaDD02 = gxz(p.I);
  const CCTK_REAL gammaDD11 = gyy(p.I);
  const CCTK_REAL gammaDD12 = gyz(p.I);
  const CCTK_REAL gammaDD22 = gzz(p.I);

  const CCTK_REAL Bx = Bvecx(p.I);
  const CCTK_REAL By = Bvecy(p.I);
  const CCTK_REAL Bz = Bvecz(p.I);

  double rho_starL = compute_rho_star(p, w_lorentz, rho0, gxx, gxy, gxz, gyy,
                                      gyz, gzz, gamma_lim);

  double norm_B_sq = gammaDD00 * Bx * Bx + 2. * gammaDD01 * Bx * By +
                     2. * gammaDD02 * Bx * Bz + gammaDD11 * By * By +
                     2. * gammaDD12 * By * Bz + gammaDD22 * Bz * Bz;
  VolIntegrand1(p.I) = rho_starL * sqrt(norm_B_sq);
  VolIntegrand2(p.I) = rho_starL;
}

/* Volume averaged norm of B in region (dens_1,dens_2]: */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void volume_norm_B_12(
    const Loop::GF3D2<double> VolIntegrand1,
    const Loop::GF3D2<double> VolIntegrand2,
    const Loop::GF3D2<double> VolIntegrand3,
    const Loop::GF3D2<double> VolIntegrand4, const Loop::PointDesc &p,
    const Loop::GF3D2<const double> rho0, const CCTK_REAL dens_1,
    const CCTK_REAL dens_2, const Loop::GF3D2<const double> Bvecx,
    const Loop::GF3D2<const double> Bvecy,
    const Loop::GF3D2<const double> Bvecz, const Loop::GF3D2<const double> gxx,
    const Loop::GF3D2<const double> gxy, const Loop::GF3D2<const double> gxz,
    const Loop::GF3D2<const double> gyy, const Loop::GF3D2<const double> gyz,
    const Loop::GF3D2<const double> gzz, const double cms_x,
    const double cms_y) {

  const CCTK_REAL my_rho = rho0(p.I);

  if (my_rho > dens_1 && my_rho <= dens_2) {

    const CCTK_REAL gammaDD00 = gxx(p.I);
    const CCTK_REAL gammaDD01 = gxy(p.I);
    const CCTK_REAL gammaDD02 = gxz(p.I);
    const CCTK_REAL gammaDD11 = gyy(p.I);
    const CCTK_REAL gammaDD12 = gyz(p.I);
    const CCTK_REAL gammaDD22 = gzz(p.I);

    const CCTK_REAL B_contra[3]{Bvecx(p.I), Bvecy(p.I), Bvecz(p.I)};

    CCTK_REAL posx = p.x - cms_x;
    CCTK_REAL posy = p.y - cms_y;
    CCTK_REAL posz = p.z;

    // Cylindrical coordinates and some
    // quantities for transformation of vectors.

    const double posr2 = posx * posx + posy * posy;
    const double posr = sqrt(posr2);

    double xy_over_r2 = posx * posy / posr2;
    double x2_over_r2 = posx * posx / posr2;
    double y2_over_r2 = posy * posy / posr2;

    if (posr <= 1e-15) {

      xy_over_r2 = 0.0;
      x2_over_r2 = 0.5;
      y2_over_r2 = 0.5;
    }

    const double sqrtgamma = compute_sqrtgamma(p, gxx, gxy, gxz, gyy, gyz, gzz);

    // Covariant B

    const CCTK_REAL B_cov[3]{gammaDD00 * B_contra[0] + gammaDD01 * B_contra[1] +
                                 gammaDD02 * B_contra[2],
                             gammaDD01 * B_contra[0] + gammaDD11 * B_contra[1] +
                                 gammaDD12 * B_contra[2],
                             gammaDD02 * B_contra[0] + gammaDD12 * B_contra[1] +
                                 gammaDD22 * B_contra[2]};

    const CCTK_REAL B_sq = B_cov[0] * B_contra[0] + B_cov[1] * B_contra[1] +
                           B_cov[2] * B_contra[2];

    // Finally, compute B_phi*B^phi

    const CCTK_REAL B_2_tor = B_cov[0] * B_contra[0] * y2_over_r2 -
                              B_cov[0] * B_contra[1] * xy_over_r2 -
                              B_cov[1] * B_contra[0] * xy_over_r2 +
                              B_cov[1] * B_contra[1] * x2_over_r2;

    // Finally, B_r*B^r + B_z*B^z

    const CCTK_REAL B_2_pol = B_cov[0] * B_contra[0] * x2_over_r2 +
                              B_cov[0] * B_contra[1] * xy_over_r2 +
                              B_cov[1] * B_contra[0] * xy_over_r2 +
                              B_cov[1] * B_contra[1] * y2_over_r2 +
                              B_cov[2] * B_contra[2];

    VolIntegrand1(p.I) = sqrtgamma * sqrt(abs(B_sq));
    VolIntegrand2(p.I) = sqrtgamma * sqrt(abs(B_2_tor));
    VolIntegrand3(p.I) = sqrtgamma * sqrt(abs(B_2_pol));
    VolIntegrand4(p.I) = sqrtgamma;

  } else {

    VolIntegrand1(p.I) = 0.0;
    VolIntegrand2(p.I) = 0.0;
    VolIntegrand3(p.I) = 0.0;
    VolIntegrand4(p.I) = 0.0;
  }
}

/* Kinetic energy as given by https://arxiv.org/pdf/1206.5911.pdf */
/* See also https://arxiv.org/pdf/2102.01346.pdf for rotational energy */
/* See also https://arxiv.org/abs/astro-ph/0605331 */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void kinetic_shibata(
    const Loop::GF3D2<double> VolIntegrand1,
    const Loop::GF3D2<double> VolIntegrand2,
    const Loop::GF3D2<double> VolIntegrand3, const Loop::PointDesc &p,
    const Loop::GF3D2<const double> velxGF,
    const Loop::GF3D2<const double> velyGF,
    const Loop::GF3D2<const double> velzGF,
    const Loop::GF3D2<const double> w_lorentzGF,
    const Loop::GF3D2<const double> rho0, const Loop::GF3D2<const double> epsGF,
    const Loop::GF3D2<const double> pressGF,
    const Loop::GF3D2<const double> alpGF,
    const Loop::GF3D2<const double> betaxGF,
    const Loop::GF3D2<const double> betayGF,
    const Loop::GF3D2<const double> betazGF,
    const Loop::GF3D2<const double> gxx, const Loop::GF3D2<const double> gxy,
    const Loop::GF3D2<const double> gxz, const Loop::GF3D2<const double> gyy,
    const Loop::GF3D2<const double> gyz, const Loop::GF3D2<const double> gzz,
    const double cms_x, const double cms_y, const double gamma_lim) {

  const CCTK_REAL gammaDD00 = gxx(p.I);
  const CCTK_REAL gammaDD01 = gxy(p.I);
  const CCTK_REAL gammaDD02 = gxz(p.I);
  const CCTK_REAL gammaDD11 = gyy(p.I);
  const CCTK_REAL gammaDD12 = gyz(p.I);
  const CCTK_REAL gammaDD22 = gzz(p.I);

  const CCTK_REAL my_rho = rho0(p.I);
  const CCTK_REAL lfac = w_lorentzGF(p.I);
  const CCTK_REAL my_press = pressGF(p.I);
  const CCTK_REAL my_eps = epsGF(p.I);

  const CCTK_REAL vx = velxGF(p.I);
  const CCTK_REAL vy = velyGF(p.I);
  const CCTK_REAL vz = velzGF(p.I);

  const CCTK_REAL my_lapse = alpGF(p.I);
  const CCTK_REAL my_shiftx = betaxGF(p.I);
  const CCTK_REAL my_shifty = betayGF(p.I);
  const CCTK_REAL my_shiftz = betazGF(p.I);

  const CCTK_REAL posx = p.x - cms_x;
  const CCTK_REAL posy = p.y - cms_y;
  const CCTK_REAL posz = p.z;

  const double sqrtgamma = compute_sqrtgamma(p, gxx, gxy, gxz, gyy, gyz, gzz);

  CCTK_REAL w_lorentz_limited = lfac;

  if (lfac > gamma_lim)
    w_lorentz_limited = gamma_lim;

  // Covariant 3 velocity

  const CCTK_REAL v_cov[3]{gammaDD00 * vx + gammaDD01 * vy + gammaDD02 * vz,
                           gammaDD01 * vx + gammaDD11 * vy + gammaDD12 * vz,
                           gammaDD02 * vx + gammaDD12 * vy + gammaDD22 * vz};

  const CCTK_REAL v_2 = v_cov[0] * vx + v_cov[1] * vy + v_cov[2] * vz;

  const CCTK_REAL vz_2 = v_cov[2] * vz;

  // Scalar product with shift

  const CCTK_REAL shift_v =
      v_cov[0] * my_shiftx + v_cov[1] * my_shifty + v_cov[2] * my_shiftz;

  const CCTK_REAL shiftz_vz = v_cov[2] * my_shiftz;

  // This is just the total kinetic energy

  VolIntegrand1(p.I) = 0.5 * sqrtgamma * (my_rho * (1.0 + my_eps) + my_press) *
                       w_lorentz_limited * w_lorentz_limited *
                       (my_lapse * v_2 - shift_v);

  // Cylindrical coordinates and some
  // quantities for transformation of vectors.

  const double posr2 = posx * posx + posy * posy;
  const double posr = sqrt(posr2);

  double xy_over_r2 = posx * posy / posr2;
  double x2_over_r2 = posx * posx / posr2;
  double y2_over_r2 = posy * posy / posr2;

  if (posr <= 1e-15) {

    xy_over_r2 = 0.0;
    x2_over_r2 = 0.5;
    y2_over_r2 = 0.5;
  }

  // Finally, compute v_phi*v^phi

  CCTK_REAL v_phi2{0.0};

  v_phi2 = v_cov[0] * vx * y2_over_r2 - v_cov[0] * vy * xy_over_r2 -
           v_cov[1] * vx * xy_over_r2 + v_cov[1] * vy * x2_over_r2;

  // Finally, compute v_phi*shift^phi

  CCTK_REAL v_phi_shift_phi{0.0};

  v_phi_shift_phi =
      v_cov[0] * my_shiftx * y2_over_r2 - v_cov[0] * my_shifty * xy_over_r2 -
      v_cov[1] * my_shiftx * xy_over_r2 + v_cov[1] * my_shifty * x2_over_r2;

  // This is the rotational kinetic energy wrt to the z-axis

  VolIntegrand2(p.I) = 0.5 * sqrtgamma * (my_rho * (1.0 + my_eps) + my_press) *
                       w_lorentz_limited * w_lorentz_limited *
                       (my_lapse * v_phi2 - v_phi_shift_phi);

  // Kinetic energy in z-direction

  VolIntegrand3(p.I) = 0.5 * sqrtgamma * (my_rho * (1.0 + my_eps) + my_press) *
                       w_lorentz_limited * w_lorentz_limited *
                       (my_lapse * vz_2 - shiftz_vz);
}

/* Kinetic energy as given by https://arxiv.org/pdf/2112.08413.pdf */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
kinetic_palenzuela(
    const Loop::GF3D2<double> VolIntegrand1,
    const Loop::GF3D2<double> VolIntegrand2,
    const Loop::GF3D2<double> VolIntegrand3, const Loop::PointDesc &p,
    const Loop::GF3D2<const double> velxGF,
    const Loop::GF3D2<const double> velyGF,
    const Loop::GF3D2<const double> velzGF,
    const Loop::GF3D2<const double> w_lorentzGF,
    const Loop::GF3D2<const double> rho0, const Loop::GF3D2<const double> epsGF,
    const Loop::GF3D2<const double> pressGF,
    const Loop::GF3D2<const double> alpGF,
    const Loop::GF3D2<const double> betaxGF,
    const Loop::GF3D2<const double> betayGF,
    const Loop::GF3D2<const double> betazGF,
    const Loop::GF3D2<const double> gxx, const Loop::GF3D2<const double> gxy,
    const Loop::GF3D2<const double> gxz, const Loop::GF3D2<const double> gyy,
    const Loop::GF3D2<const double> gyz, const Loop::GF3D2<const double> gzz,
    const double cms_x, const double cms_y, const double gamma_lim) {

  const CCTK_REAL gammaDD00 = gxx(p.I);
  const CCTK_REAL gammaDD01 = gxy(p.I);
  const CCTK_REAL gammaDD02 = gxz(p.I);
  const CCTK_REAL gammaDD11 = gyy(p.I);
  const CCTK_REAL gammaDD12 = gyz(p.I);
  const CCTK_REAL gammaDD22 = gzz(p.I);

  const CCTK_REAL my_rho = rho0(p.I);
  const CCTK_REAL lfac = w_lorentzGF(p.I);
  const CCTK_REAL my_press = pressGF(p.I);
  const CCTK_REAL my_eps = epsGF(p.I);

  const CCTK_REAL vx = velxGF(p.I);
  const CCTK_REAL vy = velyGF(p.I);
  const CCTK_REAL vz = velzGF(p.I);

  const CCTK_REAL my_lapse = alpGF(p.I);
  const CCTK_REAL my_shiftx = betaxGF(p.I);
  const CCTK_REAL my_shifty = betayGF(p.I);
  const CCTK_REAL my_shiftz = betazGF(p.I);

  const CCTK_REAL posx = p.x - cms_x;
  const CCTK_REAL posy = p.y - cms_y;
  const CCTK_REAL posz = p.z;

  const double sqrtgamma = compute_sqrtgamma(p, gxx, gxy, gxz, gyy, gyz, gzz);

  CCTK_REAL w_lorentz_limited = lfac;

  if (lfac > gamma_lim)
    w_lorentz_limited = gamma_lim;

  // Covariant 3 velocity

  const CCTK_REAL v_cov[3]{gammaDD00 * vx + gammaDD01 * vy + gammaDD02 * vz,
                           gammaDD01 * vx + gammaDD11 * vy + gammaDD12 * vz,
                           gammaDD02 * vx + gammaDD12 * vy + gammaDD22 * vz};

  const CCTK_REAL v_2 = v_cov[0] * vx + v_cov[1] * vy + v_cov[2] * vz;

  const CCTK_REAL vz_2 = v_cov[2] * vz;

  // Scalar product with shift

  const CCTK_REAL shift_v =
      v_cov[0] * my_shiftx + v_cov[1] * my_shifty + v_cov[2] * my_shiftz;

  const CCTK_REAL shiftz_vz = v_cov[2] * my_shiftz;

  // This is just the total kinetic energy

  VolIntegrand1(p.I) = 0.5 * sqrtgamma * (my_rho * (1.0 + my_eps) + my_press) *
                       w_lorentz_limited * w_lorentz_limited *
                       (v_2 - shift_v / my_lapse);

  // Cylindrical coordinates and some
  // quantities for transformation of vectors.

  const double posr2 = posx * posx + posy * posy;
  const double posr = sqrt(posr2);

  double xy_over_r2 = posx * posy / posr2;
  double x2_over_r2 = posx * posx / posr2;
  double y2_over_r2 = posy * posy / posr2;

  if (posr <= 1e-15) {

    xy_over_r2 = 0.0;
    x2_over_r2 = 0.5;
    y2_over_r2 = 0.5;
  }

  // Finally, compute v_phi*v^phi

  CCTK_REAL v_phi2{0.0};

  v_phi2 = v_cov[0] * vx * y2_over_r2 - v_cov[0] * vy * xy_over_r2 -
           v_cov[1] * vx * xy_over_r2 + v_cov[1] * vy * x2_over_r2;

  // Finally, compute v_phi*shift^phi

  CCTK_REAL v_phi_shift_phi{0.0};

  v_phi_shift_phi =
      v_cov[0] * my_shiftx * y2_over_r2 - v_cov[0] * my_shifty * xy_over_r2 -
      v_cov[1] * my_shiftx * xy_over_r2 + v_cov[1] * my_shifty * x2_over_r2;

  // This is the rotational kinetic energy wrt to the z-axis

  VolIntegrand2(p.I) = 0.5 * sqrtgamma * (my_rho * (1.0 + my_eps) + my_press) *
                       w_lorentz_limited * w_lorentz_limited *
                       (v_phi2 - v_phi_shift_phi / my_lapse);

  // Kinetic energy in z-direction

  VolIntegrand3(p.I) = 0.5 * sqrtgamma * (my_rho * (1.0 + my_eps) + my_press) *
                       w_lorentz_limited * w_lorentz_limited *
                       (vz_2 - shiftz_vz / my_lapse);
}

/* Kinetic energy density given by 0.5*rhoh*W^2v^2 */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void kinetic(
    const Loop::GF3D2<double> VolIntegrand1,
    const Loop::GF3D2<double> VolIntegrand2,
    const Loop::GF3D2<double> VolIntegrand3, const Loop::PointDesc &p,
    const Loop::GF3D2<const double> velxGF,
    const Loop::GF3D2<const double> velyGF,
    const Loop::GF3D2<const double> velzGF,
    const Loop::GF3D2<const double> w_lorentzGF,
    const Loop::GF3D2<const double> rho0, const Loop::GF3D2<const double> epsGF,
    const Loop::GF3D2<const double> pressGF,
    const Loop::GF3D2<const double> gxx, const Loop::GF3D2<const double> gxy,
    const Loop::GF3D2<const double> gxz, const Loop::GF3D2<const double> gyy,
    const Loop::GF3D2<const double> gyz, const Loop::GF3D2<const double> gzz,
    const double cms_x, const double cms_y, const double gamma_lim) {

  const CCTK_REAL gammaDD00 = gxx(p.I);
  const CCTK_REAL gammaDD01 = gxy(p.I);
  const CCTK_REAL gammaDD02 = gxz(p.I);
  const CCTK_REAL gammaDD11 = gyy(p.I);
  const CCTK_REAL gammaDD12 = gyz(p.I);
  const CCTK_REAL gammaDD22 = gzz(p.I);

  const CCTK_REAL my_rho = rho0(p.I);
  const CCTK_REAL lfac = w_lorentzGF(p.I);
  const CCTK_REAL my_press = pressGF(p.I);
  const CCTK_REAL my_eps = epsGF(p.I);

  const CCTK_REAL vx = velxGF(p.I);
  const CCTK_REAL vy = velyGF(p.I);
  const CCTK_REAL vz = velzGF(p.I);

  const CCTK_REAL posx = p.x - cms_x;
  const CCTK_REAL posy = p.y - cms_y;
  const CCTK_REAL posz = p.z;

  const double sqrtgamma = compute_sqrtgamma(p, gxx, gxy, gxz, gyy, gyz, gzz);

  CCTK_REAL w_lorentz_limited = lfac;

  if (lfac > gamma_lim)
    w_lorentz_limited = gamma_lim;

  // Covariant 3 velocity

  const CCTK_REAL v_cov[3]{gammaDD00 * vx + gammaDD01 * vy + gammaDD02 * vz,
                           gammaDD01 * vx + gammaDD11 * vy + gammaDD12 * vz,
                           gammaDD02 * vx + gammaDD12 * vy + gammaDD22 * vz};

  const CCTK_REAL v_2 = v_cov[0] * vx + v_cov[1] * vy + v_cov[2] * vz;

  const CCTK_REAL vz_2 = v_cov[2] * vz;

  // This is just the total kinetic energy

  VolIntegrand1(p.I) = 0.5 * sqrtgamma * (my_rho * (1.0 + my_eps) + my_press) *
                       (w_lorentz_limited * w_lorentz_limited * v_2);

  // Cylindrical coordinates and some
  // quantities for transformation of vectors.

  const double posr2 = posx * posx + posy * posy;
  const double posr = sqrt(posr2);

  double xy_over_r2 = posx * posy / posr2;
  double x2_over_r2 = posx * posx / posr2;
  double y2_over_r2 = posy * posy / posr2;

  if (posr <= 1e-15) {

    xy_over_r2 = 0.0;
    x2_over_r2 = 0.5;
    y2_over_r2 = 0.5;
  }

  // Finally, compute v_phi*v^phi

  CCTK_REAL v_phi2{0.0};

  v_phi2 = v_cov[0] * vx * y2_over_r2 - v_cov[0] * vy * xy_over_r2 -
           v_cov[1] * vx * xy_over_r2 + v_cov[1] * vy * x2_over_r2;

  // This is the rotational kinetic energy wrt to the z-axis

  VolIntegrand2(p.I) = 0.5 * sqrtgamma * (my_rho * (1.0 + my_eps) + my_press) *
                       w_lorentz_limited * w_lorentz_limited * v_phi2;

  // Kinetic energy in z-direction

  VolIntegrand3(p.I) = 0.5 * sqrtgamma * (my_rho * (1.0 + my_eps) + my_press) *
                       (w_lorentz_limited * w_lorentz_limited * vz_2);
}

/* Total kinetic energy as defined by total energy in hydro sector - rest mass -
 * thermal energy */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void kinetic_tot(
    const Loop::GF3D2<double> VolIntegrand1, const Loop::PointDesc &p,
    const Loop::GF3D2<const double> velxGF,
    const Loop::GF3D2<const double> velyGF,
    const Loop::GF3D2<const double> velzGF,
    const Loop::GF3D2<const double> w_lorentzGF,
    const Loop::GF3D2<const double> rho0, const Loop::GF3D2<const double> epsGF,
    const Loop::GF3D2<const double> pressGF,
    const Loop::GF3D2<const double> gxx, const Loop::GF3D2<const double> gxy,
    const Loop::GF3D2<const double> gxz, const Loop::GF3D2<const double> gyy,
    const Loop::GF3D2<const double> gyz, const Loop::GF3D2<const double> gzz,
    const double gamma_lim) {

  const CCTK_REAL gammaDD00 = gxx(p.I);
  const CCTK_REAL gammaDD01 = gxy(p.I);
  const CCTK_REAL gammaDD02 = gxz(p.I);
  const CCTK_REAL gammaDD11 = gyy(p.I);
  const CCTK_REAL gammaDD12 = gyz(p.I);
  const CCTK_REAL gammaDD22 = gzz(p.I);

  const CCTK_REAL my_rho = rho0(p.I);
  const CCTK_REAL my_press = pressGF(p.I);
  const CCTK_REAL my_eps = epsGF(p.I);

  CCTK_REAL w_lorentz = w_lorentzGF(p.I);

  const CCTK_REAL vx = velxGF(p.I);
  const CCTK_REAL vy = velyGF(p.I);
  const CCTK_REAL vz = velzGF(p.I);

  const double sqrtgamma = compute_sqrtgamma(p, gxx, gxy, gxz, gyy, gyz, gzz);

  // Covariant 3 velocity

  const CCTK_REAL v_cov[3]{gammaDD00 * vx + gammaDD01 * vy + gammaDD02 * vz,
                           gammaDD01 * vx + gammaDD11 * vy + gammaDD12 * vz,
                           gammaDD02 * vx + gammaDD12 * vy + gammaDD22 * vz};

  const CCTK_REAL v_2 = v_cov[0] * vx + v_cov[1] * vy + v_cov[2] * vz;

  const CCTK_REAL w_lorentz_2_inv = 1.0 - v_2;

  if (w_lorentz_2_inv > 0.0) {
    w_lorentz = sqrt(1. / w_lorentz_2_inv);
  }

  CCTK_REAL w_lorentz_limited = w_lorentz;

  if (w_lorentz > gamma_lim)
    w_lorentz_limited = gamma_lim;

  CCTK_REAL wm1 = w_lorentz_limited - 1.0;

  if (wm1 < 0.0)
    wm1 = 0.0;

  VolIntegrand1(p.I) =
      sqrtgamma * (my_press * (w_lorentz_limited * w_lorentz_limited * v_2) +
                   (my_rho * (1.0 + my_eps)) * w_lorentz_limited * wm1);
}

/* Thermal energy:
 * First integrand is just sqrtgamma*W*rho*eps
 * Second integrand is adiabatic heat transfer defined by the following formula
 *                  sqrtgamma*W*rho* ( eps - eps(rho,Ye,s_adiabat) )
 * s_adiabat is the entropy given by pure advection from previous to current
 * state, hence no heat transfer Important: Second integrand is not working for
 * now !!! Third integrand is sqrtgamma*W*entropy_density
 */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void thermal(
    const Loop::GF3D2<double> VolIntegrand1,
    const Loop::GF3D2<double> VolIntegrand2,
    const Loop::GF3D2<double> VolIntegrand3, const Loop::PointDesc &p,
    const Loop::GF3D2<const double> w_lorentzGF,
    const Loop::GF3D2<const double> rho0GF,
    const Loop::GF3D2<const double> epsGF,
    const Loop::GF3D2<const double> entropyGF,
    const Loop::GF3D2<const double> gxx, const Loop::GF3D2<const double> gxy,
    const Loop::GF3D2<const double> gxz, const Loop::GF3D2<const double> gyy,
    const Loop::GF3D2<const double> gyz, const Loop::GF3D2<const double> gzz,
    const double my_baryon_mass, const double gamma_lim) {

  const CCTK_REAL w_lorentz = w_lorentzGF(p.I);
  const CCTK_REAL rho0 = rho0GF(p.I);
  const CCTK_REAL eps = epsGF(p.I);
  const CCTK_REAL entropy =
      entropyGF(p.I); // this is supposed to be the entropy per baryon
                      // hence, the entropy density is given by
                      // entropy*rho0/baryon_mass for now we use neutron mass

  CCTK_REAL w_lorentz_limited = w_lorentz;
  if (w_lorentz > gamma_lim)
    w_lorentz_limited = gamma_lim;

  double sqrtgamma = compute_sqrtgamma(p, gxx, gxy, gxz, gyy, gyz, gzz);

  VolIntegrand1(p.I) = sqrtgamma * rho0 * eps * w_lorentz_limited;

  VolIntegrand2(p.I) = sqrtgamma * rho0 * w_lorentz_limited * (eps);

  // Total entropy
  VolIntegrand3(p.I) =
      sqrtgamma * rho0 * entropy * w_lorentz_limited / my_baryon_mass;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void magnetic_co(
    const Loop::GF3D2<double> VolIntegrand1, const Loop::PointDesc &p,
    const Loop::GF3D2<const double> smallb2GF,
    const Loop::GF3D2<const double> w_lorentzGF,
    const Loop::GF3D2<const double> gxx, const Loop::GF3D2<const double> gxy,
    const Loop::GF3D2<const double> gxz, const Loop::GF3D2<const double> gyy,
    const Loop::GF3D2<const double> gyz, const Loop::GF3D2<const double> gzz,
    const double gamma_lim) {

  const CCTK_REAL smallb2 = smallb2GF(p.I);
  const CCTK_REAL w_lorentz = w_lorentzGF(p.I);

  CCTK_REAL w_lorentz_limited = w_lorentz;

  if (w_lorentz > gamma_lim)
    w_lorentz_limited = gamma_lim;

  double sqrtgamma = compute_sqrtgamma(p, gxx, gxy, gxz, gyy, gyz, gzz);

  VolIntegrand1(p.I) = sqrtgamma * w_lorentz_limited * smallb2 / 2.0;
}

/* Magnetic energy split into azimuthal energy component and rest: */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void magnetic_tot(
    const Loop::GF3D2<double> VolIntegrand1,
    const Loop::GF3D2<double> VolIntegrand2, const Loop::PointDesc &p,
    const Loop::GF3D2<const double> velx, const Loop::GF3D2<const double> vely,
    const Loop::GF3D2<const double> velz, const Loop::GF3D2<const double> Bvecx,
    const Loop::GF3D2<const double> Bvecy,
    const Loop::GF3D2<const double> Bvecz, const Loop::GF3D2<const double> gxx,
    const Loop::GF3D2<const double> gxy, const Loop::GF3D2<const double> gxz,
    const Loop::GF3D2<const double> gyy, const Loop::GF3D2<const double> gyz,
    const Loop::GF3D2<const double> gzz, const double cms_x,
    const double cms_y) {

  // Notice that CoM_integrand_GAMMA_SPEED_LIMIT is applied to the integrands
  // above involving the Lorentz factor w_lorentz. The be consistent we probably
  // would have to limit the velocity vel, too.
  // This is not done so far.

  const CCTK_REAL gammaDD00 = gxx(p.I);
  const CCTK_REAL gammaDD01 = gxy(p.I);
  const CCTK_REAL gammaDD02 = gxz(p.I);
  const CCTK_REAL gammaDD11 = gyy(p.I);
  const CCTK_REAL gammaDD12 = gyz(p.I);
  const CCTK_REAL gammaDD22 = gzz(p.I);

  const CCTK_REAL B_contra[3]{Bvecx(p.I), Bvecy(p.I), Bvecz(p.I)};

  const CCTK_REAL vx = velx(p.I);
  const CCTK_REAL vy = vely(p.I);
  const CCTK_REAL vz = velz(p.I);

  CCTK_REAL posx = p.x - cms_x;
  CCTK_REAL posy = p.y - cms_y;
  CCTK_REAL posz = p.z;

  // Cylindrical coordinates and some
  // quantities for transformation of vectors.

  const double posr2 = posx * posx + posy * posy;
  const double posr = sqrt(posr2);

  double xy_over_r2 = posx * posy / posr2;
  double x2_over_r2 = posx * posx / posr2;
  double y2_over_r2 = posy * posy / posr2;

  if (posr <= 1e-15) {

    xy_over_r2 = 0.0;
    x2_over_r2 = 0.5;
    y2_over_r2 = 0.5;
  }

  const double sqrtgamma = compute_sqrtgamma(p, gxx, gxy, gxz, gyy, gyz, gzz);
  const double gamma = sqrtgamma * sqrtgamma;

  // Covariant quantities

  const CCTK_REAL v_cov[3]{gammaDD00 * vx + gammaDD01 * vy + gammaDD02 * vz,
                           gammaDD01 * vx + gammaDD11 * vy + gammaDD12 * vz,
                           gammaDD02 * vx + gammaDD12 * vy + gammaDD22 * vz};

  const CCTK_REAL B_cov[3]{gammaDD00 * B_contra[0] + gammaDD01 * B_contra[1] +
                               gammaDD02 * B_contra[2],
                           gammaDD01 * B_contra[0] + gammaDD11 * B_contra[1] +
                               gammaDD12 * B_contra[2],
                           gammaDD02 * B_contra[0] + gammaDD12 * B_contra[1] +
                               gammaDD22 * B_contra[2]};

  // Contravariant electric field * sqrtgamma

  const CCTK_REAL sqrtgammaE_contra[3]{
      B_cov[1] * v_cov[2] - B_cov[2] * v_cov[1],
      B_cov[2] * v_cov[0] - B_cov[0] * v_cov[2],
      B_cov[0] * v_cov[1] - B_cov[1] * v_cov[0]};

  // Covariant electric field * sqrtgamma

  const CCTK_REAL sqrtgammaE_cov[3]{
      gammaDD00 * sqrtgammaE_contra[0] + gammaDD01 * sqrtgammaE_contra[1] +
          gammaDD02 * sqrtgammaE_contra[2],
      gammaDD01 * sqrtgammaE_contra[0] + gammaDD11 * sqrtgammaE_contra[1] +
          gammaDD12 * sqrtgammaE_contra[2],
      gammaDD02 * sqrtgammaE_contra[0] + gammaDD12 * sqrtgammaE_contra[1] +
          gammaDD22 * sqrtgammaE_contra[2]};

  // Finally, compute E_phi*E^phi*sqrtgamma and B_phi*B_phi*sqrtgamma

  CCTK_REAL E_2{0.0};
  CCTK_REAL B_2{0.0};

  E_2 = sqrtgammaE_cov[0] * sqrtgammaE_contra[0] * y2_over_r2 -
        sqrtgammaE_cov[0] * sqrtgammaE_contra[1] * xy_over_r2 -
        sqrtgammaE_cov[1] * sqrtgammaE_contra[0] * xy_over_r2 +
        sqrtgammaE_cov[1] * sqrtgammaE_contra[1] * x2_over_r2;
  E_2 /= sqrtgamma;

  B_2 = B_cov[0] * B_contra[0] * y2_over_r2 -
        B_cov[0] * B_contra[1] * xy_over_r2 -
        B_cov[1] * B_contra[0] * xy_over_r2 +
        B_cov[1] * B_contra[1] * x2_over_r2;
  B_2 *= sqrtgamma;

  VolIntegrand1(p.I) = (E_2 + B_2) / 2.0;

  // Finally, compute E_r*E^r*sqrtgamma and B_r*B_r*sqrtgamma

  E_2 = 0.0;
  B_2 = 0.0;

  E_2 = sqrtgammaE_cov[0] * sqrtgammaE_contra[0] * x2_over_r2 +
        sqrtgammaE_cov[0] * sqrtgammaE_contra[1] * xy_over_r2 +
        sqrtgammaE_cov[1] * sqrtgammaE_contra[0] * xy_over_r2 +
        sqrtgammaE_cov[1] * sqrtgammaE_contra[1] * y2_over_r2;
  E_2 /= sqrtgamma;

  B_2 = B_cov[0] * B_contra[0] * x2_over_r2 +
        B_cov[0] * B_contra[1] * xy_over_r2 +
        B_cov[1] * B_contra[0] * xy_over_r2 +
        B_cov[1] * B_contra[1] * y2_over_r2;
  B_2 *= sqrtgamma;

  VolIntegrand2(p.I) =
      (E_2 + sqrtgammaE_cov[2] * sqrtgammaE_contra[2] / sqrtgamma + B_2 +
       B_cov[2] * B_contra[2] * sqrtgamma) /
      2.0;
}

/* Magnetic energy split into azimuthal energy component and rest: */
/* Only for specific density region: from dens_1 to dens_2 */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void magnetic_tot_12(
    const Loop::GF3D2<double> VolIntegrand1,
    const Loop::GF3D2<double> VolIntegrand2, const Loop::PointDesc &p,
    const Loop::GF3D2<const double> rho0, const CCTK_REAL dens_1,
    const CCTK_REAL dens_2, const Loop::GF3D2<const double> velx,
    const Loop::GF3D2<const double> vely, const Loop::GF3D2<const double> velz,
    const Loop::GF3D2<const double> Bvecx,
    const Loop::GF3D2<const double> Bvecy,
    const Loop::GF3D2<const double> Bvecz, const Loop::GF3D2<const double> gxx,
    const Loop::GF3D2<const double> gxy, const Loop::GF3D2<const double> gxz,
    const Loop::GF3D2<const double> gyy, const Loop::GF3D2<const double> gyz,
    const Loop::GF3D2<const double> gzz, const double cms_x,
    const double cms_y) {

  // Notice that CoM_integrand_GAMMA_SPEED_LIMIT is applied to the integrands
  // above involving the Lorentz factor w_lorentz. The be consistent we probably
  // would have to limit the velocity vel, too.
  // This is not done so far.

  const CCTK_REAL my_rho = rho0(p.I);

  if (my_rho > dens_1 && my_rho <= dens_2) {

    const CCTK_REAL gammaDD00 = gxx(p.I);
    const CCTK_REAL gammaDD01 = gxy(p.I);
    const CCTK_REAL gammaDD02 = gxz(p.I);
    const CCTK_REAL gammaDD11 = gyy(p.I);
    const CCTK_REAL gammaDD12 = gyz(p.I);
    const CCTK_REAL gammaDD22 = gzz(p.I);

    const CCTK_REAL B_contra[3]{Bvecx(p.I), Bvecy(p.I), Bvecz(p.I)};

    const CCTK_REAL vx = velx(p.I);
    const CCTK_REAL vy = vely(p.I);
    const CCTK_REAL vz = velz(p.I);

    CCTK_REAL posx = p.x - cms_x;
    CCTK_REAL posy = p.y - cms_y;
    CCTK_REAL posz = p.z;

    // Cylindrical coordinates and some
    // quantities for transformation of vectors.

    const double posr2 = posx * posx + posy * posy;
    const double posr = sqrt(posr2);

    double xy_over_r2 = posx * posy / posr2;
    double x2_over_r2 = posx * posx / posr2;
    double y2_over_r2 = posy * posy / posr2;

    if (posr <= 1e-15) {

      xy_over_r2 = 0.0;
      x2_over_r2 = 0.5;
      y2_over_r2 = 0.5;
    }

    const double sqrtgamma = compute_sqrtgamma(p, gxx, gxy, gxz, gyy, gyz, gzz);
    const double gamma = sqrtgamma * sqrtgamma;

    // Covariant quantities

    const CCTK_REAL v_cov[3]{gammaDD00 * vx + gammaDD01 * vy + gammaDD02 * vz,
                             gammaDD01 * vx + gammaDD11 * vy + gammaDD12 * vz,
                             gammaDD02 * vx + gammaDD12 * vy + gammaDD22 * vz};

    const CCTK_REAL B_cov[3]{gammaDD00 * B_contra[0] + gammaDD01 * B_contra[1] +
                                 gammaDD02 * B_contra[2],
                             gammaDD01 * B_contra[0] + gammaDD11 * B_contra[1] +
                                 gammaDD12 * B_contra[2],
                             gammaDD02 * B_contra[0] + gammaDD12 * B_contra[1] +
                                 gammaDD22 * B_contra[2]};

    // Contravariant electric field * sqrtgamma

    const CCTK_REAL sqrtgammaE_contra[3]{
        B_cov[1] * v_cov[2] - B_cov[2] * v_cov[1],
        B_cov[2] * v_cov[0] - B_cov[0] * v_cov[2],
        B_cov[0] * v_cov[1] - B_cov[1] * v_cov[0]};

    // Covariant electric field * sqrtgamma

    const CCTK_REAL sqrtgammaE_cov[3]{
        gammaDD00 * sqrtgammaE_contra[0] + gammaDD01 * sqrtgammaE_contra[1] +
            gammaDD02 * sqrtgammaE_contra[2],
        gammaDD01 * sqrtgammaE_contra[0] + gammaDD11 * sqrtgammaE_contra[1] +
            gammaDD12 * sqrtgammaE_contra[2],
        gammaDD02 * sqrtgammaE_contra[0] + gammaDD12 * sqrtgammaE_contra[1] +
            gammaDD22 * sqrtgammaE_contra[2]};

    // Finally, compute E_phi*E^phi*sqrtgamma and B_phi*B_phi*sqrtgamma

    CCTK_REAL E_2{0.0};
    CCTK_REAL B_2{0.0};

    E_2 = sqrtgammaE_cov[0] * sqrtgammaE_contra[0] * y2_over_r2 -
          sqrtgammaE_cov[0] * sqrtgammaE_contra[1] * xy_over_r2 -
          sqrtgammaE_cov[1] * sqrtgammaE_contra[0] * xy_over_r2 +
          sqrtgammaE_cov[1] * sqrtgammaE_contra[1] * x2_over_r2;
    E_2 /= sqrtgamma;

    B_2 = B_cov[0] * B_contra[0] * y2_over_r2 -
          B_cov[0] * B_contra[1] * xy_over_r2 -
          B_cov[1] * B_contra[0] * xy_over_r2 +
          B_cov[1] * B_contra[1] * x2_over_r2;
    B_2 *= sqrtgamma;

    VolIntegrand1(p.I) = (E_2 + B_2) / 2.0;

    // Finally, compute E_r*E^r*sqrtgamma and B_r*B_r*sqrtgamma

    E_2 = 0.0;
    B_2 = 0.0;

    E_2 = sqrtgammaE_cov[0] * sqrtgammaE_contra[0] * x2_over_r2 +
          sqrtgammaE_cov[0] * sqrtgammaE_contra[1] * xy_over_r2 +
          sqrtgammaE_cov[1] * sqrtgammaE_contra[0] * xy_over_r2 +
          sqrtgammaE_cov[1] * sqrtgammaE_contra[1] * y2_over_r2;
    E_2 /= sqrtgamma;

    B_2 = B_cov[0] * B_contra[0] * x2_over_r2 +
          B_cov[0] * B_contra[1] * xy_over_r2 +
          B_cov[1] * B_contra[0] * xy_over_r2 +
          B_cov[1] * B_contra[1] * y2_over_r2;
    B_2 *= sqrtgamma;

    VolIntegrand2(p.I) =
        (E_2 + sqrtgammaE_cov[2] * sqrtgammaE_contra[2] / sqrtgamma + B_2 +
         B_cov[2] * B_contra[2] * sqrtgamma) /
        2.0;

  } else {

    VolIntegrand1(p.I) = 0.0;
    VolIntegrand2(p.I) = 0.0;
  }
}
#endif
