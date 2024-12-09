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
    const Loop::PointDesc &p, const Loop::GF3D2<const double> velx,
    const Loop::GF3D2<const double> vely, const Loop::GF3D2<const double> velz,
    const Loop::GF3D2<const double> rho0GF, const Loop::GF3D2<const double> gxx,
    const Loop::GF3D2<const double> gxy, const Loop::GF3D2<const double> gxz,
    const Loop::GF3D2<const double> gyy, const Loop::GF3D2<const double> gyz,
    const Loop::GF3D2<const double> gzz, const CCTK_REAL gamma_speed_limit) {
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

  const CCTK_REAL vU0 = velx(p.I);
  const CCTK_REAL vU1 = vely(p.I);
  const CCTK_REAL vU2 = velz(p.I);

  const CCTK_REAL vD0 = gammaDD00 * vU0 + gammaDD01 * vU1 + gammaDD02 * vU2;
  const CCTK_REAL vD1 = gammaDD01 * vU0 + gammaDD11 * vU1 + gammaDD12 * vU2;
  const CCTK_REAL vD2 = gammaDD02 * vU0 + gammaDD12 * vU1 + gammaDD22 * vU2;

  const CCTK_REAL w_lorentz =
      1.0 / sqrt(1.0 - (vU0 * vD0 + vU1 * vD1 + vU2 * vD2));
  const CCTK_REAL rho0 = rho0GF(p.I);

  CCTK_REAL w_lorentz_limited = w_lorentz;
  if (w_lorentz > gamma_speed_limit)
    w_lorentz_limited = gamma_speed_limit;
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

/* Center of Mass: */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void CoM_integrand(
    const Loop::GF3D2<double> VolIntegrand1,
    const Loop::GF3D2<double> VolIntegrand2,
    const Loop::GF3D2<double> VolIntegrand3,
    const Loop::GF3D2<double> VolIntegrand4, const Loop::PointDesc &p,
    const Loop::GF3D2<const double> velx, const Loop::GF3D2<const double> vely,
    const Loop::GF3D2<const double> velz, const Loop::GF3D2<const double> rho0,
    const Loop::GF3D2<const double> gxx, const Loop::GF3D2<const double> gxy,
    const Loop::GF3D2<const double> gxz, const Loop::GF3D2<const double> gyy,
    const Loop::GF3D2<const double> gyz, const Loop::GF3D2<const double> gzz,
    const CCTK_REAL gamma_speed_limit) {
  double rho_starL = compute_rho_star(p, velx, vely, velz, rho0, gxx, gxy, gxz,
                                      gyy, gyz, gzz, gamma_speed_limit);
  VolIntegrand1(p.I) = rho_starL * p.x;
  VolIntegrand2(p.I) = rho_starL * p.y;
  VolIntegrand3(p.I) = rho_starL * p.z;
  VolIntegrand4(p.I) = rho_starL;
}

/* Rest Mass: */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void M0_integrand(
    const Loop::GF3D2<double> VolIntegrand1, const Loop::PointDesc &p,
    const Loop::GF3D2<const double> velx, const Loop::GF3D2<const double> vely,
    const Loop::GF3D2<const double> velz, const Loop::GF3D2<const double> rho0,
    const Loop::GF3D2<const double> gxx, const Loop::GF3D2<const double> gxy,
    const Loop::GF3D2<const double> gxz, const Loop::GF3D2<const double> gyy,
    const Loop::GF3D2<const double> gyz, const Loop::GF3D2<const double> gzz,
    const CCTK_REAL gamma_speed_limit) {
  double rho_starL = compute_rho_star(p, velx, vely, velz, rho0, gxx, gxy, gxz,
                                      gyy, gyz, gzz, gamma_speed_limit);
  VolIntegrand1(p.I) = rho_starL;
}

#endif
