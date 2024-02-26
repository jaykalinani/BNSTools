#ifndef INTEGRANDS_CXX
#define INTEGRANDS_CXX

/*
  NRPy+ Python code that generates below C code:
  (Note: used commit 6a86d39aa3e4c863caa51c2fc7a04f3d06727167 in
         https://github.com/zachetienne/nrpytutorial)

import sympy as sp
import grid as gri
import indexedexp as ixp
from outputC import *

DIM=3

# Need to read from memory & write to memory manually due to speed limiter if() statement.
with open("compute_rhostar.h","w") as file:
    file.write("""
// This function computes
//  rho_* = alpha*u^0*sqrt(gamma)*rho_0
//        = w_lorentz*sqrt(gamma)*rho_0 ,
// where gamma = determinant of 3-metric,
// rho_0 is the rest-mass density, and
// w_lorentz is the Lorentz factor\n\n""")

    for i in range(DIM):
        for j in range(i,DIM):
            file.write("const CCTK_REAL gammaDD"+str(i)+str(j)+" = g"+chr(ord('x')+i)+chr(ord('x')+j)+"[index];\n")
    file.write("""
const CCTK_REAL w_lorentz = w_lorentz[index];
const CCTK_REAL rho0      = rhozero[  index];

CCTK_REAL w_lorentz_limited = w_lorentz;
if(w_lorentz > CoM_integrand_GAMMA_SPEED_LIMIT) w_lorentz_limited = CoM_integrand_GAMMA_SPEED_LIMIT;
""")

    rho0,w_lorentz_limited = sp.symbols("rho0 w_lorentz_limited")
    gammaDD                = ixp.declarerank2("gammaDD", "sym01",DIM=3)

    dummy, detgamma = ixp.symm_matrix_inverter3x3(gammaDD)

    rhostar = w_lorentz_limited*sp.sqrt(detgamma)*rho0

    rho_star_str_ugly = outputC(rhostar, "const CCTK_REAL rhostar", filename="returnstring",params="includebraces=False")
    # Beautify rho_star_str_ugly
    rho_star_str = rho_star_str_ugly.replace(
        "pow(gammaDD02, 2)","(gammaDD02*gammaDD02)").replace(
        "pow(gammaDD12, 2)","(gammaDD12*gammaDD12)").replace(
        "pow(gammaDD01, 2)","(gammaDD01*gammaDD01)")

    file.write(rho_star_str)
*/

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL compute_rho_star(const int index, const CCTK_REAL *restrict w_lorentzGF, const CCTK_REAL *restrict rho0GF,
                                  const CCTK_REAL *restrict gxx,const CCTK_REAL *restrict gxy,const CCTK_REAL *restrict gxz,
                                  const CCTK_REAL *restrict gyy,const CCTK_REAL *restrict gyz,const CCTK_REAL *restrict gzz) {
  DECLARE_CCTK_PARAMETERS;

  // This function computes
  //  rho_* = alpha*u^0*sqrt(gamma)*rho_0
  //        = w_lorentz*sqrt(gamma)*rho_0 ,
  // where gamma = determinant of 3-metric,
  // rho_0 is the rest-mass density, and
  // w_lorentz is the Lorentz factor

  const CCTK_REAL gammaDD00 = gxx[index];
  const CCTK_REAL gammaDD01 = gxy[index];
  const CCTK_REAL gammaDD02 = gxz[index];
  const CCTK_REAL gammaDD11 = gyy[index];
  const CCTK_REAL gammaDD12 = gyz[index];
  const CCTK_REAL gammaDD22 = gzz[index];

  const CCTK_REAL w_lorentz = w_lorentzGF[index];
  const CCTK_REAL rho0      = rho0GF[     index];

  CCTK_REAL w_lorentz_limited = w_lorentz;
  if(w_lorentz > CoM_integrand_GAMMA_SPEED_LIMIT) w_lorentz_limited = CoM_integrand_GAMMA_SPEED_LIMIT;
  /*
   *  Original SymPy expression:
   *  "rhostar = rho0*w_lorentz_limited*sqrt(gammaDD00*gammaDD11*gammaDD22 - gammaDD00*gammaDD12**2 - gammaDD01**2*gammaDD22 + 2*gammaDD01*gammaDD02*gammaDD12 - gammaDD02**2*gammaDD11)"
   */
  const CCTK_REAL rhostar = rho0*w_lorentz_limited*sqrt(gammaDD00*gammaDD11*gammaDD22 - gammaDD00*(gammaDD12*gammaDD12) - (gammaDD01*gammaDD01)*gammaDD22 + 2*gammaDD01*gammaDD02*gammaDD12 - (gammaDD02*gammaDD02)*gammaDD11);

  return rhostar;
}


/* Center of Mass: */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void CoM_integrand(double *VolIntegrand1,double *VolIntegrand2,double *VolIntegrand3,double *VolIntegrand4, const int index,
                          const CCTK_REAL *restrict w_lorentz, const CCTK_REAL *restrict rho0,
                          const CCTK_REAL *restrict gxx,const CCTK_REAL *restrict gxy,const CCTK_REAL *restrict gxz,
                          const CCTK_REAL *restrict gyy,const CCTK_REAL *restrict gyz,const CCTK_REAL *restrict gzz,
                          const CCTK_REAL *restrict x,const CCTK_REAL *restrict y,const CCTK_REAL *restrict z) {
  double rho_starL = compute_rho_star(index, w_lorentz, rho0, gxx,gxy,gxz, gyy,gyz,gzz);
  VolIntegrand1[index] = rho_starL*x[index];
  VolIntegrand2[index] = rho_starL*y[index];
  VolIntegrand3[index] = rho_starL*z[index];
  VolIntegrand4[index] = rho_starL;
}

/* Rest Mass: */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void M0_integrand(double *VolIntegrand1, const int index,
                          const CCTK_REAL *restrict w_lorentz, const CCTK_REAL *restrict rho0,
                          const CCTK_REAL *restrict gxx,const CCTK_REAL *restrict gxy,const CCTK_REAL *restrict gxz,
                          const CCTK_REAL *restrict gyy,const CCTK_REAL *restrict gyz,const CCTK_REAL *restrict gzz) {
  double rho_starL = compute_rho_star(index, w_lorentz, rho0, gxx,gxy,gxz, gyy,gyz,gzz);
  VolIntegrand1[index] = rho_starL;
}

#endif
