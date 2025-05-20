#include "cctk.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include <loop_device.hxx>

#include "aster_utils.hxx"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define COMPUTE_DERIV_4(METRICm2,METRICm1,METRICp1,METRICp2) ((1./12.)*(METRICm2) + (-8./12.)*(METRICm1) + (8./12.)*(METRICp1) + (-1./12.)*(METRICp2))

namespace MHD_Diagnostics {

using namespace Loop;
using namespace AsterUtils;

template<typename T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T SQR(T t) { return t*t; }

extern "C" void compute_mhd_diagnostics(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_compute_mhd_diagnostics;
  DECLARE_CCTK_PARAMETERS;

  if(diagnostics_compute_every<=0 || cctk_iteration%diagnostics_compute_every!=0) return;

  const double ONE_OVER_SQRT_4PI = 1/sqrt(4*M_PI);
  double dX[3] = { CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1), CCTK_DELTA_SPACE(2) };
  double dxi[4] = { 1e100,1.0/dX[0],1.0/dX[1],1.0/dX[2] };
        
  double current_time      = cctkGH->cctk_time;
  double local_stepsize    = cctkGH->cctk_delta_time / cctkGH->cctk_timefac;
  double time_at_p         = current_time - local_stepsize;

  /* Gridfunctions */
  const vec<GF3D2<const CCTK_REAL>, dim> gf_beta{betax, betay, betaz};
  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

	// Load metric
  /* Get covariant metric */
  const smat<CCTK_REAL, 3> glo(
    [&](int i, int j) ARITH_INLINE { return calc_avg_v2c(gf_g(i, j), p); });

  CCTK_REAL gxxL = glo(0, 0);
  CCTK_REAL gxyL = glo(0, 1);
  CCTK_REAL gxzL = glo(0, 2);
  CCTK_REAL gyyL = glo(1, 1);
  CCTK_REAL gyzL = glo(1, 2);
  CCTK_REAL gzzL = glo(2, 2);

	// Further geometric quantities 
	double det = -gxzL*gxzL*gyyL + 2*gxyL*gxzL*gyzL - gxxL*gyzL*gyzL - gxyL*gxyL*gzzL + gxxL*gyyL*gzzL;
  sqrt_gamma(p.I) = sqrt(det);

	double invdet = 1.0 / det;
	double gupxxL=(-gyzL*gyzL + gyyL*gzzL)*invdet;
	double gupxyL=( gxzL*gyzL - gxyL*gzzL)*invdet;
	double gupyyL=(-gxzL*gxzL + gxxL*gzzL)*invdet;
	double gupxzL=(-gxzL*gyyL + gxyL*gyzL)*invdet;
	double gupyzL=( gxyL*gxzL - gxxL*gyzL)*invdet;
	double gupzzL=(-gxyL*gxyL + gxxL*gyyL)*invdet;

  double lapse = calc_avg_v2c(alp, p);

  const vec<CCTK_REAL, 3> betas_avg(
      [&](int i) ARITH_INLINE { return calc_avg_v2c(gf_beta(i), p); });
	double shiftx = betas_avg(0);
	double shifty = betas_avg(1);
	double shiftz = betas_avg(2);

	double ONE_OVER_LAPSE = 1.0/lapse;
	double ONE_OVER_LAPSE_SQRT_4PI = ONE_OVER_LAPSE*ONE_OVER_SQRT_4PI;
	double ONE_OVER_ALP_SQRTGAMMA_4PI = ONE_OVER_LAPSE*sqrt(invdet)*ONE_OVER_SQRT_4PI*ONE_OVER_SQRT_4PI; 

	// Fluid 3 velocity in normal frame and Lorentz factor
	double ETvx = velx(p.I);
	double ETvy = vely(p.I);
	double ETvz = velz(p.I);
	double lfac = w_lorentz(p.I);

// ===========================================
//
//
// Etienne's code for smallbs and Poynting: Start
//
//
// ===========================================

	// IllinoisGRMHD defines v^i = u^i/u^0.
        
	// Meanwhile, the ET/HydroBase formalism, called the Valencia 
	// formalism, splits the 4 velocity into a purely spatial part
	// and a part that is normal to the spatial hypersurface:
	// u^a = G (n^a + U^a), (Eq. 14 of arXiv:1304.5544; G=W, U^a=v^a)
	// where n^a is the unit normal vector to the spatial hypersurface,
	// n_a = {-\alpha,0,0,0}, and U^a is the purely spatial part, which
	// is defined in HydroBase as the vel[] vector gridfunction.
	// Then u^a n_a = - \alpha u^0 = G n^a n_a = -G, and
	// of course \alpha u^0 = 1/sqrt(1+γ^ij u_i u_j) = \Gamma,
	// the standard Lorentz factor.

	// Note that n^i = - \beta^i / \alpha, so 
	// u^a = \Gamma (n^a + U^a) 
	// -> u^i = \Gamma ( U^i - \beta^i / \alpha )
	// which implies
	// v^i = u^i/u^0
	//     = \Gamma/u^0 ( U^i - \beta^i / \alpha ) <- \Gamma = \alpha u^0
	//     = \alpha ( U^i - \beta^i / \alpha )
	//     = \alpha U^i - \beta^i

	double vxL = lapse*ETvx - shiftx;
	double vyL = lapse*ETvy - shifty;
	double vzL = lapse*ETvz - shiftz;

	// Derivation of first equation:
	// \gamma_{ij} (v^i + \beta^i)(v^j + \beta^j)/(\alpha)^2 
	//   = \gamma_{ij} 1/(u^0)^2 ( \gamma^{ik} u_k \gamma^{jl} u_l /(\alpha)^2 <- Using Eq. 53 of arXiv:astro-ph/0503420
	//   = 1/(u^0 \alpha)^2 u_j u_l \gamma^{jl}  <- Since \gamma_{ij} \gamma^{ik} = \delta^k_j
	//   = 1/(u^0 \alpha)^2 ( (u^0 \alpha)^2 - 1 ) <- Using Eq. 56 of arXiv:astro-ph/0503420
	//   = 1 - 1/(u^0 \alpha)^2 <= 1
	double one_minus_one_over_alpha_u0_squared = (gxxL* SQR(vxL + shiftx) +
						      2.0*gxyL*(vxL + shiftx)*(vyL + shifty) +
						      2.0*gxzL*(vxL + shiftx)*(vzL + shiftz) +
						      gyyL* SQR(vyL + shifty) +
						      2.0*gyzL*(vyL + shifty)*(vzL + shiftz) +
						      gzzL* SQR(vzL + shiftz) )*SQR(ONE_OVER_LAPSE);
	/*** Check for superluminal velocity ***/
	/* Don't bother with speed limits in this analysis thorn.
	//FIXME: Instead of >1.0, should be one_minus_one_over_alpha_u0_squared > ONE_MINUS_ONE_OVER_GAMMA_SPEED_LIMIT_SQUARED, for consistency with conserv_to_prims routines
	if(one_minus_one_over_alpha_u0_squared > 1.0) {
	double ONE_MINUS_ONE_OVER_GAMMA_SPEED_LIMIT_SQUARED = 1.0-1.0/SQR(GAMMA_SPEED_LIMIT);
	double correction_fac = sqrt(ONE_MINUS_ONE_OVER_GAMMA_SPEED_LIMIT_SQUARED/one_minus_one_over_alpha_u0_squared);
	vxL = (vxL + shiftx)*correction_fac-shiftx;
	vyL = (vyL + shifty)*correction_fac-shifty;
	vzL = (vzL + shiftz)*correction_fac-shiftz;
	one_minus_one_over_alpha_u0_squared=ONE_MINUS_ONE_OVER_GAMMA_SPEED_LIMIT_SQUARED;
	stats.failure_checker+=1000;
	}
	*/
	// A = 1.0-one_minus_one_over_alpha_u0_squared = 1-(1-1/(al u0)^2) = 1/(al u0)^2
	// 1/sqrt(A) = al u0
	//double alpha_u0_minus_one = 1.0/sqrt(1.0-one_minus_one_over_alpha_u0_squared)-1.0;
	//u0_out          = (alpha_u0_minus_one + 1.0)*ONE_OVER_LAPSE;
	double alpha_u0 = 1.0/sqrt(1.0-one_minus_one_over_alpha_u0_squared);
	double u0L = alpha_u0*ONE_OVER_LAPSE;
 
	double Bx_center = Bvecx(p.I);
	double By_center = Bvecy(p.I);
	double Bz_center = Bvecz(p.I);

	// NOW COMPUTE b^{\mu} and b^2 = b^{\mu} b^{\nu} g_{\mu \nu}
	double ONE_OVER_U0 = 1.0/u0L;
	double shiftx_plus_vx = (shiftx+vxL);
	double shifty_plus_vy = (shifty+vyL);
	double shiftz_plus_vz = (shiftz+vzL);

	// Eq. 56 in http://arxiv.org/pdf/astro-ph/0503420.pdf:
	//  u_i = gamma_{ij} u^0 (v^j + beta^j), gamma_{ij} is the physical metric
	double u_x_over_u0 =  gxxL*shiftx_plus_vx + gxyL*shifty_plus_vy + gxzL*shiftz_plus_vz;
	double u_y_over_u0 =  gxyL*shiftx_plus_vx + gyyL*shifty_plus_vy + gyzL*shiftz_plus_vz;
	double u_z_over_u0 =  gxzL*shiftx_plus_vx + gyzL*shifty_plus_vy + gzzL*shiftz_plus_vz;

	// Eqs. 23 and 31 in http://arxiv.org/pdf/astro-ph/0503420.pdf:
	//   Compute alpha sqrt(4 pi) b^t = u_i B^i
	double alpha_sqrt_4pi_bt = ( u_x_over_u0*Bx_center + u_y_over_u0*By_center + u_z_over_u0*Bz_center ) * u0L;
	// Eq. 24 in http://arxiv.org/pdf/astro-ph/0503420.pdf:
	// b^i = B^i_u / sqrt(4 pi)
	// b^i = ( B^i/alpha + B^0_u u^i ) / ( u^0 sqrt(4 pi) )
	// b^i = ( B^i/alpha +  sqrt(4 pi) b^t u^i ) / ( u^0 sqrt(4 pi) )
	// b^i = ( B^i +  alpha sqrt(4 pi) b^t u^i ) / ( alpha u^0 sqrt(4 pi) )
	// b^i = ( B^i/u^0 +  alpha sqrt(4 pi) b^t u^i/u^0 ) / ( alpha sqrt(4 pi) )
	// b^i = ( B^i/u^0 +  alpha sqrt(4 pi) b^t v^i ) / ( alpha sqrt(4 pi) )
	double smallbxL = (Bx_center*ONE_OVER_U0 + vxL*alpha_sqrt_4pi_bt)*ONE_OVER_LAPSE_SQRT_4PI;
	double smallbyL = (By_center*ONE_OVER_U0 + vyL*alpha_sqrt_4pi_bt)*ONE_OVER_LAPSE_SQRT_4PI;
	double smallbzL = (Bz_center*ONE_OVER_U0 + vzL*alpha_sqrt_4pi_bt)*ONE_OVER_LAPSE_SQRT_4PI;
	// Eq. 23 in http://arxiv.org/pdf/astro-ph/0503420.pdf, with alpha sqrt (4 pi) b^2 = u_i B^i already computed above
	double smallbtL = alpha_sqrt_4pi_bt * ONE_OVER_LAPSE_SQRT_4PI;

	/* POYNTING FLUX: */
	// S^i = -alp*Tem^i_0 = -alp*(b^2 u^i u_0 + b^2/2 g^i_0 - b^i b_0) 
	double LAPSE_SQUARED = lapse*lapse;
	double ONE_OVER_LAPSE_SQUARED = ONE_OVER_LAPSE*ONE_OVER_LAPSE;
	double g4dn[4][4],g4up[4][4];
    
	// g^{\mu \nu} = upper four-metric.
	g4up[0][0]              = -ONE_OVER_LAPSE_SQUARED;
	g4up[0][1] = g4up[1][0] = ONE_OVER_LAPSE_SQUARED*shiftx;
	g4up[0][2] = g4up[2][0] = ONE_OVER_LAPSE_SQUARED*shifty;
	g4up[0][3] = g4up[3][0] = ONE_OVER_LAPSE_SQUARED*shiftz;
	g4up[1][1]              = gupxxL - ONE_OVER_LAPSE_SQUARED*shiftx*shiftx;
	g4up[1][2] = g4up[2][1] = gupxyL - ONE_OVER_LAPSE_SQUARED*shiftx*shifty;
	g4up[1][3] = g4up[3][1] = gupxzL - ONE_OVER_LAPSE_SQUARED*shiftx*shiftz;
	g4up[2][2]              = gupyyL - ONE_OVER_LAPSE_SQUARED*shifty*shifty;
	g4up[2][3] = g4up[3][2] = gupyzL - ONE_OVER_LAPSE_SQUARED*shifty*shiftz;
	g4up[3][3]              = gupzzL - ONE_OVER_LAPSE_SQUARED*shiftz*shiftz;

	/* Compute beta_i */
	double shift_x = ( gxxL*(shiftx) + gxyL*(shifty) +
			   gxzL*(shiftz) );
	double shift_y = ( gxyL*(shiftx) + gyyL*(shifty) +
			   gyzL*(shiftz) );
	double shift_z = ( gxzL*(shiftx) + gyzL*(shifty) +
			   gzzL*(shiftz) );
	
	// g_{00}               = - alpha^2 + gamma_{ij} beta^i beta^j = - alpha^2 beta_i beta^i
	g4dn[0][0]              = -LAPSE_SQUARED + (shiftx*shift_x + shifty*shift_y + shiftz*shift_z);
	// g_{0i} =  gamma_{ij} beta^j = beta_i
	g4dn[0][1] = g4dn[1][0] = shift_x;
	g4dn[0][2] = g4dn[2][0] = shift_y;
	g4dn[0][3] = g4dn[3][0] = shift_z;
	// g_{ij} =  gamma_{ij} <- 3 metric
	g4dn[1][1] =              gxxL;
	g4dn[1][2] = g4dn[2][1] = gxyL;
	g4dn[1][3] = g4dn[3][1] = gxzL;
	g4dn[2][2] =              gyyL;
	g4dn[2][3] = g4dn[3][2] = gyzL;
	g4dn[3][3] =              gzzL;

	// S^i = -alp*Tem^i_0 = -alp*(b^2 u^i u_0 + b^2/2 g^i_0 - b^i b_0) 

	// First compute u_0.
	double uup[4] = {u0L,u0L*vxL,u0L*vyL,u0L*vzL};
	double u_0=0.0; for(int ii=0;ii<4;ii++) u_0 += g4dn[0][ii]*uup[ii];

	// Next compute b_0.
	double sbup[4] = {smallbtL,smallbxL,smallbyL,smallbzL};
	double smallb_0L=0.0; for(int ii=0;ii<4;ii++) smallb_0L += g4dn[0][ii]*sbup[ii];

	// Next compute g^i_0:
	double g4upx_0=0.0; for(int ii=0;ii<4;ii++) g4upx_0 += g4up[1][ii]*g4dn[0][ii];
	double g4upy_0=0.0; for(int ii=0;ii<4;ii++) g4upy_0 += g4up[2][ii]*g4dn[0][ii];
	double g4upz_0=0.0; for(int ii=0;ii<4;ii++) g4upz_0 += g4up[3][ii]*g4dn[0][ii];

	// Next compute b^2:
	double smallb[4] = {smallbtL,smallbxL,smallbyL,smallbzL};
	double smallb2L=0.0; for(int ii=0;ii<4;ii++) for(int jj=0;jj<4;jj++) smallb2L += g4dn[ii][jj]*smallb[ii]*smallb[jj];

	// S^i = -alp*Tem^i_0 = -alp*(b^2 u^i u_0 + b^2/2 g^i_0 - b^i b_0)
	// Poyn2^i = -alp*Tem^i^0 = -alp*(b^2 u^i u^0 + b^2/2 g^i^0 - b^i b^0)
	double PoynxL = -lapse*(smallb2L*(u0L*vxL)*u_0 + 0.5*smallb2L*g4upx_0 - smallbxL*smallb_0L);
	double PoynyL = -lapse*(smallb2L*(u0L*vyL)*u_0 + 0.5*smallb2L*g4upy_0 - smallbyL*smallb_0L);
	double PoynzL = -lapse*(smallb2L*(u0L*vzL)*u_0 + 0.5*smallb2L*g4upz_0 - smallbzL*smallb_0L);

	double Poyn2xL = -lapse*(smallb2L*(u0L*vxL)*u0L + 0.5*smallb2L*g4up[1][0] - smallbxL*smallbtL);
	double Poyn2yL = -lapse*(smallb2L*(u0L*vyL)*u0L + 0.5*smallb2L*g4up[2][0] - smallbyL*smallbtL);
	double Poyn2zL = -lapse*(smallb2L*(u0L*vzL)*u0L + 0.5*smallb2L*g4up[3][0] - smallbzL*smallbtL);

	minus_one_minus_u_0(p.I) = -1.0 - u_0;

	smallbt(p.I) = smallbtL;
	smallbx(p.I) = smallbxL;
	smallby(p.I) = smallbyL;
	smallbz(p.I) = smallbzL;
	smallb2(p.I) = smallb2L;

	Poynx(p.I) = PoynxL;
	Poyny(p.I) = PoynyL;
	Poynz(p.I) = PoynzL;

	Poyn2x(p.I) = Poyn2xL;
	Poyn2y(p.I) = Poyn2yL;
	Poyn2z(p.I) = Poyn2zL;

// ===========================================
//
//
// Etienne's code for smallbs and Poynting: End
//
//
// ===========================================

// Before turning to derivatives, calculate norm B

        double Bcov[3] = {gxxL*Bx_center + gxyL*By_center + gxzL*Bz_center,
	                  gxyL*Bx_center + gyyL*By_center + gyzL*Bz_center,
                          gxzL*Bx_center + gyzL*By_center + gzzL*Bz_center};

	double B2 = Bcov[0]*Bx_center + Bcov[1]*By_center + Bcov[2]*Bz_center;

	normB(p.I) = sqrt(B2);
  });
}

extern "C" void compute_mhd_derivs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_compute_mhd_derivs;
  DECLARE_CCTK_PARAMETERS;

  if(diagnostics_compute_every<=0 || cctk_iteration%diagnostics_compute_every!=0) return;

  const double ONE_OVER_SQRT_4PI = 1/sqrt(4*M_PI);
  double dX[3] = { CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1), CCTK_DELTA_SPACE(2) };
  double dxi[4] = { 1e100,1.0/dX[0],1.0/dX[1],1.0/dX[2] };
        
  double current_time      = cctkGH->cctk_time;
  double local_stepsize    = cctkGH->cctk_delta_time / cctkGH->cctk_timefac;
  double time_at_p         = current_time - local_stepsize;

  constexpr unsigned kronecker_delta[4][4] = { { 0,0,0,0 },
                                               { 0,1,0,0 },
                                               { 0,0,1,0 },
                                               { 0,0,0,1 } };

  /* Gridfunctions */
  const vec<GF3D2<const CCTK_REAL>, dim> gf_beta{betax, betay, betaz};
  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};
  const smat<GF3D2<const CCTK_REAL>, 3> gf_k{kxx, kxy, kxz, kyy, kyz, kzz};

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

	double ETvx = velx(p.I);
	double ETvy = vely(p.I);
	double ETvz = velz(p.I);
	double lfac   = w_lorentz(p.I);

// Calculate also divB from AsterX staggered field
// Use fourth-order ECHO expressions
// TODO: check if face index is correct

	double Bx_stagg = 13.*Bx_stag(p.I + p.DI[0])/12. - (Bx_stag(p.I)+Bx_stag(p.I + 2 * p.DI[0]))/24.;
	double By_stagg = 13.*By_stag(p.I + p.DI[1])/12. - (By_stag(p.I)+By_stag(p.I + 2 * p.DI[1]))/24.;
	double Bz_stagg = 13.*Bz_stag(p.I + p.DI[2])/12. - (Bz_stag(p.I)+Bz_stag(p.I + 2 * p.DI[2]))/24.;

	double Bx_stagg_m1 = 13.*Bx_stag(p.I)/12. - (Bx_stag(p.I + p.DI[0])+Bx_stag(p.I - p.DI[0]))/24.;
	double By_stagg_m1 = 13.*By_stag(p.I)/12. - (By_stag(p.I + p.DI[1])+By_stag(p.I - p.DI[1]))/24.;
	double Bz_stagg_m1 = 13.*Bz_stag(p.I)/12. - (Bz_stag(p.I + p.DI[2])+Bz_stag(p.I - p.DI[2]))/24.;

	divB(p.I) = dxi[1]*(Bx_stagg-Bx_stagg_m1) + dxi[2]*(By_stagg-By_stagg_m1) + dxi[3]*(Bz_stagg-Bz_stagg_m1);

// ===========================================

	double my_rho = rho(p.I);

	if(derivatives_compute_every<=0 || cctk_iteration%derivatives_compute_every!=0 || my_rho<rho_cutoff){

	expansion_scalar(p.I)               = 0.0;
	expansion_no_time_der(p.I)               = 0.0;
	shear_spatial_tensor_xx(p.I)        = 0.0;
	shear_spatial_tensor_xy(p.I)        = 0.0;
	shear_spatial_tensor_xz(p.I)        = 0.0;
	shear_spatial_tensor_yy(p.I)        = 0.0;
	shear_spatial_tensor_yz(p.I)        = 0.0;
	shear_spatial_tensor_zz(p.I)        = 0.0;
	kin_vorticity_spatial_xy(p.I)       = 0.0;
	kin_vorticity_spatial_xz(p.I)       = 0.0;
	kin_vorticity_spatial_yz(p.I)       = 0.0;
	kin_acceleration_spatial_x(p.I)     = 0.0;
	kin_acceleration_spatial_y(p.I)     = 0.0;
        kin_acceleration_spatial_z(p.I)     = 0.0;

        sigma4bb(p.I) = 0.0;
	sigma4Ut(p.I) = 0.0;
	sigma4Ux(p.I) = 0.0;
	sigma4Uy(p.I) = 0.0;
	sigma4Uz(p.I) = 0.0;
	sigma4Trace(p.I) = 0.0;
	omega4Ut(p.I) = 0.0;
	omega4Ux(p.I) = 0.0;
	omega4Uy(p.I) = 0.0;
	omega4Uz(p.I) = 0.0;

	normcurlB(p.I) = 0.0;
	a4sq(p.I) = 0.0;

	} else {

	// ---------------------------------
	// Fluid derivative code starts here
        // ---------------------------------

	// Load metric
  /* Get covariant metric */
  const smat<CCTK_REAL, 3> glo(
    [&](int i, int j) ARITH_INLINE { return calc_avg_v2c(gf_g(i, j), p); });

  CCTK_REAL gxxL = glo(0, 0);
  CCTK_REAL gxyL = glo(0, 1);
  CCTK_REAL gxzL = glo(0, 2);
  CCTK_REAL gyyL = glo(1, 1);
  CCTK_REAL gyzL = glo(1, 2);
  CCTK_REAL gzzL = glo(2, 2);

  double lapse = calc_avg_v2c(alp, p);
  const vec<CCTK_REAL, 3> betas_avg(
      [&](int i) ARITH_INLINE { return calc_avg_v2c(gf_beta(i), p); });
	double shiftx = betas_avg(0);
	double shifty = betas_avg(1);
	double shiftz = betas_avg(2);

	// Further geometric quantities -- Much of this is repeated from the previous routine, but doing it this way minimizes syncs
	// while properly filling ghost zones?
	double sqrtdet = sqrt_gamma(p.I);
	double det = sqrtdet * sqrtdet;

	double invdet = 1.0 / det;
	double gupxxL=(-gyzL*gyzL + gyyL*gzzL)*invdet;
	double gupxyL=( gxzL*gyzL - gxyL*gzzL)*invdet;
	double gupyyL=(-gxzL*gxzL + gxxL*gzzL)*invdet;
	double gupxzL=(-gxzL*gyyL + gxyL*gyzL)*invdet;
	double gupyzL=( gxyL*gxzL - gxxL*gyzL)*invdet;
	double gupzzL=(-gxyL*gxyL + gxxL*gyyL)*invdet;

	double LAPSE_SQUARED = lapse*lapse;
	double ONE_OVER_LAPSE = 1.0/lapse;
	double ONE_OVER_ALP_SQRTGAMMA_4PI = ONE_OVER_LAPSE*sqrt(invdet)*ONE_OVER_SQRT_4PI*ONE_OVER_SQRT_4PI; 

	/* Compute beta_i */
	double shift_x = ( gxxL*(shiftx) + gxyL*(shifty) +
			   gxzL*(shiftz) );
	double shift_y = ( gxyL*(shiftx) + gyyL*(shifty) +
			   gyzL*(shiftz) );
	double shift_z = ( gxzL*(shiftx) + gyzL*(shifty) +
			   gzzL*(shiftz) );

  double g4dn[4][4];
	// g_{00}               = - alpha^2 + gamma_{ij} beta^i beta^j = - alpha^2 beta_i beta^i
	g4dn[0][0]              = -LAPSE_SQUARED + (shiftx*shift_x + shifty*shift_y + shiftz*shift_z);
	// g_{0i} =  gamma_{ij} beta^j = beta_i
	g4dn[0][1] = g4dn[1][0] = shift_x;
	g4dn[0][2] = g4dn[2][0] = shift_y;
	g4dn[0][3] = g4dn[3][0] = shift_z;
	// g_{ij} =  gamma_{ij} <- 3 metric
	g4dn[1][1] =              gxxL;
	g4dn[1][2] = g4dn[2][1] = gxyL;
	g4dn[1][3] = g4dn[3][1] = gxzL;
	g4dn[2][2] =              gyyL;
	g4dn[2][3] = g4dn[3][2] = gyzL;
	g4dn[3][3] =              gzzL;

	double ETvx_p = velxold(p.I);
	double ETvy_p = velyold(p.I);
	double ETvz_p = velzold(p.I);

	double lfac_p = w_lorentzold(p.I);

  double vxL = lapse*ETvx - shiftx;
  double vyL = lapse*ETvy - shifty;
  double vzL = lapse*ETvz - shiftz;

	double one_minus_one_over_alpha_u0_squared = (gxxL* SQR(vxL + shiftx) +
						      2.0*gxyL*(vxL + shiftx)*(vyL + shifty) +
						      2.0*gxzL*(vxL + shiftx)*(vzL + shiftz) +
						      gyyL* SQR(vyL + shifty) +
						      2.0*gyzL*(vyL + shifty)*(vzL + shiftz) +
						      gzzL* SQR(vzL + shiftz) )*SQR(ONE_OVER_LAPSE);
	double alpha_u0 = 1.0/sqrt(1.0-one_minus_one_over_alpha_u0_squared);
	double u0L = alpha_u0*ONE_OVER_LAPSE;
	double uup[4] = {u0L,u0L*vxL,u0L*vyL,u0L*vzL};
	double smallb[4] = {smallbt(p.I),smallbx(p.I),smallby(p.I),smallbz(p.I)};

	// Calculate temporal partial derivatives

	double partial_t_Wv_xU = (lfac*ETvx-lfac_p*ETvx_p)/local_stepsize;
	double partial_t_Wv_yU = (lfac*ETvy-lfac_p*ETvy_p)/local_stepsize;
	double partial_t_Wv_zU = (lfac*ETvz-lfac_p*ETvz_p)/local_stepsize;

	double partial_t_W = (lfac-lfac_p)/local_stepsize;

	// Load extrinsic curvature
  const smat<CCTK_REAL, 3> klo(
    [&](int i, int j) ARITH_INLINE { return calc_avg_v2c(gf_k(i, j), p); });

	double kxxL = klo(0, 0);
	double kxyL = klo(0, 1);
	double kxzL = klo(0, 2);
	double kyyL = klo(1, 1);
	double kyzL = klo(1, 2);
	double kzzL = klo(2, 2);

	// Get extrinsic curvature with contravariant indices
	double k_xU_xD = gupxxL*kxxL+gupxyL*kxyL+gupxzL*kxzL;
        double k_xU_yD = gupxxL*kxyL+gupxyL*kyyL+gupxzL*kyzL;
        double k_xU_zD = gupxxL*kxzL+gupxyL*kyzL+gupxzL*kzzL;
        double k_yU_yD = gupxyL*kxyL+gupyyL*kyyL+gupyzL*kyzL;
        double k_yU_zD = gupxyL*kxzL+gupyyL*kyzL+gupyzL*kzzL;
        double k_zU_zD = gupxzL*kxzL+gupyzL*kyzL+gupzzL*kzzL;

        double k_xD_yU = gupxyL*kxxL+gupyyL*kxyL+gupyzL*kxzL;
        double k_xD_zU = gupxzL*kxxL+gupyzL*kxyL+gupzzL*kxzL;
        double k_yD_zU = gupxzL*kxyL+gupyzL*kyyL+gupzzL*kyzL;

	double kupxxL = gupxxL*k_xU_xD+gupxyL*k_xU_yD+gupxzL*k_xU_zD;
        double kupxyL = gupxyL*k_xU_xD+gupyyL*k_xU_yD+gupyzL*k_xU_zD;
        double kupxzL = gupxzL*k_xU_xD+gupyzL*k_xU_yD+gupzzL*k_xU_zD;
        double kupyyL = gupxyL*k_xD_yU+gupyyL*k_yU_yD+gupyzL*k_yU_zD;
        double kupyzL = gupxzL*k_xD_yU+gupyzL*k_yU_yD+gupzzL*k_yU_zD;
        double kupzzL = gupxzL*k_xD_zU+gupyzL*k_yD_zU+gupzzL*k_zU_zD;


	// Init Christoffel symbols from 3 metric
	double christoffel_xU_xxD{0.0};
        double christoffel_xU_xyD{0.0};
        double christoffel_xU_xzD{0.0};
        double christoffel_xU_yyD{0.0};
        double christoffel_xU_yzD{0.0};
        double christoffel_xU_zzD{0.0};

        double christoffel_yU_xxD{0.0};
        double christoffel_yU_xyD{0.0};
        double christoffel_yU_xzD{0.0};
	double christoffel_yU_yyD{0.0};
        double christoffel_yU_yzD{0.0};
        double christoffel_yU_zzD{0.0};

        double christoffel_zU_xxD{0.0};
        double christoffel_zU_xyD{0.0};
        double christoffel_zU_xzD{0.0};
        double christoffel_zU_yyD{0.0};
        double christoffel_zU_yzD{0.0};
        double christoffel_zU_zzD{0.0};

	double christoffel_xD_xxD{0.0};
        double christoffel_xD_xyD{0.0};
        double christoffel_xD_xzD{0.0};
        double christoffel_xD_yyD{0.0};
        double christoffel_xD_yzD{0.0};
        double christoffel_xD_zzD{0.0};

        double christoffel_yD_xxD{0.0};
        double christoffel_yD_xyD{0.0};
        double christoffel_yD_xzD{0.0};
	double christoffel_yD_yyD{0.0};
        double christoffel_yD_yzD{0.0};
        double christoffel_yD_zzD{0.0};

        double christoffel_zD_xxD{0.0};
        double christoffel_zD_xyD{0.0};
        double christoffel_zD_xzD{0.0};
        double christoffel_zD_yyD{0.0};
        double christoffel_zD_yzD{0.0};
        double christoffel_zD_zzD{0.0};

	// Init partial derivatives of lapse and shift	
	double partial_xD_lapse{0.};
        double partial_yD_lapse{0.};
        double partial_zD_lapse{0.};

	double partial_xD_shift_xU{0.};
        double partial_xD_shift_yU{0.};
        double partial_xD_shift_zU{0.};
        double partial_yD_shift_xU{0.};
        double partial_yD_shift_yU{0.};
        double partial_yD_shift_zU{0.};
        double partial_zD_shift_xU{0.};
        double partial_zD_shift_yU{0.};
        double partial_zD_shift_zU{0.};

	// Init covariant derivatives of W*v^i
	double cov_xD_Wv_xU{0.};
        double cov_xD_Wv_yU{0.};
        double cov_xD_Wv_zU{0.};
        double cov_yD_Wv_xU{0.};
        double cov_yD_Wv_yU{0.};
        double cov_yD_Wv_zU{0.};
        double cov_zD_Wv_xU{0.};
        double cov_zD_Wv_yU{0.};
        double cov_zD_Wv_zU{0.};

	double cov_xU_Wv_xU{0.};
        double cov_xU_Wv_yU{0.};
        double cov_xU_Wv_zU{0.};
        double cov_yU_Wv_xU{0.};
        double cov_yU_Wv_yU{0.};
        double cov_yU_Wv_zU{0.};
        double cov_zU_Wv_xU{0.};
        double cov_zU_Wv_yU{0.};
        double cov_zU_Wv_zU{0.};

	// Init partial derivatives of Lorentz factor	
	double partial_xD_W{0.};
	double partial_yD_W{0.};
	double partial_zD_W{0.};

	// Init partial derivatives of W*v^i
	double partial_xD_Wv_xU{0.};
        double partial_xD_Wv_yU{0.};
        double partial_xD_Wv_zU{0.};
        double partial_yD_Wv_xU{0.};
        double partial_yD_Wv_yU{0.};
        double partial_yD_Wv_zU{0.};
        double partial_zD_Wv_xU{0.};
        double partial_zD_Wv_yU{0.};
        double partial_zD_Wv_zU{0.};

	// Init partial derivatives of magnetic field
	double partial_xD_alphaB_yD{0.};
	double partial_xD_alphaB_zD{0.};
	double partial_yD_alphaB_xD{0.};
	double partial_yD_alphaB_zD{0.};
	double partial_zD_alphaB_xD{0.};
	double partial_zD_alphaB_yD{0.};


	// ######################################################
        // Compute derivatives with respect to x,y,z
	// ######################################################

        // #pragma unroll
        for(int dir=1;dir<4;dir++){

        	// int indexm2 = CCTK_GFINDEX3D(cctkGH,i-2*kronecker_delta[dir][1],j-2*kronecker_delta[dir][2],k-2*kronecker_delta[dir][3]);
        	// int indexm1 = CCTK_GFINDEX3D(cctkGH,i-  kronecker_delta[dir][1],j-  kronecker_delta[dir][2],k-  kronecker_delta[dir][3]);
        	// int indexp1 = CCTK_GFINDEX3D(cctkGH,i+  kronecker_delta[dir][1],j+  kronecker_delta[dir][2],k+  kronecker_delta[dir][3]);
        	// int indexp2 = CCTK_GFINDEX3D(cctkGH,i+2*kronecker_delta[dir][1],j+2*kronecker_delta[dir][2],k+2*kronecker_delta[dir][3]);

        	const auto indexm2 = p.I - 2*p.DI[dir];
        	const auto indexm1 = p.I -   p.DI[dir];
        	const auto indexp1 = p.I +   p.DI[dir];
        	const auto indexp2 = p.I + 2*p.DI[dir];

         	// Metric stencil
          const smat<CCTK_REAL, 3> gm2(
            [&](int i, int j) ARITH_INLINE { return calc_avg_v2c(gf_k(i, j), p, indexm2); });
        	double gxxLm2 = gm2(0, 0);
        	double gxyLm2 = gm2(0, 1); 
        	double gxzLm2 = gm2(0, 2); 
        	double gyyLm2 = gm2(1, 1); 
        	double gyzLm2 = gm2(1, 2); 
        	double gzzLm2 = gm2(2, 2); 

          const smat<CCTK_REAL, 3> gm1(
            [&](int i, int j) ARITH_INLINE { return calc_avg_v2c(gf_k(i, j), p, indexm1); });
        	double gxxLm1 = gm1(0, 0);
        	double gxyLm1 = gm1(0, 1); 
        	double gxzLm1 = gm1(0, 2); 
        	double gyyLm1 = gm1(1, 1); 
        	double gyzLm1 = gm1(1, 2); 
        	double gzzLm1 = gm1(2, 2); 

          const smat<CCTK_REAL, 3> gp1(
            [&](int i, int j) ARITH_INLINE { return calc_avg_v2c(gf_k(i, j), p, indexp1); });
        	double gxxLp1 = gp1(0, 0);
        	double gxyLp1 = gp1(0, 1); 
        	double gxzLp1 = gp1(0, 2); 
        	double gyyLp1 = gp1(1, 1); 
        	double gyzLp1 = gp1(1, 2); 
        	double gzzLp1 = gp1(2, 2); 

          const smat<CCTK_REAL, 3> gp2(
            [&](int i, int j) ARITH_INLINE { return calc_avg_v2c(gf_k(i, j), p, indexp2); });
        	double gxxLp2 = gp2(0, 0);
        	double gxyLp2 = gp2(0, 1); 
        	double gxzLp2 = gp2(0, 2); 
        	double gyyLp2 = gp2(1, 1); 
        	double gyzLp2 = gp2(1, 2); 
        	double gzzLp2 = gp2(2, 2); 

		// Lapse and shift stencil
          double lapsem2 = calc_avg_v2c(alp, p, indexm2);
          const vec<CCTK_REAL, 3> betasm2(
              [&](int i) ARITH_INLINE { return calc_avg_v2c(gf_beta(i), p, indexm2); });
        	double shiftxm2 = betasm2(0);
        	double shiftym2 = betasm2(1);
        	double shiftzm2 = betasm2(2);
	
          double lapsem1 = calc_avg_v2c(alp, p, indexm1);
          const vec<CCTK_REAL, 3> betasm1(
              [&](int i) ARITH_INLINE { return calc_avg_v2c(gf_beta(i), p, indexm1); });
        	double shiftxm1 = betasm1(0);
        	double shiftym1 = betasm1(1);
        	double shiftzm1 = betasm1(2);
        
          double lapsep1 = calc_avg_v2c(alp, p, indexp1);
          const vec<CCTK_REAL, 3> betasp1(
              [&](int i) ARITH_INLINE { return calc_avg_v2c(gf_beta(i), p, indexp1); });
        	double shiftxp1 = betasp1(0);
        	double shiftyp1 = betasp1(1);
        	double shiftzp1 = betasp1(2);

          double lapsep2 = calc_avg_v2c(alp, p, indexp2);
          const vec<CCTK_REAL, 3> betasp2(
              [&](int i) ARITH_INLINE { return calc_avg_v2c(gf_beta(i), p, indexp2); });
        	double shiftxp2 = betasp2(0);
        	double shiftyp2 = betasp2(1);
        	double shiftzp2 = betasp2(2);

                // Fluid 3 velocity in normal frame and Lorentz factor stencil
	        // double ETvxm2 = vel[CCTK_GFINDEX4D(cctkGH,i-2*kronecker_delta[dir][1],j-2*kronecker_delta[dir][2],k-2*kronecker_delta[dir][3],0)];
	        // double ETvym2 = vel[CCTK_GFINDEX4D(cctkGH,i-2*kronecker_delta[dir][1],j-2*kronecker_delta[dir][2],k-2*kronecker_delta[dir][3],1)];
	        // double ETvzm2 = vel[CCTK_GFINDEX4D(cctkGH,i-2*kronecker_delta[dir][1],j-2*kronecker_delta[dir][2],k-2*kronecker_delta[dir][3],2)];
	        double ETvxm2 = velx(indexm2);
	        double ETvym2 = vely(indexm2);
	        double ETvzm2 = velz(indexm2);
          double lfacm2 = w_lorentz(indexm2);

	        double ETvxm1 = velx(indexm1);
	        double ETvym1 = vely(indexm1);
	        double ETvzm1 = velz(indexm1);
          double lfacm1 = w_lorentz(indexm1);

	        double ETvxp1 = velx(indexp1);
	        double ETvyp1 = vely(indexp1);
	        double ETvzp1 = velz(indexp1);
          double lfacp1 = w_lorentz(indexp1);

	        double ETvxp2 = velx(indexp2);
	        double ETvyp2 = vely(indexp2);
	        double ETvzp2 = velz(indexp2);
          double lfacp2 = w_lorentz(indexp2);

		// Central magnetic field in normal frame stencil, CONTRAVARIANT vector
	        double Bx_m2 = Bvecx(indexm2);
	        double By_m2 = Bvecy(indexm2);
	        double Bz_m2 = Bvecz(indexm2);

	        double Bx_m1 = Bvecx(indexm1);
	        double By_m1 = Bvecy(indexm1);
	        double Bz_m1 = Bvecz(indexm1);

	        double Bx_p1 = Bvecx(indexp1);
	        double By_p1 = Bvecy(indexp1);
	        double Bz_p1 = Bvecz(indexp1);

	        double Bx_p2 = Bvecx(indexp2);
	        double By_p2 = Bvecy(indexp2);
	        double Bz_p2 = Bvecz(indexp2);

		//Central magnetic field in normal frame, covariant vector
                double Bcov_m2[3] = {gxxLm2*Bx_m2 + gxyLm2*By_m2 + gxzLm2*Bz_m2,
	                             gxyLm2*Bx_m2 + gyyLm2*By_m2 + gyzLm2*Bz_m2,
                                     gxzLm2*Bx_m2 + gyzLm2*By_m2 + gzzLm2*Bz_m2};

                double Bcov_m1[3] = {gxxLm1*Bx_m1 + gxyLm1*By_m1 + gxzLm1*Bz_m1,
	                             gxyLm1*Bx_m1 + gyyLm1*By_m1 + gyzLm1*Bz_m1,
                                     gxzLm1*Bx_m1 + gyzLm1*By_m1 + gzzLm1*Bz_m1};

                double Bcov_p1[3] = {gxxLp1*Bx_p1 + gxyLp1*By_p1 + gxzLp1*Bz_p1,
	                             gxyLp1*Bx_p1 + gyyLp1*By_p1 + gyzLp1*Bz_p1,
                                     gxzLp1*Bx_p1 + gyzLp1*By_p1 + gzzLp1*Bz_p1};

                double Bcov_p2[3] = {gxxLp2*Bx_p2 + gxyLp2*By_p2 + gxzLp2*Bz_p2,
	                             gxyLp2*Bx_p2 + gyyLp2*By_p2 + gyzLp2*Bz_p2,
                                     gxzLp2*Bx_p2 + gyzLp2*By_p2 + gzzLp2*Bz_p2};

		// #####################################################
                // Calculate Christoffel symbols having all indices down
		// First index x
          	christoffel_xD_xxD +=     0.5*COMPUTE_DERIV_4(gxxLm2,gxxLm1,gxxLp1,gxxLp2)*dxi[1]*kronecker_delta[dir][1]; 
                christoffel_xD_xyD +=     0.5*COMPUTE_DERIV_4(gxxLm2,gxxLm1,gxxLp1,gxxLp2)*dxi[2]*kronecker_delta[dir][2];
                christoffel_xD_xzD +=     0.5*COMPUTE_DERIV_4(gxxLm2,gxxLm1,gxxLp1,gxxLp2)*dxi[3]*kronecker_delta[dir][3];
                christoffel_xD_yyD += 0.5*(2.*COMPUTE_DERIV_4(gxyLm2,gxyLm1,gxyLp1,gxyLp2)*dxi[2]*kronecker_delta[dir][2]
				             -COMPUTE_DERIV_4(gyyLm2,gyyLm1,gyyLp1,gyyLp2)*dxi[1]*kronecker_delta[dir][1]);
                christoffel_xD_yzD +=    0.5*(COMPUTE_DERIV_4(gxyLm2,gxyLm1,gxyLp1,gxyLp2)*dxi[3]*kronecker_delta[dir][3]
			                     +COMPUTE_DERIV_4(gxzLm2,gxzLm1,gxzLp1,gxzLp2)*dxi[2]*kronecker_delta[dir][2]
				             -COMPUTE_DERIV_4(gyzLm2,gyzLm1,gyzLp1,gyzLp2)*dxi[1]*kronecker_delta[dir][1]);
                christoffel_xD_zzD += 0.5*(2.*COMPUTE_DERIV_4(gxzLm2,gxzLm1,gxzLp1,gxzLp2)*dxi[3]*kronecker_delta[dir][3]
				             -COMPUTE_DERIV_4(gzzLm2,gzzLm1,gzzLp1,gzzLp2)*dxi[1]*kronecker_delta[dir][1]);

		// First index y
          	christoffel_yD_xxD += 0.5*(2.*COMPUTE_DERIV_4(gxyLm2,gxyLm1,gxyLp1,gxyLp2)*dxi[1]*kronecker_delta[dir][1]
				             -COMPUTE_DERIV_4(gxxLm2,gxxLm1,gxxLp1,gxxLp2)*dxi[2]*kronecker_delta[dir][2]); 
                christoffel_yD_xyD +=     0.5*COMPUTE_DERIV_4(gyyLm2,gyyLm1,gyyLp1,gyyLp2)*dxi[1]*kronecker_delta[dir][1];
                christoffel_yD_xzD +=    0.5*(COMPUTE_DERIV_4(gxyLm2,gxyLm1,gxyLp1,gxyLp2)*dxi[3]*kronecker_delta[dir][3]
				             +COMPUTE_DERIV_4(gyzLm2,gyzLm1,gyzLp1,gyzLp2)*dxi[1]*kronecker_delta[dir][1]
				             -COMPUTE_DERIV_4(gxzLm2,gxzLm1,gxzLp1,gxzLp2)*dxi[2]*kronecker_delta[dir][2]);
                christoffel_yD_yyD +=     0.5*COMPUTE_DERIV_4(gyyLm2,gyyLm1,gyyLp1,gyyLp2)*dxi[2]*kronecker_delta[dir][2];
                christoffel_yD_yzD +=     0.5*COMPUTE_DERIV_4(gyyLm2,gyyLm1,gyyLp1,gyyLp2)*dxi[3]*kronecker_delta[dir][3];
                christoffel_yD_zzD += 0.5*(2.*COMPUTE_DERIV_4(gyzLm2,gyzLm1,gyzLp1,gyzLp2)*dxi[3]*kronecker_delta[dir][3]
				             -COMPUTE_DERIV_4(gzzLm2,gzzLm1,gzzLp1,gzzLp2)*dxi[2]*kronecker_delta[dir][2]);

                // First index z
          	christoffel_zD_xxD += 0.5*(2.*COMPUTE_DERIV_4(gxzLm2,gxzLm1,gxzLp1,gxzLp2)*dxi[1]*kronecker_delta[dir][1]
				             -COMPUTE_DERIV_4(gxxLm2,gxxLm1,gxxLp1,gxxLp2)*dxi[3]*kronecker_delta[dir][3]);  
                christoffel_zD_xyD +=    0.5*(COMPUTE_DERIV_4(gxzLm2,gxzLm1,gxzLp1,gxzLp2)*dxi[2]*kronecker_delta[dir][2]
				             +COMPUTE_DERIV_4(gyzLm2,gyzLm1,gyzLp1,gyzLp2)*dxi[1]*kronecker_delta[dir][1]
				             -COMPUTE_DERIV_4(gxyLm2,gxyLm1,gxyLp1,gxyLp2)*dxi[3]*kronecker_delta[dir][3]);
                christoffel_zD_xzD +=     0.5*COMPUTE_DERIV_4(gzzLm2,gzzLm1,gzzLp1,gzzLp2)*dxi[1]*kronecker_delta[dir][1];
                christoffel_zD_yyD += 0.5*(2.*COMPUTE_DERIV_4(gyzLm2,gyzLm1,gyzLp1,gyzLp2)*dxi[2]*kronecker_delta[dir][2]
				             -COMPUTE_DERIV_4(gyyLm2,gyyLm1,gyyLp1,gyyLp2)*dxi[3]*kronecker_delta[dir][3]);
                christoffel_zD_yzD +=     0.5*COMPUTE_DERIV_4(gzzLm2,gzzLm1,gzzLp1,gzzLp2)*dxi[2]*kronecker_delta[dir][2];
                christoffel_zD_zzD +=     0.5*COMPUTE_DERIV_4(gzzLm2,gzzLm1,gzzLp1,gzzLp2)*dxi[3]*kronecker_delta[dir][3]; 

		// ###############################################################
	        // Calculate partial derivatives of lapse and shift in x direction
		partial_xD_shift_xU += COMPUTE_DERIV_4(shiftxm2,shiftxm1,shiftxp1,shiftxp2)*dxi[1]*kronecker_delta[dir][1];
                partial_xD_shift_yU += COMPUTE_DERIV_4(shiftym2,shiftym1,shiftyp1,shiftyp2)*dxi[1]*kronecker_delta[dir][1];
                partial_xD_shift_zU += COMPUTE_DERIV_4(shiftzm2,shiftzm1,shiftzp1,shiftzp2)*dxi[1]*kronecker_delta[dir][1];
                partial_xD_lapse    += COMPUTE_DERIV_4(lapsem2,lapsem1,lapsep1,lapsep2)*dxi[1]*kronecker_delta[dir][1];

	        // Calculate partial derivatives of lapse and shift in y direction
		partial_yD_shift_xU += COMPUTE_DERIV_4(shiftxm2,shiftxm1,shiftxp1,shiftxp2)*dxi[2]*kronecker_delta[dir][2];
                partial_yD_shift_yU += COMPUTE_DERIV_4(shiftym2,shiftym1,shiftyp1,shiftyp2)*dxi[2]*kronecker_delta[dir][2];
                partial_yD_shift_zU += COMPUTE_DERIV_4(shiftzm2,shiftzm1,shiftzp1,shiftzp2)*dxi[2]*kronecker_delta[dir][2];
                partial_yD_lapse    += COMPUTE_DERIV_4(lapsem2,lapsem1,lapsep1,lapsep2)*dxi[2]*kronecker_delta[dir][2];

	        // Calculate partial derivatives of lapse and shift in z direction
		partial_zD_shift_xU += COMPUTE_DERIV_4(shiftxm2,shiftxm1,shiftxp1,shiftxp2)*dxi[3]*kronecker_delta[dir][3];
                partial_zD_shift_yU += COMPUTE_DERIV_4(shiftym2,shiftym1,shiftyp1,shiftyp2)*dxi[3]*kronecker_delta[dir][3];
                partial_zD_shift_zU += COMPUTE_DERIV_4(shiftzm2,shiftzm1,shiftzp1,shiftzp2)*dxi[3]*kronecker_delta[dir][3];
                partial_zD_lapse    += COMPUTE_DERIV_4(lapsem2,lapsem1,lapsep1,lapsep2)*dxi[3]*kronecker_delta[dir][3];

		// #####################################################
	        // Calculate partial derivatives of W*v^i in x direction
		partial_xD_Wv_xU += COMPUTE_DERIV_4(lfacm2*ETvxm2,lfacm1*ETvxm1,lfacp1*ETvxp1,lfacp2*ETvxp2)*dxi[1]*kronecker_delta[dir][1];
                partial_xD_Wv_yU += COMPUTE_DERIV_4(lfacm2*ETvym2,lfacm1*ETvym1,lfacp1*ETvyp1,lfacp2*ETvyp2)*dxi[1]*kronecker_delta[dir][1];
                partial_xD_Wv_zU += COMPUTE_DERIV_4(lfacm2*ETvzm2,lfacm1*ETvzm1,lfacp1*ETvzp1,lfacp2*ETvzp2)*dxi[1]*kronecker_delta[dir][1];

	        // Calculate partial derivatives of W*v^i in y direction
		partial_yD_Wv_xU += COMPUTE_DERIV_4(lfacm2*ETvxm2,lfacm1*ETvxm1,lfacp1*ETvxp1,lfacp2*ETvxp2)*dxi[2]*kronecker_delta[dir][2];
                partial_yD_Wv_yU += COMPUTE_DERIV_4(lfacm2*ETvym2,lfacm1*ETvym1,lfacp1*ETvyp1,lfacp2*ETvyp2)*dxi[2]*kronecker_delta[dir][2];
                partial_yD_Wv_zU += COMPUTE_DERIV_4(lfacm2*ETvzm2,lfacm1*ETvzm1,lfacp1*ETvzp1,lfacp2*ETvzp2)*dxi[2]*kronecker_delta[dir][2];

	        // Calculate partial derivatives of W*v^i in z direction
		partial_zD_Wv_xU += COMPUTE_DERIV_4(lfacm2*ETvxm2,lfacm1*ETvxm1,lfacp1*ETvxp1,lfacp2*ETvxp2)*dxi[3]*kronecker_delta[dir][3];
                partial_zD_Wv_yU += COMPUTE_DERIV_4(lfacm2*ETvym2,lfacm1*ETvym1,lfacp1*ETvyp1,lfacp2*ETvyp2)*dxi[3]*kronecker_delta[dir][3];
                partial_zD_Wv_zU += COMPUTE_DERIV_4(lfacm2*ETvzm2,lfacm1*ETvzm1,lfacp1*ETvzp1,lfacp2*ETvzp2)*dxi[3]*kronecker_delta[dir][3];


		// ####################################
		// Partial derivatives of Lorentz factor
		partial_xD_W += COMPUTE_DERIV_4(lfacm2,lfacm1,lfacp1,lfacp2)*dxi[1]*kronecker_delta[dir][1];
		partial_yD_W += COMPUTE_DERIV_4(lfacm2,lfacm1,lfacp1,lfacp2)*dxi[2]*kronecker_delta[dir][2];
		partial_zD_W += COMPUTE_DERIV_4(lfacm2,lfacm1,lfacp1,lfacp2)*dxi[3]*kronecker_delta[dir][3];

		// #####################################################
	        // Calculate partial derivatives of alpha*B_i in x direction
                partial_xD_alphaB_yD += COMPUTE_DERIV_4(lapsem2*Bcov_m2[1],lapsem1*Bcov_m1[1],lapsep1*Bcov_p1[1],lapsep2*Bcov_p2[1])*dxi[1]*kronecker_delta[dir][1];
                partial_xD_alphaB_zD += COMPUTE_DERIV_4(lapsem2*Bcov_m2[2],lapsem1*Bcov_m1[2],lapsep1*Bcov_p1[2],lapsep2*Bcov_p2[2])*dxi[1]*kronecker_delta[dir][1];

	        // Calculate partial derivatives of alpha*B_i in y direction
		partial_yD_alphaB_xD += COMPUTE_DERIV_4(lapsem2*Bcov_m2[0],lapsem1*Bcov_m1[0],lapsep1*Bcov_p1[0],lapsep2*Bcov_p2[0])*dxi[2]*kronecker_delta[dir][2];
                partial_yD_alphaB_zD += COMPUTE_DERIV_4(lapsem2*Bcov_m2[2],lapsem1*Bcov_m1[2],lapsep1*Bcov_p1[2],lapsep2*Bcov_p2[2])*dxi[2]*kronecker_delta[dir][2];

	        // Calculate partial derivatives of alpha*B_i in z direction
		partial_zD_alphaB_xD += COMPUTE_DERIV_4(lapsem2*Bcov_m2[0],lapsem1*Bcov_m1[0],lapsep1*Bcov_p1[0],lapsep2*Bcov_p2[0])*dxi[3]*kronecker_delta[dir][3];
                partial_zD_alphaB_yD += COMPUTE_DERIV_4(lapsem2*Bcov_m2[1],lapsem1*Bcov_m1[1],lapsep1*Bcov_p1[1],lapsep2*Bcov_p2[1])*dxi[3]*kronecker_delta[dir][3];


	}

	// ######################################################
	// Now: Calculate Christoffel symbols having first index up

	christoffel_xU_xxD = gupxxL*christoffel_xD_xxD+gupxyL*christoffel_yD_xxD+gupxzL*christoffel_zD_xxD; 
        christoffel_xU_xyD = gupxxL*christoffel_xD_xyD+gupxyL*christoffel_yD_xyD+gupxzL*christoffel_zD_xyD;
        christoffel_xU_xzD = gupxxL*christoffel_xD_xzD+gupxyL*christoffel_yD_xzD+gupxzL*christoffel_zD_xzD;
        christoffel_xU_yyD = gupxxL*christoffel_xD_yyD+gupxyL*christoffel_yD_yyD+gupxzL*christoffel_zD_yyD;
        christoffel_xU_yzD = gupxxL*christoffel_xD_yzD+gupxyL*christoffel_yD_yzD+gupxzL*christoffel_zD_yzD;
        christoffel_xU_zzD = gupxxL*christoffel_xD_zzD+gupxyL*christoffel_yD_zzD+gupxzL*christoffel_zD_zzD;

        christoffel_yU_xxD = gupxyL*christoffel_xD_xxD+gupyyL*christoffel_yD_xxD+gupyzL*christoffel_zD_xxD;
        christoffel_yU_xyD = gupxyL*christoffel_xD_xyD+gupyyL*christoffel_yD_xyD+gupyzL*christoffel_zD_xyD;
        christoffel_yU_xzD = gupxyL*christoffel_xD_xzD+gupyyL*christoffel_yD_xzD+gupyzL*christoffel_zD_xzD;
	christoffel_yU_yyD = gupxyL*christoffel_xD_yyD+gupyyL*christoffel_yD_yyD+gupyzL*christoffel_zD_yyD;
        christoffel_yU_yzD = gupxyL*christoffel_xD_yzD+gupyyL*christoffel_yD_yzD+gupyzL*christoffel_zD_yzD;
        christoffel_yU_zzD = gupxyL*christoffel_xD_zzD+gupyyL*christoffel_yD_zzD+gupyzL*christoffel_zD_zzD;

        christoffel_zU_xxD = gupxzL*christoffel_xD_xxD+gupyzL*christoffel_yD_xxD+gupzzL*christoffel_zD_xxD;
        christoffel_zU_xyD = gupxzL*christoffel_xD_xyD+gupyzL*christoffel_yD_xyD+gupzzL*christoffel_zD_xyD;
        christoffel_zU_xzD = gupxzL*christoffel_xD_xzD+gupyzL*christoffel_yD_xzD+gupzzL*christoffel_zD_xzD;
        christoffel_zU_yyD = gupxzL*christoffel_xD_yyD+gupyzL*christoffel_yD_yyD+gupzzL*christoffel_zD_yyD;
        christoffel_zU_yzD = gupxzL*christoffel_xD_yzD+gupyzL*christoffel_yD_yzD+gupzzL*christoffel_zD_yzD;
        christoffel_zU_zzD = gupxzL*christoffel_xD_zzD+gupyzL*christoffel_yD_zzD+gupzzL*christoffel_zD_zzD;

	// ####################################################
        // Now: Calculate spatial covariant derivatives of Wv^i

        // First index down	
        cov_xD_Wv_xU = partial_xD_Wv_xU + christoffel_xU_xxD*lfac*ETvx + christoffel_xU_xyD*lfac*ETvy + christoffel_xU_xzD*lfac*ETvz;  
        cov_xD_Wv_yU = partial_xD_Wv_yU + christoffel_yU_xxD*lfac*ETvx + christoffel_yU_xyD*lfac*ETvy + christoffel_yU_xzD*lfac*ETvz;
        cov_xD_Wv_zU = partial_xD_Wv_zU + christoffel_zU_xxD*lfac*ETvx + christoffel_zU_xyD*lfac*ETvy + christoffel_zU_xzD*lfac*ETvz;

	cov_yD_Wv_xU = partial_yD_Wv_xU + christoffel_xU_xyD*lfac*ETvx + christoffel_xU_yyD*lfac*ETvy + christoffel_xU_yzD*lfac*ETvz;  
        cov_yD_Wv_yU = partial_yD_Wv_yU + christoffel_yU_xyD*lfac*ETvx + christoffel_yU_yyD*lfac*ETvy + christoffel_yU_yzD*lfac*ETvz;
        cov_yD_Wv_zU = partial_yD_Wv_zU + christoffel_zU_xyD*lfac*ETvx + christoffel_zU_yyD*lfac*ETvy + christoffel_zU_yzD*lfac*ETvz;

        cov_zD_Wv_xU = partial_zD_Wv_xU + christoffel_xU_xzD*lfac*ETvx + christoffel_xU_yzD*lfac*ETvy + christoffel_xU_zzD*lfac*ETvz;  
        cov_zD_Wv_yU = partial_zD_Wv_yU + christoffel_yU_xzD*lfac*ETvx + christoffel_yU_yzD*lfac*ETvy + christoffel_yU_zzD*lfac*ETvz;
        cov_zD_Wv_zU = partial_zD_Wv_zU + christoffel_zU_xzD*lfac*ETvx + christoffel_zU_yzD*lfac*ETvy + christoffel_zU_zzD*lfac*ETvz;
	
	// First index up
        cov_xU_Wv_xU = gupxxL*cov_xD_Wv_xU+gupxyL*cov_yD_Wv_xU+gupxzL*cov_zD_Wv_xU; 
        cov_xU_Wv_yU = gupxxL*cov_xD_Wv_yU+gupxyL*cov_yD_Wv_yU+gupxzL*cov_zD_Wv_yU;
        cov_xU_Wv_zU = gupxxL*cov_xD_Wv_zU+gupxyL*cov_yD_Wv_zU+gupxzL*cov_zD_Wv_zU;

	cov_yU_Wv_xU = gupxyL*cov_xD_Wv_xU+gupyyL*cov_yD_Wv_xU+gupyzL*cov_zD_Wv_xU;
        cov_yU_Wv_yU = gupxyL*cov_xD_Wv_yU+gupyyL*cov_yD_Wv_yU+gupyzL*cov_zD_Wv_yU;
        cov_yU_Wv_zU = gupxyL*cov_xD_Wv_zU+gupyyL*cov_yD_Wv_zU+gupyzL*cov_zD_Wv_zU;

        cov_zU_Wv_xU = gupxzL*cov_xD_Wv_xU+gupyzL*cov_yD_Wv_xU+gupzzL*cov_zD_Wv_xU; 
        cov_zU_Wv_yU = gupxzL*cov_xD_Wv_yU+gupyzL*cov_yD_Wv_yU+gupzzL*cov_zD_Wv_yU;
        cov_zU_Wv_zU = gupxzL*cov_xD_Wv_zU+gupyzL*cov_yD_Wv_zU+gupzzL*cov_zD_Wv_zU;


	// ######################################################
	// Now: Calculate all intermediate tensors from 3+1 paper

	double expansion_3             = cov_xD_Wv_xU+cov_yD_Wv_yU+cov_zD_Wv_zU;

	double lambda_no_time_der      = -(shiftx*partial_xD_W+shifty*partial_yD_W+shiftz*partial_zD_W)/lapse;
	lambda_no_time_der += lfac*(ETvx*partial_xD_lapse+ETvy*partial_yD_lapse+ETvz*partial_zD_lapse)/lapse;
        // Add time derivative
	double lambda = lambda_no_time_der + partial_t_W/lapse;

        double k_trace = gupxxL*kxxL+gupyyL*kyyL+gupzzL*kzzL
	                +2.*gupxyL*kxyL+2.*gupxzL*kxzL+2.*gupyzL*kyzL;

	double lambda_vec_xU           = -(shiftx*partial_xD_Wv_xU+shifty*partial_yD_Wv_xU+shiftz*partial_zD_Wv_xU
			                   -lfac*ETvx*partial_xD_shift_xU-lfac*ETvy*partial_yD_shift_xU-lfac*ETvz*partial_zD_shift_xU)/lapse
		                           +lfac*(gupxxL*partial_xD_lapse+gupxyL*partial_yD_lapse+gupxzL*partial_zD_lapse)/lapse;
        // Add time derivative
	lambda_vec_xU += partial_t_Wv_xU/lapse;

	double lambda_vec_yU           = -(shiftx*partial_xD_Wv_yU+shifty*partial_yD_Wv_yU+shiftz*partial_zD_Wv_yU
			                   -lfac*ETvx*partial_xD_shift_yU-lfac*ETvy*partial_yD_shift_yU-lfac*ETvz*partial_zD_shift_yU)/lapse
		                           +lfac*(gupxyL*partial_xD_lapse+gupyyL*partial_yD_lapse+gupyzL*partial_zD_lapse)/lapse;
        // Add time derivative
	lambda_vec_yU += partial_t_Wv_yU/lapse;

	double lambda_vec_zU           = -(shiftx*partial_xD_Wv_zU+shifty*partial_yD_Wv_zU+shiftz*partial_zD_Wv_zU
			                   -lfac*ETvx*partial_xD_shift_zU-lfac*ETvy*partial_yD_shift_zU-lfac*ETvz*partial_zD_shift_zU)/lapse
		                           +lfac*(gupxzL*partial_xD_lapse+gupyzL*partial_yD_lapse+gupzzL*partial_zD_lapse)/lapse;
        // Add time derivative
	lambda_vec_zU += partial_t_Wv_zU/lapse;

	double acc_3_vec_xU            = lfac*(ETvx*cov_xD_Wv_xU+ETvy*cov_yD_Wv_xU+ETvz*cov_zD_Wv_xU);
	double acc_3_vec_yU            = lfac*(ETvx*cov_xD_Wv_yU+ETvy*cov_yD_Wv_yU+ETvz*cov_zD_Wv_yU);
	double acc_3_vec_zU            = lfac*(ETvx*cov_xD_Wv_zU+ETvy*cov_yD_Wv_zU+ETvz*cov_zD_Wv_zU);

	double shear_3_xxU             = cov_xU_Wv_xU + lfac*acc_3_vec_xU*ETvx - expansion_3*(gupxxL+lfac*lfac*ETvx*ETvx)/3.;

	double shear_3_xyU             = (cov_xU_Wv_yU + cov_yU_Wv_xU)/2. + lfac*(acc_3_vec_xU*ETvy + acc_3_vec_yU*ETvx)/2.
	                           	- expansion_3*(gupxyL+lfac*lfac*ETvx*ETvy)/3.;
	
	double shear_3_xzU             = (cov_xU_Wv_zU + cov_zU_Wv_xU)/2. + lfac*(acc_3_vec_xU*ETvz + acc_3_vec_zU*ETvx)/2.
	                           	- expansion_3*(gupxzL+lfac*lfac*ETvx*ETvz)/3.;

	double shear_3_yyU             = cov_yU_Wv_yU + lfac*acc_3_vec_yU*ETvy - expansion_3*(gupyyL+lfac*lfac*ETvy*ETvy)/3.;

	double shear_3_yzU             = (cov_yU_Wv_zU + cov_zU_Wv_yU)/2. + lfac*(acc_3_vec_yU*ETvz + acc_3_vec_zU*ETvy)/2.
	                           	- expansion_3*(gupyzL+lfac*lfac*ETvy*ETvz)/3.;

	double shear_3_zzU             = cov_zU_Wv_zU + lfac*acc_3_vec_zU*ETvz - expansion_3*(gupzzL+lfac*lfac*ETvz*ETvz)/3.;

	double omega_3_xyU             = (cov_xU_Wv_yU - cov_yU_Wv_xU)/2. - lfac*(acc_3_vec_xU*ETvy - acc_3_vec_yU*ETvx)/2.;

	double omega_3_xzU             = (cov_xU_Wv_zU - cov_zU_Wv_xU)/2. - lfac*(acc_3_vec_xU*ETvz - acc_3_vec_zU*ETvx)/2.;

	double omega_3_yzU             = (cov_yU_Wv_zU - cov_zU_Wv_yU)/2. - lfac*(acc_3_vec_yU*ETvz - acc_3_vec_zU*ETvy)/2.;

	double lambda_tensor_xxU       = lfac*lfac*lambda_vec_xU*ETvx - lambda*(gupxxL+lfac*lfac*ETvx*ETvx)/3.;

	double lambda_tensor_xyU       = lfac*lfac*(lambda_vec_xU*ETvy + lambda_vec_yU*ETvx)/2. - lambda*(gupxyL+lfac*lfac*ETvx*ETvy)/3.;

	double lambda_tensor_xzU       = lfac*lfac*(lambda_vec_xU*ETvz + lambda_vec_zU*ETvx)/2. - lambda*(gupxzL+lfac*lfac*ETvx*ETvz)/3.;

	double lambda_tensor_yyU       = lfac*lfac*lambda_vec_yU*ETvy - lambda*(gupyyL+lfac*lfac*ETvy*ETvy)/3.;

	double lambda_tensor_yzU       = lfac*lfac*(lambda_vec_yU*ETvz + lambda_vec_zU*ETvy)/2. - lambda*(gupyzL+lfac*lfac*ETvy*ETvz)/3.;

	double lambda_tensor_zzU       = lfac*lfac*lambda_vec_zU*ETvz - lambda*(gupzzL+lfac*lfac*ETvz*ETvz)/3.;

	double k_dot_v_xU              = k_xU_xD*ETvx+k_xU_yD*ETvy+k_xU_zD*ETvz;
       	double k_dot_v_yU              = k_xD_yU*ETvx+k_yU_yD*ETvy+k_yU_zD*ETvz;
	double k_dot_v_zU              = k_xD_zU*ETvx+k_yD_zU*ETvy+k_zU_zD*ETvz;

	double k_tensor_xxU            = kupxxL + 2.*lfac*lfac*k_dot_v_xU*ETvx - k_trace*(gupxxL+lfac*lfac*ETvx*ETvx)/3.;
	double k_tensor_xyU            = kupxyL + lfac*lfac*(k_dot_v_xU*ETvy+k_dot_v_yU*ETvx) - k_trace*(gupxyL+lfac*lfac*ETvx*ETvy)/3.;
	double k_tensor_xzU            = kupxzL + lfac*lfac*(k_dot_v_xU*ETvz+k_dot_v_zU*ETvx) - k_trace*(gupxzL+lfac*lfac*ETvx*ETvz)/3.;
	double k_tensor_yyU            = kupyyL + 2.*lfac*lfac*k_dot_v_yU*ETvy-k_trace*(gupyyL+lfac*lfac*ETvy*ETvy)/3.;
	double k_tensor_yzU            = kupyzL + lfac*lfac*(k_dot_v_yU*ETvz+k_dot_v_zU*ETvy) - k_trace*(gupyzL+lfac*lfac*ETvy*ETvz)/3.;
	double k_tensor_zzU            = kupzzL + 2.*lfac*lfac*k_dot_v_zU*ETvz-k_trace*(gupzzL+lfac*lfac*ETvz*ETvz)/3.;

	// ################################
        // Finally, let's calculate the GFs

	expansion_scalar(p.I)               = expansion_3+lambda-k_trace*lfac;
	expansion_no_time_der(p.I)          = expansion_3+lambda_no_time_der-k_trace*lfac; 

	double shear_spatial_tensor_xx_tmp        = shear_3_xxU+lambda_tensor_xxU-lfac*k_tensor_xxU;
	double shear_spatial_tensor_xy_tmp        = shear_3_xyU+lambda_tensor_xyU-lfac*k_tensor_xyU;
	double shear_spatial_tensor_xz_tmp        = shear_3_xzU+lambda_tensor_xzU-lfac*k_tensor_xzU;
	double shear_spatial_tensor_yy_tmp        = shear_3_yyU+lambda_tensor_yyU-lfac*k_tensor_yyU;
	double shear_spatial_tensor_yz_tmp        = shear_3_yzU+lambda_tensor_yzU-lfac*k_tensor_yzU;
	double shear_spatial_tensor_zz_tmp        = shear_3_zzU+lambda_tensor_zzU-lfac*k_tensor_zzU;

	shear_spatial_tensor_xx(p.I)        = shear_spatial_tensor_xx_tmp;
	shear_spatial_tensor_xy(p.I)        = shear_spatial_tensor_xy_tmp;
	shear_spatial_tensor_xz(p.I)        = shear_spatial_tensor_xz_tmp;
	shear_spatial_tensor_yy(p.I)        = shear_spatial_tensor_yy_tmp;
	shear_spatial_tensor_yz(p.I)        = shear_spatial_tensor_yz_tmp;
	shear_spatial_tensor_zz(p.I)        = shear_spatial_tensor_zz_tmp;

	double ET_v_xD                        = gxxL*ETvx+gxyL*ETvy+gxzL*ETvz;
	double ET_v_yD                        = gxyL*ETvx+gyyL*ETvy+gyzL*ETvz;
	double ET_v_zD                        = gxzL*ETvx+gyzL*ETvy+gzzL*ETvz;

        double shear_spatial_vector_x       = ET_v_xD*shear_spatial_tensor_xx_tmp+ET_v_yD*shear_spatial_tensor_xy_tmp+ET_v_zD*shear_spatial_tensor_xz_tmp;
	double shear_spatial_vector_y       = ET_v_xD*shear_spatial_tensor_xy_tmp+ET_v_yD*shear_spatial_tensor_yy_tmp+ET_v_zD*shear_spatial_tensor_yz_tmp;
	double shear_spatial_vector_z       = ET_v_xD*shear_spatial_tensor_xz_tmp+ET_v_yD*shear_spatial_tensor_yz_tmp+ET_v_zD*shear_spatial_tensor_zz_tmp;

        double shear_spatial_trace                     = ET_v_xD*shear_spatial_vector_x+ET_v_yD*shear_spatial_vector_y+ET_v_zD*shear_spatial_vector_z;

	double kin_vorticity_spatial_xy_tmp       = omega_3_xyU - lfac*lfac*(lambda_vec_xU*ETvy-lambda_vec_yU*ETvx)/2. 
		                              + lfac*lfac*lfac*(k_dot_v_xU*ETvy-k_dot_v_yU*ETvx);

	double kin_vorticity_spatial_xz_tmp       = omega_3_xzU - lfac*lfac*(lambda_vec_xU*ETvz-lambda_vec_zU*ETvx)/2. 
		                              + lfac*lfac*lfac*(k_dot_v_xU*ETvz-k_dot_v_zU*ETvx);

	double kin_vorticity_spatial_yz_tmp       = omega_3_yzU - lfac*lfac*(lambda_vec_yU*ETvz-lambda_vec_zU*ETvy)/2. 
		                              + lfac*lfac*lfac*(k_dot_v_yU*ETvz-k_dot_v_zU*ETvy);

	double kin_acceleration_spatial_x_tmp     = acc_3_vec_xU + lfac*lambda_vec_xU - 2.*lfac*lfac*k_dot_v_xU;
	double kin_acceleration_spatial_y_tmp     = acc_3_vec_yU + lfac*lambda_vec_yU - 2.*lfac*lfac*k_dot_v_yU;
        double kin_acceleration_spatial_z_tmp     = acc_3_vec_zU + lfac*lambda_vec_zU - 2.*lfac*lfac*k_dot_v_zU;

        kin_acceleration_spatial_x(p.I) = kin_acceleration_spatial_x_tmp;
        kin_acceleration_spatial_y(p.I) = kin_acceleration_spatial_y_tmp;
        kin_acceleration_spatial_z(p.I) = kin_acceleration_spatial_z_tmp;

        kin_vorticity_spatial_xy(p.I) = kin_vorticity_spatial_xy_tmp;
	kin_vorticity_spatial_xz(p.I) = kin_vorticity_spatial_xz_tmp;
	kin_vorticity_spatial_yz(p.I) = kin_vorticity_spatial_yz_tmp;

        double spatial_vort_vec_x = ET_v_yD*kin_vorticity_spatial_xy_tmp+ET_v_zD*kin_vorticity_spatial_xz_tmp;
        double spatial_vort_vec_y = -ET_v_xD*kin_vorticity_spatial_xy_tmp+ET_v_zD*kin_vorticity_spatial_yz_tmp;
        double spatial_vort_vec_z = -ET_v_xD*kin_vorticity_spatial_xz_tmp-ET_v_yD*kin_vorticity_spatial_yz_tmp;

	// Calculate constraints_derivative_tensors and diagnostics_mhd groups of GFs 

	double shear4_ttU=shear_spatial_trace/lapse/lapse;
        double shear4_txU=-shear_spatial_trace*shiftx/lapse/lapse+shear_spatial_vector_x/lapse;
        double shear4_tyU=-shear_spatial_trace*shifty/lapse/lapse+shear_spatial_vector_y/lapse;
        double shear4_tzU=-shear_spatial_trace*shiftz/lapse/lapse+shear_spatial_vector_z/lapse;

        double shear4_xxU=shear_spatial_trace*shiftx*shiftx/lapse/lapse-shear_spatial_vector_x*shiftx/lapse-shear_spatial_vector_x*shiftx/lapse+shear_spatial_tensor_xx_tmp;

	double shear4_xyU=shear_spatial_trace*shiftx*shifty/lapse/lapse-shear_spatial_vector_x*shifty/lapse-shear_spatial_vector_y*shiftx/lapse+shear_spatial_tensor_xy_tmp;

	double shear4_xzU=shear_spatial_trace*shiftx*shiftz/lapse/lapse-shear_spatial_vector_x*shiftz/lapse-shear_spatial_vector_z*shiftx/lapse+shear_spatial_tensor_xz_tmp;

	double shear4_yyU=shear_spatial_trace*shifty*shifty/lapse/lapse-shear_spatial_vector_y*shifty/lapse-shear_spatial_vector_y*shifty/lapse+shear_spatial_tensor_yy_tmp;

	double shear4_yzU=shear_spatial_trace*shifty*shiftz/lapse/lapse-shear_spatial_vector_y*shiftz/lapse-shear_spatial_vector_z*shifty/lapse+shear_spatial_tensor_yz_tmp;

	double shear4_zzU=shear_spatial_trace*shiftz*shiftz/lapse/lapse-shear_spatial_vector_z*shiftz/lapse-shear_spatial_vector_z*shiftz/lapse+shear_spatial_tensor_zz_tmp;

                         
	double omega4_txU = -spatial_vort_vec_x/lapse;

        double omega4_tyU = -spatial_vort_vec_y/lapse;

        double omega4_tzU = -spatial_vort_vec_z/lapse;

        double omega4_xyU = -spatial_vort_vec_x*shifty/lapse + spatial_vort_vec_y*shiftx/lapse + kin_vorticity_spatial_xy_tmp; 	
        double omega4_xzU = -spatial_vort_vec_x*shiftz/lapse + spatial_vort_vec_z*shiftx/lapse + kin_vorticity_spatial_xz_tmp;

        double omega4_yzU = -spatial_vort_vec_y*shiftz/lapse + spatial_vort_vec_z*shifty/lapse + kin_vorticity_spatial_yz_tmp;


        double trace4_sigma = g4dn[0][0]*shear4_ttU+g4dn[0][1]*shear4_txU+g4dn[0][2]*shear4_tyU+g4dn[0][3]*shear4_tzU
                             +g4dn[1][0]*shear4_txU+g4dn[1][1]*shear4_xxU+g4dn[1][2]*shear4_xyU+g4dn[1][3]*shear4_xzU
			     +g4dn[2][0]*shear4_tyU+g4dn[2][1]*shear4_xyU+g4dn[2][2]*shear4_yyU+g4dn[2][3]*shear4_yzU
                             +g4dn[3][0]*shear4_tzU+g4dn[3][1]*shear4_xzU+g4dn[3][2]*shear4_yzU+g4dn[3][3]*shear4_zzU;

	double udown[4]{0,0,0,0};
	double smallbdown[4]{0,0,0,0};

        for(int crd=0;crd<4;crd++){
		udown[crd] = g4dn[crd][0]*uup[0]+g4dn[crd][1]*uup[1]+g4dn[crd][2]*uup[2]+g4dn[crd][3]*uup[3]; 
		smallbdown[crd] = g4dn[crd][0]*smallb[0]+g4dn[crd][1]*smallb[1]+g4dn[crd][2]*smallb[2]+g4dn[crd][3]*smallb[3];
	}

        double sigmadotBdotB = smallbdown[0]*smallbdown[0]*shear4_ttU
		              +smallbdown[0]*smallbdown[1]*shear4_txU
			      +smallbdown[0]*smallbdown[2]*shear4_tyU
			      +smallbdown[0]*smallbdown[3]*shear4_tzU  
                              +smallbdown[1]*smallbdown[0]*shear4_txU
			      +smallbdown[1]*smallbdown[1]*shear4_xxU
			      +smallbdown[1]*smallbdown[2]*shear4_xyU
			      +smallbdown[1]*smallbdown[3]*shear4_xzU
                              +smallbdown[2]*smallbdown[0]*shear4_tyU
			      +smallbdown[2]*smallbdown[1]*shear4_xyU
			      +smallbdown[2]*smallbdown[2]*shear4_yyU
			      +smallbdown[2]*smallbdown[3]*shear4_yzU
                              +smallbdown[3]*smallbdown[0]*shear4_tzU
			      +smallbdown[3]*smallbdown[1]*shear4_xzU
			      +smallbdown[3]*smallbdown[2]*shear4_yzU
			      +smallbdown[3]*smallbdown[3]*shear4_zzU;


	double sigmadotu_tU = udown[0]*shear4_ttU+udown[1]*shear4_txU+udown[2]*shear4_tyU+udown[3]*shear4_tzU;
        double sigmadotu_xU = udown[0]*shear4_txU+udown[1]*shear4_xxU+udown[2]*shear4_xyU+udown[3]*shear4_xzU;
        double sigmadotu_yU = udown[0]*shear4_tyU+udown[1]*shear4_xyU+udown[2]*shear4_yyU+udown[3]*shear4_yzU;
        double sigmadotu_zU = udown[0]*shear4_tzU+udown[1]*shear4_xzU+udown[2]*shear4_yzU+udown[3]*shear4_zzU;

	double omegadotu_tU =  udown[1]*omega4_txU+udown[2]*omega4_tyU+udown[3]*omega4_tzU;
        double omegadotu_xU = -udown[0]*omega4_txU+udown[2]*omega4_xyU+udown[3]*omega4_xzU;
        double omegadotu_yU = -udown[0]*omega4_tyU-udown[1]*omega4_xyU+udown[3]*omega4_yzU;
        double omegadotu_zU = -udown[0]*omega4_tzU-udown[1]*omega4_xzU-udown[2]*omega4_yzU;

        sigma4bb(p.I) = sigmadotBdotB;
	sigma4Ut(p.I) = sigmadotu_tU;
	sigma4Ux(p.I) = sigmadotu_xU;
	sigma4Uy(p.I) = sigmadotu_yU;
	sigma4Uz(p.I) = sigmadotu_zU;
	sigma4Trace(p.I) = trace4_sigma;
	omega4Ut(p.I) = omegadotu_tU;
	omega4Ux(p.I) = omegadotu_xU;
	omega4Uy(p.I) = omegadotu_yU;
	omega4Uz(p.I) = omegadotu_zU;

	// Calculate normcurlB
	
	double J_current[3] = { ONE_OVER_ALP_SQRTGAMMA_4PI*(partial_yD_alphaB_zD - partial_zD_alphaB_yD),
	                        ONE_OVER_ALP_SQRTGAMMA_4PI*(partial_zD_alphaB_xD - partial_xD_alphaB_zD),
	                        ONE_OVER_ALP_SQRTGAMMA_4PI*(partial_xD_alphaB_yD - partial_yD_alphaB_xD)};

        double Jcov[3] = {gxxL*J_current[0] + gxyL*J_current[1] + gxzL*J_current[2],
	                  gxyL*J_current[0] + gyyL*J_current[1] + gyzL*J_current[2],
                          gxzL*J_current[0] + gyzL*J_current[1] + gzzL*J_current[2]};

	double J2 = Jcov[0]*J_current[0] + Jcov[1]*J_current[1] + Jcov[2]*J_current[2];

	normcurlB(p.I) = sqrt(J2);

	// Calculate square of 4 acceleration
	
	double a_zero = ET_v_xD*kin_acceleration_spatial_x_tmp 
		      + ET_v_yD*kin_acceleration_spatial_y_tmp 
		      + ET_v_zD*kin_acceleration_spatial_z_tmp;

	double a4[4] = { a_zero/lapse,
                        -a_zero*shiftx/lapse+kin_acceleration_spatial_x_tmp,
	                -a_zero*shifty/lapse+kin_acceleration_spatial_y_tmp,
	                -a_zero*shiftz/lapse+kin_acceleration_spatial_z_tmp};

	double a4sqL=0.0; for(int ii=0;ii<4;ii++) for(int jj=0;jj<4;jj++) a4sqL += g4dn[ii][jj]*a4[ii]*a4[jj];

	a4sq(p.I) = a4sqL;
	} 

  // Cycle velocity timelevels
  velxold(p.I) = ETvx;
  velyold(p.I) = ETvy;
  velzold(p.I) = ETvz;
  w_lorentzold(p.I) = lfac;
    });
}

} // namespace
