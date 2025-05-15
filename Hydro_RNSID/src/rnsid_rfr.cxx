/* Changed by N. Stergioulas, 24/10/2001: include latest source files from cactus 3 */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "SpaceMask.h"

#include "cctk_DefineThorn.h" 
#include "rnsid.h"
#include <loop.hxx>
#include <loop_device.hxx>
#include <rnsreader_utils.hxx>


#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velx_p (&vel_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely_p (&vel_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz_p (&vel_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velx_p_p (&vel_p_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely_p_p (&vel_p_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz_p_p (&vel_p_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

#define sx (&scon[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sy (&scon[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sz (&scon[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sx_p (&scon_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sy_p (&scon_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sz_p (&scon_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sx_p_p (&scon_p_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sy_p_p (&scon_p_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sz_p_p (&scon_p_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

namespace RNSReader {
using namespace Loop;

static const char *rcsid="$Header$";
CCTK_FILEVERSION(Hydro_RNSID_rnsid_rfr_c)

//void Hydro_rnsid_init(CCTK_ARGUMENTS);

extern "C" void Hydro_rnsid_init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_Hydro_rnsid_init;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_REAL  Omega,                    /* angular velocity */
    R_e,                               /* circumferential radius */
    r_e,                               /* coordinate radius */
    mass0;                             /* rest mass */
    
  CCTK_REAL eos_k, eos_ideal_fluid_gamma, rnsid_rho_min;

  CCTK_REAL *x_coord=0, *y_coord=0, *z_coord=0;

  CCTK_REAL gamma_center;  /* central value of 3-determinant of metric,
                              needed for applying perturbations */

  if ( (RNS_K < 0.0) || (RNS_Gamma < 0.0)) {
    CCTK_WARN(0,"RNS_K and RNS_Gamma must be greater than 0: using 100.0 and 2!");
    eos_k                 = 100.0;
    eos_ideal_fluid_gamma = 2.0;
  } else { 
    eos_k                 = RNS_K;
    eos_ideal_fluid_gamma = RNS_Gamma;
  }

  rnsid_rho_min = RNS_rho_min;
 

 /* SET-UP INITIAL DATA */
  
  x_coord = x;
  y_coord = y;
  z_coord = z;
  
  Hydro_rnsid(cctkGH, x_coord, y_coord, z_coord,eos_k, eos_ideal_fluid_gamma, rnsid_rho_min, 
              &Omega, &R_e, &r_e, &mass0, &gamma_center);  
  
	/* Filling time levels. Not used.
 
  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::metric") > 1)
    {
      #pragma omp parallel for
      for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	{
	  gxx_p[i] = gxx[i];
	  gyy_p[i] = gyy[i];
	  gzz_p[i] = gzz[i];
	  gxy_p[i] = gxy[i];
	  gxz_p[i] = gxz[i];
	  gyz_p[i] = gyz[i];
	}
      if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::metric") > 2)
	{
	  for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	    {
	      gxx_p_p[i] = gxx[i];
	      gyy_p_p[i] = gyy[i];
	      gzz_p_p[i] = gzz[i];
	      gxy_p_p[i] = gxy[i];
	      gxz_p_p[i] = gxz[i];
	      gyz_p_p[i] = gyz[i];
	    }
	}
    }
  
  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::curv") > 1)
    {
      #pragma omp parallel for
      for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	{
	  kxx_p[i] = kxx[i];
	  kyy_p[i] = kyy[i];
	  kzz_p[i] = kzz[i];
	  kxy_p[i] = kxy[i];
	  kxz_p[i] = kxz[i];
	  kyz_p[i] = kyz[i];
	}
      if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::curv") > 2)
	{
          #pragma omp parallel for
	  for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	    {
	      kxx_p_p[i] = kxx[i];
	      kyy_p_p[i] = kyy[i];
	      kzz_p_p[i] = kzz[i];
	      kxy_p_p[i] = kxy[i];
	      kxz_p_p[i] = kxz[i];
	      kyz_p_p[i] = kyz[i];
	    }
	}
    }
  
  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::lapse") > 1)
    {
      #pragma omp parallel for
      for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	{
	  alp_p[i] = alp[i];
	}
      if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::lapse") > 2)
	{
          #pragma omp parallel for
	  for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	    {
	      alp_p_p[i] = alp[i];
	    }
	}
    }
  
  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::shift") > 1)
    {
      #pragma omp parallel for
      for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	{
	  betax_p[i] = betax[i];
	  betay_p[i] = betay[i];
	  betaz_p[i] = betaz[i];
	}
      if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::shift") > 2)
	{
          #pragma omp parallel for
	  for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	    {
	      betax_p_p[i] = betax[i];
	      betay_p_p[i] = betay[i];
	      betaz_p_p[i] = betaz[i];
	    }
	}
    }
  
  if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::rho") > 1)
    {
      #pragma omp parallel for
      for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	{
	  rho_p[i] = rho[i];
	  velx_p[i] = velx[i];
	  vely_p[i] = vely[i];
	  velz_p[i] = velz[i];
	  press_p[i] = press[i];
	  eps_p[i] = eps[i];
	  w_lorentz_p[i] = w_lorentz[i];
	}
      if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::rho") > 2)
	{
          #pragma omp parallel for
	  for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
	    {
	      rho_p_p[i] = rho[i];
	      velx_p_p[i] = velx[i];
	      vely_p_p[i] = vely[i];
	      velz_p_p[i] = velz[i];
	      press_p_p[i] = press[i];
	      eps_p_p[i] = eps[i];
	      w_lorentz_p_p[i] = w_lorentz[i]; 
	    }
	}
    }
	*/
  
  return;
}

/*
extern "C" void Hydro_RNSID_CheckParameters(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (timelevels < 2)
  {
      CCTK_PARAMWARN("You have to set 'HydroBase::timelevels to at least 2");
  }
}
*/

extern "C" void RNSReader_Interpolation_C2V(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_RNSReader_Interpolation_C2V;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO("Starting interpolation for ADM variables.");
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

  CCTK_INFO("Done interpolation for ADM variables.");

}

} //namespace
