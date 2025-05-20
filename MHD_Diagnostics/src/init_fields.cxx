#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cmath>
#include <iostream>
#include <loop_device.hxx>

namespace MHD_Diagnostics {

using namespace Loop;
extern "C" void init_fields(CCTK_ARGUMENTS) 
{
  DECLARE_CCTK_ARGUMENTSX_init_fields;  // Declare all grid functions from interface.ccl
  DECLARE_CCTK_PARAMETERS; // Declare all parameters from param.ccl

  CCTK_INFO ("Setting  diagnostic fields to zero");
  grid.loop_all_device<1, 1, 1>(
    grid.nghostzones,
    [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

        smallbt(p.I) = 0.0;
	smallbx(p.I) = 0.0;
	smallby(p.I) = 0.0;
	smallbz(p.I) = 0.0;
	smallb2(p.I) = 0.0;
	Poynx(p.I)   = 0.0;
	Poyny(p.I)   = 0.0;
	Poynz(p.I)   = 0.0;
	Poyn2x(p.I)   = 0.0;
	Poyn2y(p.I)   = 0.0;
	Poyn2z(p.I)   = 0.0;
	sqrt_gamma(p.I) = 0.0;
	minus_one_minus_u_0(p.I) = 0.0;

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

	normB(p.I) = 0.0;
	normcurlB(p.I) = 0.0;
	a4sq(p.I) = 0.0;

  velxold(p.I) = 0.0;
  velyold(p.I) = 0.0;
  velzold(p.I) = 0.0;
  w_lorentzold(p.I) = 0.0;

	/*
	expansion_4(p.I)               = 0.0;
	shear_tt(p.I)                  = 0.0;
	shear_tx(p.I)                  = 0.0;
	shear_ty(p.I)                  = 0.0;
	shear_tz(p.I)                  = 0.0;
	shear_xx(p.I)                  = 0.0;
	shear_xy(p.I)                  = 0.0;
	shear_xz(p.I)                  = 0.0;
	shear_yy(p.I)                  = 0.0;
	shear_yz(p.I)                  = 0.0;
        shear_zz(p.I)                  = 0.0;
	omega_tx(p.I)                  = 0.0;
	omega_ty(p.I)                  = 0.0;
	omega_tz(p.I)                  = 0.0;
	omega_xy(p.I)                  = 0.0;
	omega_xz(p.I)                  = 0.0;
        omega_yz(p.I)                  = 0.0;
	acc_t(p.I)                     = 0.0;
	acc_x(p.I)                     = 0.0;
	acc_y(p.I)                     = 0.0;
	acc_z(p.I)                     = 0.0;
*/

   });

   CCTK_INFO ("Setting diagnostic fields to zero was successful");
}

} // namespace
