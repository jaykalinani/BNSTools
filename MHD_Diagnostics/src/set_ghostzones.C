#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cmath>
#include <iostream>

extern "C" void set_ghostzones(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;  // Declare all grid functions from interface.ccl
  DECLARE_CCTK_PARAMETERS; // Declare all parameters from param.ccl

// Ghost zones in x-y plane left
#pragma omp parallel for
  for (int k = 0; k < cctk_nghostzones[2]; k++) {
    for (int j = 0; j < cctk_lsh[1]; j++) {
      for (int i = 0; i < cctk_lsh[0]; i++) {
        const size_t index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        const size_t index_copy =
            CCTK_GFINDEX3D(cctkGH, i, j, cctk_nghostzones[2]);

        smallbt[index] = smallbt[index_copy];
        smallbx[index] = smallbx[index_copy];
        smallby[index] = smallby[index_copy];
        smallbz[index] = smallbz[index_copy];
        smallb2[index] = smallb2[index_copy];
        Poynx[index] = Poynx[index_copy];
        Poyny[index] = Poyny[index_copy];
        Poynz[index] = Poynz[index_copy];
        Poyn2x[index] = Poyn2x[index_copy];
        Poyn2y[index] = Poyn2y[index_copy];
        Poyn2z[index] = Poyn2z[index_copy];
        sqrt_gamma[index] = sqrt_gamma[index_copy];

        minus_one_minus_u_0[index] = minus_one_minus_u_0[index_copy];

        expansion_scalar[index] = expansion_scalar[index_copy];
        expansion_no_time_der[index] = expansion_no_time_der[index_copy];
        shear_spatial_tensor_xx[index] = shear_spatial_tensor_xx[index_copy];
        shear_spatial_tensor_xy[index] = shear_spatial_tensor_xy[index_copy];
        shear_spatial_tensor_xz[index] = shear_spatial_tensor_xz[index_copy];
        shear_spatial_tensor_yy[index] = shear_spatial_tensor_yy[index_copy];
        shear_spatial_tensor_yz[index] = shear_spatial_tensor_yz[index_copy];
        shear_spatial_tensor_zz[index] = shear_spatial_tensor_zz[index_copy];
        kin_vorticity_spatial_xy[index] = kin_vorticity_spatial_xy[index_copy];
        kin_vorticity_spatial_xz[index] = kin_vorticity_spatial_xz[index_copy];
        kin_vorticity_spatial_yz[index] = kin_vorticity_spatial_yz[index_copy];
        kin_acceleration_spatial_x[index] =
            kin_acceleration_spatial_x[index_copy];
        kin_acceleration_spatial_y[index] =
            kin_acceleration_spatial_y[index_copy];
        kin_acceleration_spatial_z[index] =
            kin_acceleration_spatial_z[index_copy];

        sigma4bb[index] = sigma4bb[index_copy];
        sigma4Ut[index] = sigma4Ut[index_copy];
        sigma4Ux[index] = sigma4Ux[index_copy];
        sigma4Uy[index] = sigma4Uy[index_copy];
        sigma4Uz[index] = sigma4Uz[index_copy];
        sigma4Trace[index] = sigma4Trace[index_copy];
        omega4Ut[index] = omega4Ut[index_copy];
        omega4Ux[index] = omega4Ux[index_copy];
        omega4Uy[index] = omega4Uy[index_copy];
        omega4Uz[index] = omega4Uz[index_copy];

        normB[index] = normB[index_copy];
        normcurlB[index] = normcurlB[index_copy];
        a4sq[index] = a4sq[index_copy];
      }
    }
  }

// Ghost zones in x-y plane right
#pragma omp parallel for
  for (int k = cctk_lsh[2] - cctk_nghostzones[2]; k < cctk_lsh[2]; k++) {
    for (int j = 0; j < cctk_lsh[1]; j++) {
      for (int i = 0; i < cctk_lsh[0]; i++) {
        const size_t index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        const size_t index_copy =
            CCTK_GFINDEX3D(cctkGH, i, j, cctk_lsh[2] - cctk_nghostzones[2] - 1);

        smallbt[index] = smallbt[index_copy];
        smallbx[index] = smallbx[index_copy];
        smallby[index] = smallby[index_copy];
        smallbz[index] = smallbz[index_copy];
        smallb2[index] = smallb2[index_copy];
        Poynx[index] = Poynx[index_copy];
        Poyny[index] = Poyny[index_copy];
        Poynz[index] = Poynz[index_copy];
        Poyn2x[index] = Poyn2x[index_copy];
        Poyn2y[index] = Poyn2y[index_copy];
        Poyn2z[index] = Poyn2z[index_copy];
        sqrt_gamma[index] = sqrt_gamma[index_copy];

        minus_one_minus_u_0[index] = minus_one_minus_u_0[index_copy];

        expansion_scalar[index] = expansion_scalar[index_copy];
        expansion_no_time_der[index] = expansion_no_time_der[index_copy];
        shear_spatial_tensor_xx[index] = shear_spatial_tensor_xx[index_copy];
        shear_spatial_tensor_xy[index] = shear_spatial_tensor_xy[index_copy];
        shear_spatial_tensor_xz[index] = shear_spatial_tensor_xz[index_copy];
        shear_spatial_tensor_yy[index] = shear_spatial_tensor_yy[index_copy];
        shear_spatial_tensor_yz[index] = shear_spatial_tensor_yz[index_copy];
        shear_spatial_tensor_zz[index] = shear_spatial_tensor_zz[index_copy];
        kin_vorticity_spatial_xy[index] = kin_vorticity_spatial_xy[index_copy];
        kin_vorticity_spatial_xz[index] = kin_vorticity_spatial_xz[index_copy];
        kin_vorticity_spatial_yz[index] = kin_vorticity_spatial_yz[index_copy];
        kin_acceleration_spatial_x[index] =
            kin_acceleration_spatial_x[index_copy];
        kin_acceleration_spatial_y[index] =
            kin_acceleration_spatial_y[index_copy];
        kin_acceleration_spatial_z[index] =
            kin_acceleration_spatial_z[index_copy];

        sigma4bb[index] = sigma4bb[index_copy];
        sigma4Ut[index] = sigma4Ut[index_copy];
        sigma4Ux[index] = sigma4Ux[index_copy];
        sigma4Uy[index] = sigma4Uy[index_copy];
        sigma4Uz[index] = sigma4Uz[index_copy];
        sigma4Trace[index] = sigma4Trace[index_copy];
        omega4Ut[index] = omega4Ut[index_copy];
        omega4Ux[index] = omega4Ux[index_copy];
        omega4Uy[index] = omega4Uy[index_copy];
        omega4Uz[index] = omega4Uz[index_copy];

        normB[index] = normB[index_copy];
        normcurlB[index] = normcurlB[index_copy];
        a4sq[index] = a4sq[index_copy];
      }
    }
  }

// Ghost zones in x-z plane left
#pragma omp parallel for
  for (int k = 0; k < cctk_lsh[2]; k++) {
    for (int j = 0; j < cctk_nghostzones[1]; j++) {
      for (int i = 0; i < cctk_lsh[0]; i++) {
        const size_t index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        const size_t index_copy =
            CCTK_GFINDEX3D(cctkGH, i, cctk_nghostzones[1], k);

        smallbt[index] = smallbt[index_copy];
        smallbx[index] = smallbx[index_copy];
        smallby[index] = smallby[index_copy];
        smallbz[index] = smallbz[index_copy];
        smallb2[index] = smallb2[index_copy];
        Poynx[index] = Poynx[index_copy];
        Poyny[index] = Poyny[index_copy];
        Poynz[index] = Poynz[index_copy];
        Poyn2x[index] = Poyn2x[index_copy];
        Poyn2y[index] = Poyn2y[index_copy];
        Poyn2z[index] = Poyn2z[index_copy];
        sqrt_gamma[index] = sqrt_gamma[index_copy];

        minus_one_minus_u_0[index] = minus_one_minus_u_0[index_copy];

        expansion_scalar[index] = expansion_scalar[index_copy];
        expansion_no_time_der[index] = expansion_no_time_der[index_copy];
        shear_spatial_tensor_xx[index] = shear_spatial_tensor_xx[index_copy];
        shear_spatial_tensor_xy[index] = shear_spatial_tensor_xy[index_copy];
        shear_spatial_tensor_xz[index] = shear_spatial_tensor_xz[index_copy];
        shear_spatial_tensor_yy[index] = shear_spatial_tensor_yy[index_copy];
        shear_spatial_tensor_yz[index] = shear_spatial_tensor_yz[index_copy];
        shear_spatial_tensor_zz[index] = shear_spatial_tensor_zz[index_copy];
        kin_vorticity_spatial_xy[index] = kin_vorticity_spatial_xy[index_copy];
        kin_vorticity_spatial_xz[index] = kin_vorticity_spatial_xz[index_copy];
        kin_vorticity_spatial_yz[index] = kin_vorticity_spatial_yz[index_copy];
        kin_acceleration_spatial_x[index] =
            kin_acceleration_spatial_x[index_copy];
        kin_acceleration_spatial_y[index] =
            kin_acceleration_spatial_y[index_copy];
        kin_acceleration_spatial_z[index] =
            kin_acceleration_spatial_z[index_copy];

        sigma4bb[index] = sigma4bb[index_copy];
        sigma4Ut[index] = sigma4Ut[index_copy];
        sigma4Ux[index] = sigma4Ux[index_copy];
        sigma4Uy[index] = sigma4Uy[index_copy];
        sigma4Uz[index] = sigma4Uz[index_copy];
        sigma4Trace[index] = sigma4Trace[index_copy];
        omega4Ut[index] = omega4Ut[index_copy];
        omega4Ux[index] = omega4Ux[index_copy];
        omega4Uy[index] = omega4Uy[index_copy];
        omega4Uz[index] = omega4Uz[index_copy];

        normB[index] = normB[index_copy];
        normcurlB[index] = normcurlB[index_copy];
        a4sq[index] = a4sq[index_copy];
      }
    }
  }

// Ghost zones in x-z plane right
#pragma omp parallel for
  for (int k = 0; k < cctk_lsh[2]; k++) {
    for (int j = cctk_lsh[1] - cctk_nghostzones[1]; j < cctk_lsh[1]; j++) {
      for (int i = 0; i < cctk_lsh[0]; i++) {
        const size_t index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        const size_t index_copy =
            CCTK_GFINDEX3D(cctkGH, i, cctk_lsh[1] - cctk_nghostzones[1] - 1, k);

        smallbt[index] = smallbt[index_copy];
        smallbx[index] = smallbx[index_copy];
        smallby[index] = smallby[index_copy];
        smallbz[index] = smallbz[index_copy];
        smallb2[index] = smallb2[index_copy];
        Poynx[index] = Poynx[index_copy];
        Poyny[index] = Poyny[index_copy];
        Poynz[index] = Poynz[index_copy];
        Poyn2x[index] = Poyn2x[index_copy];
        Poyn2y[index] = Poyn2y[index_copy];
        Poyn2z[index] = Poyn2z[index_copy];
        sqrt_gamma[index] = sqrt_gamma[index_copy];

        minus_one_minus_u_0[index] = minus_one_minus_u_0[index_copy];

        expansion_scalar[index] = expansion_scalar[index_copy];
        expansion_no_time_der[index] = expansion_no_time_der[index_copy];
        shear_spatial_tensor_xx[index] = shear_spatial_tensor_xx[index_copy];
        shear_spatial_tensor_xy[index] = shear_spatial_tensor_xy[index_copy];
        shear_spatial_tensor_xz[index] = shear_spatial_tensor_xz[index_copy];
        shear_spatial_tensor_yy[index] = shear_spatial_tensor_yy[index_copy];
        shear_spatial_tensor_yz[index] = shear_spatial_tensor_yz[index_copy];
        shear_spatial_tensor_zz[index] = shear_spatial_tensor_zz[index_copy];
        kin_vorticity_spatial_xy[index] = kin_vorticity_spatial_xy[index_copy];
        kin_vorticity_spatial_xz[index] = kin_vorticity_spatial_xz[index_copy];
        kin_vorticity_spatial_yz[index] = kin_vorticity_spatial_yz[index_copy];
        kin_acceleration_spatial_x[index] =
            kin_acceleration_spatial_x[index_copy];
        kin_acceleration_spatial_y[index] =
            kin_acceleration_spatial_y[index_copy];
        kin_acceleration_spatial_z[index] =
            kin_acceleration_spatial_z[index_copy];

        sigma4bb[index] = sigma4bb[index_copy];
        sigma4Ut[index] = sigma4Ut[index_copy];
        sigma4Ux[index] = sigma4Ux[index_copy];
        sigma4Uy[index] = sigma4Uy[index_copy];
        sigma4Uz[index] = sigma4Uz[index_copy];
        sigma4Trace[index] = sigma4Trace[index_copy];
        omega4Ut[index] = omega4Ut[index_copy];
        omega4Ux[index] = omega4Ux[index_copy];
        omega4Uy[index] = omega4Uy[index_copy];
        omega4Uz[index] = omega4Uz[index_copy];

        normB[index] = normB[index_copy];
        normcurlB[index] = normcurlB[index_copy];
        a4sq[index] = a4sq[index_copy];
      }
    }
  }

// Ghost zones in y-z plane left
#pragma omp parallel for
  for (int k = 0; k < cctk_lsh[2]; k++) {
    for (int j = 0; j < cctk_lsh[1]; j++) {
      for (int i = 0; i < cctk_nghostzones[0]; i++) {
        const size_t index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        const size_t index_copy =
            CCTK_GFINDEX3D(cctkGH, cctk_nghostzones[0], j, k);

        smallbt[index] = smallbt[index_copy];
        smallbx[index] = smallbx[index_copy];
        smallby[index] = smallby[index_copy];
        smallbz[index] = smallbz[index_copy];
        smallb2[index] = smallb2[index_copy];
        Poynx[index] = Poynx[index_copy];
        Poyny[index] = Poyny[index_copy];
        Poynz[index] = Poynz[index_copy];
        Poyn2x[index] = Poyn2x[index_copy];
        Poyn2y[index] = Poyn2y[index_copy];
        Poyn2z[index] = Poyn2z[index_copy];
        sqrt_gamma[index] = sqrt_gamma[index_copy];

        minus_one_minus_u_0[index] = minus_one_minus_u_0[index_copy];

        expansion_scalar[index] = expansion_scalar[index_copy];
        expansion_no_time_der[index] = expansion_no_time_der[index_copy];
        shear_spatial_tensor_xx[index] = shear_spatial_tensor_xx[index_copy];
        shear_spatial_tensor_xy[index] = shear_spatial_tensor_xy[index_copy];
        shear_spatial_tensor_xz[index] = shear_spatial_tensor_xz[index_copy];
        shear_spatial_tensor_yy[index] = shear_spatial_tensor_yy[index_copy];
        shear_spatial_tensor_yz[index] = shear_spatial_tensor_yz[index_copy];
        shear_spatial_tensor_zz[index] = shear_spatial_tensor_zz[index_copy];
        kin_vorticity_spatial_xy[index] = kin_vorticity_spatial_xy[index_copy];
        kin_vorticity_spatial_xz[index] = kin_vorticity_spatial_xz[index_copy];
        kin_vorticity_spatial_yz[index] = kin_vorticity_spatial_yz[index_copy];
        kin_acceleration_spatial_x[index] =
            kin_acceleration_spatial_x[index_copy];
        kin_acceleration_spatial_y[index] =
            kin_acceleration_spatial_y[index_copy];
        kin_acceleration_spatial_z[index] =
            kin_acceleration_spatial_z[index_copy];

        sigma4bb[index] = sigma4bb[index_copy];
        sigma4Ut[index] = sigma4Ut[index_copy];
        sigma4Ux[index] = sigma4Ux[index_copy];
        sigma4Uy[index] = sigma4Uy[index_copy];
        sigma4Uz[index] = sigma4Uz[index_copy];
        sigma4Trace[index] = sigma4Trace[index_copy];
        omega4Ut[index] = omega4Ut[index_copy];
        omega4Ux[index] = omega4Ux[index_copy];
        omega4Uy[index] = omega4Uy[index_copy];
        omega4Uz[index] = omega4Uz[index_copy];

        normB[index] = normB[index_copy];
        normcurlB[index] = normcurlB[index_copy];
        a4sq[index] = a4sq[index_copy];
      }
    }
  }

// Ghost zones in y-z plane right
#pragma omp parallel for
  for (int k = 0; k < cctk_lsh[2]; k++) {
    for (int j = 0; j < cctk_lsh[1]; j++) {
      for (int i = cctk_lsh[0] - cctk_nghostzones[0]; i < cctk_lsh[0]; i++) {
        const size_t index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        const size_t index_copy =
            CCTK_GFINDEX3D(cctkGH, cctk_lsh[0] - cctk_nghostzones[0] - 1, j, k);

        smallbt[index] = smallbt[index_copy];
        smallbx[index] = smallbx[index_copy];
        smallby[index] = smallby[index_copy];
        smallbz[index] = smallbz[index_copy];
        smallb2[index] = smallb2[index_copy];
        Poynx[index] = Poynx[index_copy];
        Poyny[index] = Poyny[index_copy];
        Poynz[index] = Poynz[index_copy];
        Poyn2x[index] = Poyn2x[index_copy];
        Poyn2y[index] = Poyn2y[index_copy];
        Poyn2z[index] = Poyn2z[index_copy];
        sqrt_gamma[index] = sqrt_gamma[index_copy];

        minus_one_minus_u_0[index] = minus_one_minus_u_0[index_copy];

        expansion_scalar[index] = expansion_scalar[index_copy];
        expansion_no_time_der[index] = expansion_no_time_der[index_copy];
        shear_spatial_tensor_xx[index] = shear_spatial_tensor_xx[index_copy];
        shear_spatial_tensor_xy[index] = shear_spatial_tensor_xy[index_copy];
        shear_spatial_tensor_xz[index] = shear_spatial_tensor_xz[index_copy];
        shear_spatial_tensor_yy[index] = shear_spatial_tensor_yy[index_copy];
        shear_spatial_tensor_yz[index] = shear_spatial_tensor_yz[index_copy];
        shear_spatial_tensor_zz[index] = shear_spatial_tensor_zz[index_copy];
        kin_vorticity_spatial_xy[index] = kin_vorticity_spatial_xy[index_copy];
        kin_vorticity_spatial_xz[index] = kin_vorticity_spatial_xz[index_copy];
        kin_vorticity_spatial_yz[index] = kin_vorticity_spatial_yz[index_copy];
        kin_acceleration_spatial_x[index] =
            kin_acceleration_spatial_x[index_copy];
        kin_acceleration_spatial_y[index] =
            kin_acceleration_spatial_y[index_copy];
        kin_acceleration_spatial_z[index] =
            kin_acceleration_spatial_z[index_copy];

        sigma4bb[index] = sigma4bb[index_copy];
        sigma4Ut[index] = sigma4Ut[index_copy];
        sigma4Ux[index] = sigma4Ux[index_copy];
        sigma4Uy[index] = sigma4Uy[index_copy];
        sigma4Uz[index] = sigma4Uz[index_copy];
        sigma4Trace[index] = sigma4Trace[index_copy];
        omega4Ut[index] = omega4Ut[index_copy];
        omega4Ux[index] = omega4Ux[index_copy];
        omega4Uy[index] = omega4Uy[index_copy];
        omega4Uz[index] = omega4Uz[index_copy];

        normB[index] = normB[index_copy];
        normcurlB[index] = normcurlB[index_copy];
        a4sq[index] = a4sq[index_copy];
      }
    }
  }
}
