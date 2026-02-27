#ifndef INTEGRANDS__VOLUMEINTEGRALS_VACUUMX__
#define INTEGRANDS__VOLUMEINTEGRALS_VACUUMX__

#include <cctk.h>
#include <loop_device.hxx>

#include <cmath>

using namespace Loop;

// Second-order average of vertex-centered grid functions to cell center
// at arbitrary cell index idx.
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
VI_vacuumX_avg_v2c_at(const GF3D2<const T> &gf, const PointDesc &p,
                      const vect<int, dim> &idx) {
  T gf_avg = 0;
  for (int dk = 0; dk < 2; ++dk) {
    for (int dj = 0; dj < 2; ++dj) {
      for (int di = 0; di < 2; ++di) {
        gf_avg += gf(idx + p.DI[0] * di + p.DI[1] * dj + p.DI[2] * dk);
      }
    }
  }
  return gf_avg * T(0.125);
}

/* Integrand for L2 norms */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
VI_vacuumX_L2_integrand(const GF3D2<double> VolIntegrand1,
                        const PointDesc &p,
                        const GF3D2<const CCTK_REAL> f) {
  const CCTK_REAL fL = f(p.I);
  VolIntegrand1(p.I) = fL * fL;
}

/* Center of Lapse: */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
VI_vacuumX_CoL_integrand(const GF3D2<double> VolIntegrand1,
                         const GF3D2<double> VolIntegrand2,
                         const GF3D2<double> VolIntegrand3,
                         const GF3D2<double> VolIntegrand4,
                         const PointDesc &p,
                         const GF3D2<const CCTK_REAL> lapse) {
  const CCTK_REAL lapse_cc = VI_vacuumX_avg_v2c_at(lapse, p, p.I);
  const CCTK_REAL one_minus_lapseL =
      pow(1.0 - lapse_cc, 80); // <- Yields consistent CoL results.
  VolIntegrand1(p.I) = one_minus_lapseL * p.x;
  VolIntegrand2(p.I) = one_minus_lapseL * p.y;
  VolIntegrand3(p.I) = one_minus_lapseL * p.z;
  VolIntegrand4(p.I) = one_minus_lapseL;
}

/* ADM Mass */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
VI_vacuumX_ADM_Mass_integrand_eval_derivs(
    const GF3D2<double> ADM_M_integrand_x,
    const GF3D2<double> ADM_M_integrand_y,
    const GF3D2<double> ADM_M_integrand_z, const PointDesc &p,
    const CCTK_REAL idx, const CCTK_REAL idy, const CCTK_REAL idz,
    const GF3D2<const CCTK_REAL> alp, const GF3D2<const CCTK_REAL> gxx,
    const GF3D2<const CCTK_REAL> gxy, const GF3D2<const CCTK_REAL> gxz,
    const GF3D2<const CCTK_REAL> gyy, const GF3D2<const CCTK_REAL> gyz,
    const GF3D2<const CCTK_REAL> gzz) {
  const CCTK_REAL cm1 = -0.5;

  const auto index = p.I;

  // Read in gamma_{i j} (vertex-centered -> cell-centered average)
  const CCTK_REAL g11L = VI_vacuumX_avg_v2c_at(gxx, p, index);
  const CCTK_REAL g12L = VI_vacuumX_avg_v2c_at(gxy, p, index);
  const CCTK_REAL g13L = VI_vacuumX_avg_v2c_at(gxz, p, index);
  const CCTK_REAL g22L = VI_vacuumX_avg_v2c_at(gyy, p, index);
  const CCTK_REAL g23L = VI_vacuumX_avg_v2c_at(gyz, p, index);
  const CCTK_REAL g33L = VI_vacuumX_avg_v2c_at(gzz, p, index);

  // Metric determinant
  const CCTK_REAL detgL =
      -g13L * g13L * g22L + 2 * g12L * g13L * g23L - g11L * g23L * g23L -
      g12L * g12L * g33L + g11L * g22L * g33L;
  if (!std::isfinite(detgL) || detgL <= 0.0) {
    ADM_M_integrand_x(index) = 0.0;
    ADM_M_integrand_y(index) = 0.0;
    ADM_M_integrand_z(index) = 0.0;
    return;
  }

  CCTK_REAL ginv[3][3];

  // Calculate inverse metric gamma^{i j}
  ginv[0][0] = (g22L * g33L - g23L * g23L) / detgL;
  ginv[0][1] = (g13L * g23L - g12L * g33L) / detgL;
  ginv[0][2] = (g12L * g23L - g13L * g22L) / detgL;
  ginv[1][1] = (g11L * g33L - g13L * g13L) / detgL;
  ginv[1][2] = (g12L * g13L - g11L * g23L) / detgL;
  ginv[2][2] = (g11L * g22L - g12L * g12L) / detgL;

  ginv[1][0] = ginv[0][1];
  ginv[2][0] = ginv[0][2];
  ginv[2][1] = ginv[1][2];

  CCTK_REAL g_d1[3][3][3]; // g_d1[i][j][k] = d_i g_{j k}

  g_d1[0][0][0] =
      (cm1 * (VI_vacuumX_avg_v2c_at(gxx, p, index - p.DI[0]) -
              VI_vacuumX_avg_v2c_at(gxx, p, index + p.DI[0]))) *
      idx;
  g_d1[0][0][1] =
      (cm1 * (VI_vacuumX_avg_v2c_at(gxy, p, index - p.DI[0]) -
              VI_vacuumX_avg_v2c_at(gxy, p, index + p.DI[0]))) *
      idx;
  g_d1[0][0][2] =
      (cm1 * (VI_vacuumX_avg_v2c_at(gxz, p, index - p.DI[0]) -
              VI_vacuumX_avg_v2c_at(gxz, p, index + p.DI[0]))) *
      idx;
  g_d1[0][1][1] =
      (cm1 * (VI_vacuumX_avg_v2c_at(gyy, p, index - p.DI[0]) -
              VI_vacuumX_avg_v2c_at(gyy, p, index + p.DI[0]))) *
      idx;
  g_d1[0][1][2] =
      (cm1 * (VI_vacuumX_avg_v2c_at(gyz, p, index - p.DI[0]) -
              VI_vacuumX_avg_v2c_at(gyz, p, index + p.DI[0]))) *
      idx;
  g_d1[0][2][2] =
      (cm1 * (VI_vacuumX_avg_v2c_at(gzz, p, index - p.DI[0]) -
              VI_vacuumX_avg_v2c_at(gzz, p, index + p.DI[0]))) *
      idx;

  g_d1[0][1][0] = g_d1[0][0][1];
  g_d1[0][2][0] = g_d1[0][0][2];
  g_d1[0][2][1] = g_d1[0][1][2];

  g_d1[1][0][0] =
      (cm1 * (VI_vacuumX_avg_v2c_at(gxx, p, index - p.DI[1]) -
              VI_vacuumX_avg_v2c_at(gxx, p, index + p.DI[1]))) *
      idy;
  g_d1[1][0][1] =
      (cm1 * (VI_vacuumX_avg_v2c_at(gxy, p, index - p.DI[1]) -
              VI_vacuumX_avg_v2c_at(gxy, p, index + p.DI[1]))) *
      idy;
  g_d1[1][0][2] =
      (cm1 * (VI_vacuumX_avg_v2c_at(gxz, p, index - p.DI[1]) -
              VI_vacuumX_avg_v2c_at(gxz, p, index + p.DI[1]))) *
      idy;
  g_d1[1][1][1] =
      (cm1 * (VI_vacuumX_avg_v2c_at(gyy, p, index - p.DI[1]) -
              VI_vacuumX_avg_v2c_at(gyy, p, index + p.DI[1]))) *
      idy;
  g_d1[1][1][2] =
      (cm1 * (VI_vacuumX_avg_v2c_at(gyz, p, index - p.DI[1]) -
              VI_vacuumX_avg_v2c_at(gyz, p, index + p.DI[1]))) *
      idy;
  g_d1[1][2][2] =
      (cm1 * (VI_vacuumX_avg_v2c_at(gzz, p, index - p.DI[1]) -
              VI_vacuumX_avg_v2c_at(gzz, p, index + p.DI[1]))) *
      idy;

  g_d1[1][1][0] = g_d1[1][0][1];
  g_d1[1][2][0] = g_d1[1][0][2];
  g_d1[1][2][1] = g_d1[1][1][2];

  g_d1[2][0][0] =
      (cm1 * (VI_vacuumX_avg_v2c_at(gxx, p, index - p.DI[2]) -
              VI_vacuumX_avg_v2c_at(gxx, p, index + p.DI[2]))) *
      idz;
  g_d1[2][0][1] =
      (cm1 * (VI_vacuumX_avg_v2c_at(gxy, p, index - p.DI[2]) -
              VI_vacuumX_avg_v2c_at(gxy, p, index + p.DI[2]))) *
      idz;
  g_d1[2][0][2] =
      (cm1 * (VI_vacuumX_avg_v2c_at(gxz, p, index - p.DI[2]) -
              VI_vacuumX_avg_v2c_at(gxz, p, index + p.DI[2]))) *
      idz;
  g_d1[2][1][1] =
      (cm1 * (VI_vacuumX_avg_v2c_at(gyy, p, index - p.DI[2]) -
              VI_vacuumX_avg_v2c_at(gyy, p, index + p.DI[2]))) *
      idz;
  g_d1[2][1][2] =
      (cm1 * (VI_vacuumX_avg_v2c_at(gyz, p, index - p.DI[2]) -
              VI_vacuumX_avg_v2c_at(gyz, p, index + p.DI[2]))) *
      idz;
  g_d1[2][2][2] =
      (cm1 * (VI_vacuumX_avg_v2c_at(gzz, p, index - p.DI[2]) -
              VI_vacuumX_avg_v2c_at(gzz, p, index + p.DI[2]))) *
      idz;

  g_d1[2][1][0] = g_d1[2][0][1];
  g_d1[2][2][0] = g_d1[2][0][2];
  g_d1[2][2][1] = g_d1[2][1][2];

  ADM_M_integrand_x(index) = 0;
  ADM_M_integrand_y(index) = 0;
  ADM_M_integrand_z(index) = 0;

  for (int i2 = 0; i2 < 3; i2++) {
    for (int j2 = 0; j2 < 3; j2++) {
      for (int k2 = 0; k2 < 3; k2++) {
        ADM_M_integrand_x(index) +=
            ginv[i2][j2] * ginv[k2][0] *
            (g_d1[j2][i2][k2] - g_d1[k2][i2][j2]);
        ADM_M_integrand_y(index) +=
            ginv[i2][j2] * ginv[k2][1] *
            (g_d1[j2][i2][k2] - g_d1[k2][i2][j2]);
        ADM_M_integrand_z(index) +=
            ginv[i2][j2] * ginv[k2][2] *
            (g_d1[j2][i2][k2] - g_d1[k2][i2][j2]);
      }
    }
  }

  const CCTK_REAL alp_cc = VI_vacuumX_avg_v2c_at(alp, p, index);
  ADM_M_integrand_x(index) *= alp_cc * sqrt(detgL);
  ADM_M_integrand_y(index) *= alp_cc * sqrt(detgL);
  ADM_M_integrand_z(index) *= alp_cc * sqrt(detgL);
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
VI_vacuumX_ADM_Mass_integrand(
    const GF3D2<double> ADM_M_integrand, const PointDesc &p,
    const CCTK_REAL idx, const CCTK_REAL idy, const CCTK_REAL idz,
    const GF3D2<CCTK_REAL> ADM_M_integrand_x,
    const GF3D2<CCTK_REAL> ADM_M_integrand_y,
    const GF3D2<CCTK_REAL> ADM_M_integrand_z) {
  const CCTK_REAL cm1 = -0.5;

  ADM_M_integrand(p.I) =
      0.0625 / M_PI *
      ((cm1 * (ADM_M_integrand_x(p.I - p.DI[0]) -
               ADM_M_integrand_x(p.I + p.DI[0]))) *
           idx +
       (cm1 * (ADM_M_integrand_y(p.I - p.DI[1]) -
               ADM_M_integrand_y(p.I + p.DI[1]))) *
           idy +
       (cm1 * (ADM_M_integrand_z(p.I - p.DI[2]) -
               ADM_M_integrand_z(p.I + p.DI[2]))) *
           idz);
}

/* ADM Momentum */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
VI_vacuumX_ADM_Momentum_integrand_eval_derivs(
    const GF3D2<double> ADM_Py_integrand_x,
    const GF3D2<double> ADM_Py_integrand_y,
    const GF3D2<double> ADM_Py_integrand_z, const PointDesc &p,
    const CCTK_REAL idx, const CCTK_REAL idy, const CCTK_REAL idz,
    const GF3D2<const CCTK_REAL> alp, const GF3D2<const CCTK_REAL> gxx,
    const GF3D2<const CCTK_REAL> gxy, const GF3D2<const CCTK_REAL> gxz,
    const GF3D2<const CCTK_REAL> gyy, const GF3D2<const CCTK_REAL> gyz,
    const GF3D2<const CCTK_REAL> gzz, const GF3D2<const CCTK_REAL> kxx,
    const GF3D2<const CCTK_REAL> kxy, const GF3D2<const CCTK_REAL> kxz,
    const GF3D2<const CCTK_REAL> kyy, const GF3D2<const CCTK_REAL> kyz,
    const GF3D2<const CCTK_REAL> kzz) {

  const auto index = p.I;

  // Read in gamma_{i j} (vertex-centered -> cell-centered average)
  const CCTK_REAL g11L = VI_vacuumX_avg_v2c_at(gxx, p, index);
  const CCTK_REAL g12L = VI_vacuumX_avg_v2c_at(gxy, p, index);
  const CCTK_REAL g13L = VI_vacuumX_avg_v2c_at(gxz, p, index);
  const CCTK_REAL g22L = VI_vacuumX_avg_v2c_at(gyy, p, index);
  const CCTK_REAL g23L = VI_vacuumX_avg_v2c_at(gyz, p, index);
  const CCTK_REAL g33L = VI_vacuumX_avg_v2c_at(gzz, p, index);

  // Metric determinant
  const CCTK_REAL detgL =
      -g13L * g13L * g22L + 2 * g12L * g13L * g23L - g11L * g23L * g23L -
      g12L * g12L * g33L + g11L * g22L * g33L;
  if (!std::isfinite(detgL) || detgL <= 0.0) {
    ADM_Py_integrand_x(index) = 0.0;
    ADM_Py_integrand_y(index) = 0.0;
    ADM_Py_integrand_z(index) = 0.0;
    return;
  }

  CCTK_REAL ginv[3][3];

  // Calculate inverse metric gamma^{i j}
  ginv[0][0] = (g22L * g33L - g23L * g23L) / detgL;
  ginv[0][1] = (g13L * g23L - g12L * g33L) / detgL;
  ginv[0][2] = (g12L * g23L - g13L * g22L) / detgL;
  ginv[1][1] = (g11L * g33L - g13L * g13L) / detgL;
  ginv[1][2] = (g12L * g13L - g11L * g23L) / detgL;
  ginv[2][2] = (g11L * g22L - g12L * g12L) / detgL;

  ginv[1][0] = ginv[0][1];
  ginv[2][0] = ginv[0][2];
  ginv[2][1] = ginv[1][2];

  CCTK_REAL Kdown[3][3];

  // Read in covariant extrinsic curvature K_{i j}
  Kdown[0][0] = VI_vacuumX_avg_v2c_at(kxx, p, index);
  Kdown[0][1] = VI_vacuumX_avg_v2c_at(kxy, p, index);
  Kdown[0][2] = VI_vacuumX_avg_v2c_at(kxz, p, index);
  Kdown[1][1] = VI_vacuumX_avg_v2c_at(kyy, p, index);
  Kdown[1][2] = VI_vacuumX_avg_v2c_at(kyz, p, index);
  Kdown[2][2] = VI_vacuumX_avg_v2c_at(kzz, p, index);

  Kdown[1][0] = Kdown[0][1];
  Kdown[2][0] = Kdown[0][2];
  Kdown[2][1] = Kdown[1][2];

  CCTK_REAL Kup[3][3];

  // Calculate contravariant extrinsic curvature K^{i j}
  for (int i2 = 0; i2 < 3; i2++) {
    for (int j2 = 0; j2 < 3; j2++) {
      Kup[i2][j2] = 0;

      for (int k2 = 0; k2 < 3; k2++) {
        for (int l2 = 0; l2 < 3; l2++) {
          Kup[i2][j2] += ginv[i2][k2] * ginv[j2][l2] * Kdown[k2][l2];
        }
      }
    }
  }

  CCTK_REAL K = 0;

  // Calculate mean curvature K = g^{i j} K_{i j}
  for (int i2 = 0; i2 < 3; i2++) {
    for (int j2 = 0; j2 < 3; j2++) {
      K += ginv[i2][j2] * Kdown[i2][j2];
    }
  }

  CCTK_REAL surface_integrand[3][3];

  const CCTK_REAL alp_cc = VI_vacuumX_avg_v2c_at(alp, p, index);
  for (int i2 = 0; i2 < 3; i2++) {
    for (int j2 = 0; j2 < 3; j2++) {
      surface_integrand[i2][j2] =
          alp_cc * sqrt(detgL) * (Kup[i2][j2] - ginv[i2][j2] * K);
    }
  }

  ADM_Py_integrand_x(index) = surface_integrand[1][0];
  ADM_Py_integrand_y(index) = surface_integrand[1][1];
  ADM_Py_integrand_z(index) = surface_integrand[1][2];
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
VI_vacuumX_ADM_Momentum_integrand(
    const GF3D2<double> ADM_Py_integrand, const PointDesc &p,
    const CCTK_REAL idx, const CCTK_REAL idy, const CCTK_REAL idz,
    const GF3D2<CCTK_REAL> ADM_Py_integrand_x,
    const GF3D2<CCTK_REAL> ADM_Py_integrand_y,
    const GF3D2<CCTK_REAL> ADM_Py_integrand_z) {
  const CCTK_REAL cm1 = -0.5;

  ADM_Py_integrand(p.I) =
      0.125 / M_PI *
      ((cm1 * (ADM_Py_integrand_x(p.I - p.DI[0]) -
               ADM_Py_integrand_x(p.I + p.DI[0]))) *
           idx +
       (cm1 * (ADM_Py_integrand_y(p.I - p.DI[1]) -
               ADM_Py_integrand_y(p.I + p.DI[1]))) *
           idy +
       (cm1 * (ADM_Py_integrand_z(p.I - p.DI[2]) -
               ADM_Py_integrand_z(p.I + p.DI[2]))) *
           idz);
}

/* ADM Angular Momentum */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
VI_vacuumX_ADM_Angular_Momentum_integrand_eval_derivs(
    const GF3D2<double> ADM_Jz_integrand_x,
    const GF3D2<double> ADM_Jz_integrand_y,
    const GF3D2<double> ADM_Jz_integrand_z, const PointDesc &p,
    const CCTK_REAL idx, const CCTK_REAL idy, const CCTK_REAL idz,
    const GF3D2<const CCTK_REAL> alp, const GF3D2<const CCTK_REAL> gxx,
    const GF3D2<const CCTK_REAL> gxy, const GF3D2<const CCTK_REAL> gxz,
    const GF3D2<const CCTK_REAL> gyy, const GF3D2<const CCTK_REAL> gyz,
    const GF3D2<const CCTK_REAL> gzz, const GF3D2<const CCTK_REAL> kxx,
    const GF3D2<const CCTK_REAL> kxy, const GF3D2<const CCTK_REAL> kxz,
    const GF3D2<const CCTK_REAL> kyy, const GF3D2<const CCTK_REAL> kyz,
    const GF3D2<const CCTK_REAL> kzz) {

  const auto index = p.I;

  // Read in gamma_{i j} (vertex-centered -> cell-centered average)
  const CCTK_REAL g11L = VI_vacuumX_avg_v2c_at(gxx, p, index);
  const CCTK_REAL g12L = VI_vacuumX_avg_v2c_at(gxy, p, index);
  const CCTK_REAL g13L = VI_vacuumX_avg_v2c_at(gxz, p, index);
  const CCTK_REAL g22L = VI_vacuumX_avg_v2c_at(gyy, p, index);
  const CCTK_REAL g23L = VI_vacuumX_avg_v2c_at(gyz, p, index);
  const CCTK_REAL g33L = VI_vacuumX_avg_v2c_at(gzz, p, index);

  // Metric determinant
  const CCTK_REAL detgL =
      -g13L * g13L * g22L + 2 * g12L * g13L * g23L - g11L * g23L * g23L -
      g12L * g12L * g33L + g11L * g22L * g33L;
  if (!std::isfinite(detgL) || detgL <= 0.0) {
    ADM_Jz_integrand_x(index) = 0.0;
    ADM_Jz_integrand_y(index) = 0.0;
    ADM_Jz_integrand_z(index) = 0.0;
    return;
  }

  CCTK_REAL ginv[3][3];

  // Calculate inverse metric gamma^{i j}
  ginv[0][0] = (g22L * g33L - g23L * g23L) / detgL;
  ginv[0][1] = (g13L * g23L - g12L * g33L) / detgL;
  ginv[0][2] = (g12L * g23L - g13L * g22L) / detgL;
  ginv[1][1] = (g11L * g33L - g13L * g13L) / detgL;
  ginv[1][2] = (g12L * g13L - g11L * g23L) / detgL;
  ginv[2][2] = (g11L * g22L - g12L * g12L) / detgL;

  ginv[1][0] = ginv[0][1];
  ginv[2][0] = ginv[0][2];
  ginv[2][1] = ginv[1][2];

  CCTK_REAL Kdown[3][3];

  // Read in covariant extrinsic curvature K_{i j}
  Kdown[0][0] = VI_vacuumX_avg_v2c_at(kxx, p, index);
  Kdown[0][1] = VI_vacuumX_avg_v2c_at(kxy, p, index);
  Kdown[0][2] = VI_vacuumX_avg_v2c_at(kxz, p, index);
  Kdown[1][1] = VI_vacuumX_avg_v2c_at(kyy, p, index);
  Kdown[1][2] = VI_vacuumX_avg_v2c_at(kyz, p, index);
  Kdown[2][2] = VI_vacuumX_avg_v2c_at(kzz, p, index);

  Kdown[1][0] = Kdown[0][1];
  Kdown[2][0] = Kdown[0][2];
  Kdown[2][1] = Kdown[1][2];

  CCTK_REAL Kup[3][3];

  // Calculate contravariant extrinsic curvature K^{i j}
  for (int i2 = 0; i2 < 3; i2++) {
    for (int j2 = 0; j2 < 3; j2++) {
      Kup[i2][j2] = 0;

      for (int k2 = 0; k2 < 3; k2++) {
        for (int l2 = 0; l2 < 3; l2++) {
          Kup[i2][j2] += ginv[i2][k2] * ginv[j2][l2] * Kdown[k2][l2];
        }
      }
    }
  }

  CCTK_REAL K = 0;

  // Calculate mean curvature K = g^{i j} K_{i j}
  for (int i2 = 0; i2 < 3; i2++) {
    for (int j2 = 0; j2 < 3; j2++) {
      K += ginv[i2][j2] * Kdown[i2][j2];
    }
  }

  // Levi-Civita tensor
  CCTK_REAL LCT[3][3][3] = {
      {{0, 0, 0}, {0, 0, 1}, {0, -1, 0}},
      {{0, 0, -1}, {0, 0, 0}, {1, 0, 0}},
      {{0, 1, 0}, {-1, 0, 0}, {0, 0, 0}}};

  // Coordinate vector
  CCTK_REAL xx[3] = {p.x, p.y, p.z};

  CCTK_REAL surface_integrand[3][3];

  const CCTK_REAL alp_cc = VI_vacuumX_avg_v2c_at(alp, p, index);
  for (int i2 = 0; i2 < 3; i2++) {
    for (int j2 = 0; j2 < 3; j2++) {
      surface_integrand[i2][j2] = 0;

      for (int k2 = 0; k2 < 3; k2++) {
        for (int l2 = 0; l2 < 3; l2++) {
          surface_integrand[i2][j2] +=
              LCT[i2][k2][l2] * xx[k2] * (Kup[l2][j2] - ginv[l2][j2] * K);
        }
      }

      surface_integrand[i2][j2] *= alp_cc * sqrt(detgL);
    }
  }

  ADM_Jz_integrand_x(index) = surface_integrand[2][0];
  ADM_Jz_integrand_y(index) = surface_integrand[2][1];
  ADM_Jz_integrand_z(index) = surface_integrand[2][2];
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
VI_vacuumX_ADM_Angular_Momentum_integrand(
    const GF3D2<double> ADM_Jz_integrand, const PointDesc &p,
    const CCTK_REAL idx, const CCTK_REAL idy, const CCTK_REAL idz,
    const GF3D2<CCTK_REAL> ADM_Jz_integrand_x,
    const GF3D2<CCTK_REAL> ADM_Jz_integrand_y,
    const GF3D2<CCTK_REAL> ADM_Jz_integrand_z) {
  const CCTK_REAL cm1 = -0.5;

  ADM_Jz_integrand(p.I) =
      0.125 / M_PI *
      ((cm1 * (ADM_Jz_integrand_x(p.I - p.DI[0]) -
               ADM_Jz_integrand_x(p.I + p.DI[0]))) *
           idx +
       (cm1 * (ADM_Jz_integrand_y(p.I - p.DI[1]) -
               ADM_Jz_integrand_y(p.I + p.DI[1]))) *
           idy +
       (cm1 * (ADM_Jz_integrand_z(p.I - p.DI[2]) -
               ADM_Jz_integrand_z(p.I + p.DI[2]))) *
           idz);
}

#endif // INTEGRANDS__VOLUMEINTEGRALS_VACUUMX__
