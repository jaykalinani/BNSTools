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
                        const GF3D2<const CCTK_REAL> f,
                        const GF3D2<double> VolIntegrand2,
                        const GF3D2<double> VolIntegrand3,
                        const GF3D2<double> VolIntegrand4) {
  const CCTK_REAL fL = f(p.I);
  VolIntegrand1(p.I) = fL * fL;
  VolIntegrand2(p.I) = 0.0;
  VolIntegrand3(p.I) = 0.0;
  VolIntegrand4(p.I) = 0.0;
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
  (void)alp;

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

  // The documented ADM expressions are coordinate surface integrals in
  // asymptotically Cartesian coordinates. After converting them to a volume
  // divergence, do not densitize the flux by lapse or sqrt(det(gamma)).
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

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
VI_vacuumX_ADM_Mass_surface_integrand(
    const GF3D2<double> ADM_M_surface_integrand, const PointDesc &p,
    const CCTK_REAL center_x, const CCTK_REAL center_y,
    const CCTK_REAL center_z, const CCTK_REAL radius,
    const CCTK_REAL kernel_width, const GF3D2<CCTK_REAL> ADM_M_integrand_x,
    const GF3D2<CCTK_REAL> ADM_M_integrand_y,
    const GF3D2<CCTK_REAL> ADM_M_integrand_z) {
  ADM_M_surface_integrand(p.I) = 0.0;

  if (radius <= 0.0) {
    return;
  }

  const CCTK_REAL dx = p.x - center_x;
  const CCTK_REAL dy = p.y - center_y;
  const CCTK_REAL dz = p.z - center_z;
  const CCTK_REAL r2 = dx * dx + dy * dy + dz * dz;
  if (r2 <= 0.0) {
    return;
  }

  const CCTK_REAL r = sqrt(r2);
  const CCTK_REAL max_dx = p.dx > p.dy ? (p.dx > p.dz ? p.dx : p.dz)
                                      : (p.dy > p.dz ? p.dy : p.dz);
  const CCTK_REAL width = kernel_width > 0.0 ? kernel_width : 2.0 * max_dx;
  if (width <= 0.0) {
    return;
  }

  const CCTK_REAL q = fabs(r - radius);
  if (q >= width) {
    return;
  }

  // Triangular approximation to delta(r - R), normalized to unit radial
  // integral. The volume reduction then approximates the spherical surface
  // integral of F^i n_i.
  const CCTK_REAL radial_delta = (1.0 - q / width) / width;
  const CCTK_REAL inv_r = 1.0 / r;
  const CCTK_REAL flux_dot_normal =
      ADM_M_integrand_x(p.I) * dx * inv_r +
      ADM_M_integrand_y(p.I) * dy * inv_r +
      ADM_M_integrand_z(p.I) * dz * inv_r;

  ADM_M_surface_integrand(p.I) =
      0.0625 / M_PI * flux_dot_normal * radial_delta;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
VI_vacuumX_surface_kernel(const PointDesc &p, const CCTK_REAL center_x,
                          const CCTK_REAL center_y,
                          const CCTK_REAL center_z,
                          const CCTK_REAL radius,
                          const CCTK_REAL kernel_width,
                          CCTK_REAL normal[3], CCTK_REAL &radial_delta) {
  radial_delta = 0.0;
  normal[0] = 0.0;
  normal[1] = 0.0;
  normal[2] = 0.0;

  if (radius <= 0.0) {
    return false;
  }

  const CCTK_REAL dx = p.x - center_x;
  const CCTK_REAL dy = p.y - center_y;
  const CCTK_REAL dz = p.z - center_z;
  const CCTK_REAL r2 = dx * dx + dy * dy + dz * dz;
  if (r2 <= 0.0) {
    return false;
  }

  const CCTK_REAL r = sqrt(r2);
  const CCTK_REAL max_dx = p.dx > p.dy ? (p.dx > p.dz ? p.dx : p.dz)
                                      : (p.dy > p.dz ? p.dy : p.dz);
  const CCTK_REAL width = kernel_width > 0.0 ? kernel_width : 2.0 * max_dx;
  if (width <= 0.0) {
    return false;
  }

  const CCTK_REAL q = fabs(r - radius);
  if (q >= width) {
    return false;
  }

  const CCTK_REAL inv_r = 1.0 / r;
  normal[0] = dx * inv_r;
  normal[1] = dy * inv_r;
  normal[2] = dz * inv_r;
  radial_delta = (1.0 - q / width) / width;
  return true;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
VI_vacuumX_ADM_momentum_surface_data(
    const PointDesc &p, const GF3D2<const CCTK_REAL> gxx,
    const GF3D2<const CCTK_REAL> gxy, const GF3D2<const CCTK_REAL> gxz,
    const GF3D2<const CCTK_REAL> gyy, const GF3D2<const CCTK_REAL> gyz,
    const GF3D2<const CCTK_REAL> gzz, const GF3D2<const CCTK_REAL> kxx,
    const GF3D2<const CCTK_REAL> kxy, const GF3D2<const CCTK_REAL> kxz,
    const GF3D2<const CCTK_REAL> kyy, const GF3D2<const CCTK_REAL> kyz,
    const GF3D2<const CCTK_REAL> kzz, CCTK_REAL ginv[3][3],
    CCTK_REAL Kup[3][3], CCTK_REAL &K) {
  const auto index = p.I;

  const CCTK_REAL g11L = VI_vacuumX_avg_v2c_at(gxx, p, index);
  const CCTK_REAL g12L = VI_vacuumX_avg_v2c_at(gxy, p, index);
  const CCTK_REAL g13L = VI_vacuumX_avg_v2c_at(gxz, p, index);
  const CCTK_REAL g22L = VI_vacuumX_avg_v2c_at(gyy, p, index);
  const CCTK_REAL g23L = VI_vacuumX_avg_v2c_at(gyz, p, index);
  const CCTK_REAL g33L = VI_vacuumX_avg_v2c_at(gzz, p, index);

  const CCTK_REAL detgL =
      -g13L * g13L * g22L + 2 * g12L * g13L * g23L - g11L * g23L * g23L -
      g12L * g12L * g33L + g11L * g22L * g33L;
  if (!std::isfinite(detgL) || detgL <= 0.0) {
    return false;
  }

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
  Kdown[0][0] = VI_vacuumX_avg_v2c_at(kxx, p, index);
  Kdown[0][1] = VI_vacuumX_avg_v2c_at(kxy, p, index);
  Kdown[0][2] = VI_vacuumX_avg_v2c_at(kxz, p, index);
  Kdown[1][1] = VI_vacuumX_avg_v2c_at(kyy, p, index);
  Kdown[1][2] = VI_vacuumX_avg_v2c_at(kyz, p, index);
  Kdown[2][2] = VI_vacuumX_avg_v2c_at(kzz, p, index);
  Kdown[1][0] = Kdown[0][1];
  Kdown[2][0] = Kdown[0][2];
  Kdown[2][1] = Kdown[1][2];

  K = 0.0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      K += ginv[i][j] * Kdown[i][j];

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      Kup[i][j] = 0.0;
      for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l)
          Kup[i][j] += ginv[i][k] * ginv[j][l] * Kdown[k][l];
    }
  }

  return true;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
VI_vacuumX_ADM_Momentum_surface_integrand(
    const GF3D2<double> ADM_Px_surface_integrand,
    const GF3D2<double> ADM_Py_surface_integrand,
    const GF3D2<double> ADM_Pz_surface_integrand,
    const GF3D2<double> VolIntegrand4, const PointDesc &p,
    const CCTK_REAL center_x, const CCTK_REAL center_y,
    const CCTK_REAL center_z, const CCTK_REAL radius,
    const CCTK_REAL kernel_width, const GF3D2<const CCTK_REAL> alp,
    const GF3D2<const CCTK_REAL> gxx, const GF3D2<const CCTK_REAL> gxy,
    const GF3D2<const CCTK_REAL> gxz, const GF3D2<const CCTK_REAL> gyy,
    const GF3D2<const CCTK_REAL> gyz, const GF3D2<const CCTK_REAL> gzz,
    const GF3D2<const CCTK_REAL> kxx, const GF3D2<const CCTK_REAL> kxy,
    const GF3D2<const CCTK_REAL> kxz, const GF3D2<const CCTK_REAL> kyy,
    const GF3D2<const CCTK_REAL> kyz, const GF3D2<const CCTK_REAL> kzz) {
  (void)alp;
  ADM_Px_surface_integrand(p.I) = 0.0;
  ADM_Py_surface_integrand(p.I) = 0.0;
  ADM_Pz_surface_integrand(p.I) = 0.0;
  VolIntegrand4(p.I) = 0.0;

  CCTK_REAL normal[3], radial_delta;
  if (!VI_vacuumX_surface_kernel(p, center_x, center_y, center_z, radius,
                                 kernel_width, normal, radial_delta)) {
    return;
  }

  CCTK_REAL ginv[3][3], Kup[3][3], K;
  if (!VI_vacuumX_ADM_momentum_surface_data(
          p, gxx, gxy, gxz, gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz, kzz, ginv,
          Kup, K)) {
    return;
  }

  CCTK_REAL P[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      P[i] += (Kup[i][j] - ginv[i][j] * K) * normal[j];

  ADM_Px_surface_integrand(p.I) = 0.125 / M_PI * P[0] * radial_delta;
  ADM_Py_surface_integrand(p.I) = 0.125 / M_PI * P[1] * radial_delta;
  ADM_Pz_surface_integrand(p.I) = 0.125 / M_PI * P[2] * radial_delta;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
VI_vacuumX_ADM_Angular_Momentum_surface_integrand(
    const GF3D2<double> ADM_Jx_surface_integrand,
    const GF3D2<double> ADM_Jy_surface_integrand,
    const GF3D2<double> ADM_Jz_surface_integrand,
    const GF3D2<double> VolIntegrand4, const PointDesc &p,
    const CCTK_REAL center_x, const CCTK_REAL center_y,
    const CCTK_REAL center_z, const CCTK_REAL radius,
    const CCTK_REAL kernel_width, const GF3D2<const CCTK_REAL> alp,
    const GF3D2<const CCTK_REAL> gxx, const GF3D2<const CCTK_REAL> gxy,
    const GF3D2<const CCTK_REAL> gxz, const GF3D2<const CCTK_REAL> gyy,
    const GF3D2<const CCTK_REAL> gyz, const GF3D2<const CCTK_REAL> gzz,
    const GF3D2<const CCTK_REAL> kxx, const GF3D2<const CCTK_REAL> kxy,
    const GF3D2<const CCTK_REAL> kxz, const GF3D2<const CCTK_REAL> kyy,
    const GF3D2<const CCTK_REAL> kyz, const GF3D2<const CCTK_REAL> kzz) {
  (void)alp;
  ADM_Jx_surface_integrand(p.I) = 0.0;
  ADM_Jy_surface_integrand(p.I) = 0.0;
  ADM_Jz_surface_integrand(p.I) = 0.0;
  VolIntegrand4(p.I) = 0.0;

  CCTK_REAL normal[3], radial_delta;
  if (!VI_vacuumX_surface_kernel(p, center_x, center_y, center_z, radius,
                                 kernel_width, normal, radial_delta)) {
    return;
  }

  CCTK_REAL ginv[3][3], Kup[3][3], K;
  if (!VI_vacuumX_ADM_momentum_surface_data(
          p, gxx, gxy, gxz, gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz, kzz, ginv,
          Kup, K)) {
    return;
  }

  const CCTK_REAL xx[3] = {p.x, p.y, p.z};
  const CCTK_REAL LCT[3][3][3] = {
      {{0, 0, 0}, {0, 0, 1}, {0, -1, 0}},
      {{0, 0, -1}, {0, 0, 0}, {1, 0, 0}},
      {{0, 1, 0}, {-1, 0, 0}, {0, 0, 0}}};

  CCTK_REAL J[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      CCTK_REAL surface_integrand = 0.0;
      for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l)
          surface_integrand +=
              LCT[i][k][l] * xx[k] * (Kup[l][j] - ginv[l][j] * K);
      J[i] += surface_integrand * normal[j];
    }
  }

  ADM_Jx_surface_integrand(p.I) = 0.125 / M_PI * J[0] * radial_delta;
  ADM_Jy_surface_integrand(p.I) = 0.125 / M_PI * J[1] * radial_delta;
  ADM_Jz_surface_integrand(p.I) = 0.125 / M_PI * J[2] * radial_delta;
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
  (void)alp;

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

  for (int i2 = 0; i2 < 3; i2++) {
    for (int j2 = 0; j2 < 3; j2++) {
      surface_integrand[i2][j2] = Kup[i2][j2] - ginv[i2][j2] * K;
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
  (void)alp;

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

  for (int i2 = 0; i2 < 3; i2++) {
    for (int j2 = 0; j2 < 3; j2++) {
      surface_integrand[i2][j2] = 0;

      for (int k2 = 0; k2 < 3; k2++) {
        for (int l2 = 0; l2 < 3; l2++) {
          surface_integrand[i2][j2] +=
              LCT[i2][k2][l2] * xx[k2] * (Kup[l2][j2] - ginv[l2][j2] * K);
        }
      }

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
