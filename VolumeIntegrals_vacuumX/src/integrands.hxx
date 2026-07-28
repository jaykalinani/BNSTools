#ifndef INTEGRANDS__VOLUMEINTEGRALS_VACUUMX__
#define INTEGRANDS__VOLUMEINTEGRALS_VACUUMX__

#include <cctk.h>
#include <loop_device.hxx>

#include <cmath>
#include <limits>

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

// Fourth-order first derivative of a vertex-centered field evaluated at a
// cell center. Each stencil value is first averaged from vertices to the
// corresponding cell center.
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
VI_vacuumX_diff1_v2c_at(const GF3D2<const T> &gf, const PointDesc &p,
                        const vect<int, dim> &idx, const int dir,
                        const T inv_dx) {
  const auto di = p.DI[dir];
  return inv_dx *
         (VI_vacuumX_avg_v2c_at(gf, p, idx - 2 * di) -
          8 * VI_vacuumX_avg_v2c_at(gf, p, idx - di) +
          8 * VI_vacuumX_avg_v2c_at(gf, p, idx + di) -
          VI_vacuumX_avg_v2c_at(gf, p, idx + 2 * di)) /
         12;
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
    const CCTK_REAL nan = std::numeric_limits<CCTK_REAL>::quiet_NaN();
    ADM_M_integrand_x(index) = nan;
    ADM_M_integrand_y(index) = nan;
    ADM_M_integrand_z(index) = nan;
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

  g_d1[0][0][0] = VI_vacuumX_diff1_v2c_at(gxx, p, index, 0, idx);
  g_d1[0][0][1] = VI_vacuumX_diff1_v2c_at(gxy, p, index, 0, idx);
  g_d1[0][0][2] = VI_vacuumX_diff1_v2c_at(gxz, p, index, 0, idx);
  g_d1[0][1][1] = VI_vacuumX_diff1_v2c_at(gyy, p, index, 0, idx);
  g_d1[0][1][2] = VI_vacuumX_diff1_v2c_at(gyz, p, index, 0, idx);
  g_d1[0][2][2] = VI_vacuumX_diff1_v2c_at(gzz, p, index, 0, idx);

  g_d1[0][1][0] = g_d1[0][0][1];
  g_d1[0][2][0] = g_d1[0][0][2];
  g_d1[0][2][1] = g_d1[0][1][2];

  g_d1[1][0][0] = VI_vacuumX_diff1_v2c_at(gxx, p, index, 1, idy);
  g_d1[1][0][1] = VI_vacuumX_diff1_v2c_at(gxy, p, index, 1, idy);
  g_d1[1][0][2] = VI_vacuumX_diff1_v2c_at(gxz, p, index, 1, idy);
  g_d1[1][1][1] = VI_vacuumX_diff1_v2c_at(gyy, p, index, 1, idy);
  g_d1[1][1][2] = VI_vacuumX_diff1_v2c_at(gyz, p, index, 1, idy);
  g_d1[1][2][2] = VI_vacuumX_diff1_v2c_at(gzz, p, index, 1, idy);

  g_d1[1][1][0] = g_d1[1][0][1];
  g_d1[1][2][0] = g_d1[1][0][2];
  g_d1[1][2][1] = g_d1[1][1][2];

  g_d1[2][0][0] = VI_vacuumX_diff1_v2c_at(gxx, p, index, 2, idz);
  g_d1[2][0][1] = VI_vacuumX_diff1_v2c_at(gxy, p, index, 2, idz);
  g_d1[2][0][2] = VI_vacuumX_diff1_v2c_at(gxz, p, index, 2, idz);
  g_d1[2][1][1] = VI_vacuumX_diff1_v2c_at(gyy, p, index, 2, idz);
  g_d1[2][1][2] = VI_vacuumX_diff1_v2c_at(gyz, p, index, 2, idz);
  g_d1[2][2][2] = VI_vacuumX_diff1_v2c_at(gzz, p, index, 2, idz);

  g_d1[2][1][0] = g_d1[2][0][1];
  g_d1[2][2][0] = g_d1[2][0][2];
  g_d1[2][2][1] = g_d1[2][1][2];

  const CCTK_REAL gdown[3][3] = {
      {g11L, g12L, g13L}, {g12L, g22L, g23L}, {g13L, g23L, g33L}};

  // Reconstruct a unit-determinant conformal metric from the physical ADM
  // metric. The chi floor only affects pathological points near a puncture,
  // not finite-radius extraction spheres.
  const CCTK_REAL chi_raw = pow(detgL, -1.0 / 3.0);
  const CCTK_REAL chi = chi_raw > 1.0e-4 ? chi_raw : 1.0e-4;
  const CCTK_REAL inv_chi_raw = 1.0 / chi_raw;

  CCTK_REAL h_inv[3][3];
  CCTK_REAL dchi[3] = {0.0, 0.0, 0.0};
  CCTK_REAL dh[3][3][3];

  for (int d = 0; d < 3; ++d) {
    CCTK_REAL trace = 0.0;
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k)
        trace += ginv[j][k] * g_d1[d][j][k];
    dchi[d] = -chi_raw * trace / 3.0;
  }

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      h_inv[i][j] = ginv[i][j] * inv_chi_raw;
      for (int d = 0; d < 3; ++d)
        dh[d][i][j] =
            dchi[d] * gdown[i][j] + chi_raw * g_d1[d][i][j];
    }
  }

  CCTK_REAL contracted_christoffel[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        CCTK_REAL christoffel = 0.0;
        for (int l = 0; l < 3; ++l)
          christoffel +=
              0.5 * h_inv[i][l] *
              (dh[j][l][k] + dh[k][l][j] - dh[l][j][k]);
        contracted_christoffel[i] += h_inv[j][k] * christoffel;
      }
    }
  }

  CCTK_REAL conformal_dchi[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      conformal_dchi[i] += h_inv[i][j] * dchi[j];

  const CCTK_REAL mass_factor = pow(chi, -0.5);
  ADM_M_integrand_x(index) =
      mass_factor *
      (contracted_christoffel[0] + 2.0 / chi * conformal_dchi[0]);
  ADM_M_integrand_y(index) =
      mass_factor *
      (contracted_christoffel[1] + 2.0 / chi * conformal_dchi[1]);
  ADM_M_integrand_z(index) =
      mass_factor *
      (contracted_christoffel[2] + 2.0 / chi * conformal_dchi[2]);
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
  (void)kernel_width;

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
  const CCTK_REAL inv_r = 1.0 / r;
  const CCTK_REAL flux_dot_normal =
      ADM_M_integrand_x(p.I) * dx * inv_r +
      ADM_M_integrand_y(p.I) * dy * inv_r +
      ADM_M_integrand_z(p.I) * dz * inv_r;

  ADM_M_surface_integrand(p.I) = 0.0625 / M_PI * flux_dot_normal;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
VI_vacuumX_surface_normal(const PointDesc &p, const CCTK_REAL center_x,
                          const CCTK_REAL center_y, const CCTK_REAL center_z,
                          const CCTK_REAL radius, CCTK_REAL normal[3]) {
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
  const CCTK_REAL inv_r = 1.0 / r;
  normal[0] = dx * inv_r;
  normal[1] = dy * inv_r;
  normal[2] = dz * inv_r;
  return true;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
VI_vacuumX_ADM_momentum_surface_data(
    const PointDesc &p, const vect<int, dim> &index,
    const GF3D2<const CCTK_REAL> gxx,
    const GF3D2<const CCTK_REAL> gxy, const GF3D2<const CCTK_REAL> gxz,
    const GF3D2<const CCTK_REAL> gyy, const GF3D2<const CCTK_REAL> gyz,
    const GF3D2<const CCTK_REAL> gzz, const GF3D2<const CCTK_REAL> kxx,
    const GF3D2<const CCTK_REAL> kxy, const GF3D2<const CCTK_REAL> kxz,
    const GF3D2<const CCTK_REAL> kyy, const GF3D2<const CCTK_REAL> kyz,
    const GF3D2<const CCTK_REAL> kzz, CCTK_REAL ginv[3][3],
    CCTK_REAL Atilde_UL[3][3], CCTK_REAL &K,
    CCTK_REAL &chi_to_minus_three_halves) {
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

  const CCTK_REAL chi = pow(detgL, -1.0 / 3.0);
  chi_to_minus_three_halves = pow(chi, -1.5);

  // Atilde^m_i = htilde^{ml} Atilde_li. Reconstructing it from physical ADM
  // fields gives K^m_i - delta^m_i K/3, with only the first index raised.
  for (int m = 0; m < 3; ++m) {
    for (int i = 0; i < 3; ++i) {
      Atilde_UL[m][i] = 0.0;
      for (int l = 0; l < 3; ++l)
        Atilde_UL[m][i] += ginv[m][l] * Kdown[l][i];
      if (m == i)
        Atilde_UL[m][i] -= K / 3.0;
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
  (void)kernel_width;
  ADM_Px_surface_integrand(p.I) = 0.0;
  ADM_Py_surface_integrand(p.I) = 0.0;
  ADM_Pz_surface_integrand(p.I) = 0.0;
  VolIntegrand4(p.I) = 0.0;

  CCTK_REAL normal[3];
  if (!VI_vacuumX_surface_normal(p, center_x, center_y, center_z, radius,
                                 normal)) {
    return;
  }

  CCTK_REAL ginv[3][3], Atilde_UL[3][3], K;
  CCTK_REAL chi_to_minus_three_halves;
  if (!VI_vacuumX_ADM_momentum_surface_data(
          p, p.I, gxx, gxy, gxz, gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz,
          kzz, ginv, Atilde_UL, K, chi_to_minus_three_halves)) {
    const CCTK_REAL nan = std::numeric_limits<CCTK_REAL>::quiet_NaN();
    ADM_Px_surface_integrand(p.I) = nan;
    ADM_Py_surface_integrand(p.I) = nan;
    ADM_Pz_surface_integrand(p.I) = nan;
    return;
  }

  CCTK_REAL P[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i)
    for (int m = 0; m < 3; ++m)
      P[i] += chi_to_minus_three_halves *
              (Atilde_UL[m][i] - (m == i ? 2.0 * K / 3.0 : 0.0)) *
              normal[m];

  ADM_Px_surface_integrand(p.I) = 0.125 / M_PI * P[0];
  ADM_Py_surface_integrand(p.I) = 0.125 / M_PI * P[1];
  ADM_Pz_surface_integrand(p.I) = 0.125 / M_PI * P[2];
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
  (void)kernel_width;
  ADM_Jx_surface_integrand(p.I) = 0.0;
  ADM_Jy_surface_integrand(p.I) = 0.0;
  ADM_Jz_surface_integrand(p.I) = 0.0;
  VolIntegrand4(p.I) = 0.0;

  CCTK_REAL normal[3];
  if (!VI_vacuumX_surface_normal(p, center_x, center_y, center_z, radius,
                                 normal)) {
    return;
  }

  CCTK_REAL ginv[3][3], Atilde_UL[3][3], K;
  CCTK_REAL chi_to_minus_three_halves;
  if (!VI_vacuumX_ADM_momentum_surface_data(
          p, p.I, gxx, gxy, gxz, gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz,
          kzz, ginv, Atilde_UL, K, chi_to_minus_three_halves)) {
    const CCTK_REAL nan = std::numeric_limits<CCTK_REAL>::quiet_NaN();
    ADM_Jx_surface_integrand(p.I) = nan;
    ADM_Jy_surface_integrand(p.I) = nan;
    ADM_Jz_surface_integrand(p.I) = nan;
    return;
  }

  const CCTK_REAL xx[3] = {p.x - center_x, p.y - center_y,
                           p.z - center_z};
  const CCTK_REAL LCT[3][3][3] = {
      {{0, 0, 0}, {0, 0, 1}, {0, -1, 0}},
      {{0, 0, -1}, {0, 0, 0}, {1, 0, 0}},
      {{0, 1, 0}, {-1, 0, 0}, {0, 0, 0}}};

  CCTK_REAL J[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i)
    for (int m = 0; m < 3; ++m)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
          J[i] += chi_to_minus_three_halves * LCT[i][j][k] * xx[j] *
                  Atilde_UL[m][k] * normal[m];

  ADM_Jx_surface_integrand(p.I) = 0.125 / M_PI * J[0];
  ADM_Jy_surface_integrand(p.I) = 0.125 / M_PI * J[1];
  ADM_Jz_surface_integrand(p.I) = 0.125 / M_PI * J[2];
}

/* Volume divergence of the same conformal momentum and angular-momentum
 * fluxes used by the direct spherical diagnostics. */
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
VI_vacuumX_ADM_vector_flux_at(
    const PointDesc &p, const vect<int, dim> &index,
    const CCTK_REAL point_x, const CCTK_REAL point_y,
    const CCTK_REAL point_z, const CCTK_REAL center_x,
    const CCTK_REAL center_y, const CCTK_REAL center_z,
    const bool angular_momentum, CCTK_REAL flux[3][3],
    const GF3D2<const CCTK_REAL> gxx, const GF3D2<const CCTK_REAL> gxy,
    const GF3D2<const CCTK_REAL> gxz, const GF3D2<const CCTK_REAL> gyy,
    const GF3D2<const CCTK_REAL> gyz, const GF3D2<const CCTK_REAL> gzz,
    const GF3D2<const CCTK_REAL> kxx, const GF3D2<const CCTK_REAL> kxy,
    const GF3D2<const CCTK_REAL> kxz, const GF3D2<const CCTK_REAL> kyy,
    const GF3D2<const CCTK_REAL> kyz,
    const GF3D2<const CCTK_REAL> kzz) {
  CCTK_REAL ginv[3][3], Atilde_UL[3][3], K;
  CCTK_REAL chi_to_minus_three_halves;
  if (!VI_vacuumX_ADM_momentum_surface_data(
          p, index, gxx, gxy, gxz, gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz,
          kzz, ginv, Atilde_UL, K, chi_to_minus_three_halves)) {
    const CCTK_REAL nan = std::numeric_limits<CCTK_REAL>::quiet_NaN();
    for (int i = 0; i < 3; ++i)
      for (int m = 0; m < 3; ++m)
        flux[i][m] = nan;
    return;
  }

  if (!angular_momentum) {
    for (int i = 0; i < 3; ++i)
      for (int m = 0; m < 3; ++m)
        flux[i][m] =
            0.125 / M_PI * chi_to_minus_three_halves *
            (Atilde_UL[m][i] - (m == i ? 2.0 * K / 3.0 : 0.0));
    return;
  }

  const CCTK_REAL xx[3] = {point_x - center_x, point_y - center_y,
                           point_z - center_z};
  const CCTK_REAL LCT[3][3][3] = {
      {{0, 0, 0}, {0, 0, 1}, {0, -1, 0}},
      {{0, 0, -1}, {0, 0, 0}, {1, 0, 0}},
      {{0, 1, 0}, {-1, 0, 0}, {0, 0, 0}}};
  for (int i = 0; i < 3; ++i) {
    for (int m = 0; m < 3; ++m) {
      flux[i][m] = 0.0;
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
          flux[i][m] +=
              0.125 / M_PI * chi_to_minus_three_halves * LCT[i][j][k] *
              xx[j] * Atilde_UL[m][k];
    }
  }
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
VI_vacuumX_ADM_vector_volume_integrand(
    const GF3D2<double> output_x, const GF3D2<double> output_y,
    const GF3D2<double> output_z, const GF3D2<double> output_unused,
    const PointDesc &p, const CCTK_REAL idx, const CCTK_REAL idy,
    const CCTK_REAL idz, const CCTK_REAL center_x,
    const CCTK_REAL center_y, const CCTK_REAL center_z,
    const bool angular_momentum, const GF3D2<const CCTK_REAL> gxx,
    const GF3D2<const CCTK_REAL> gxy, const GF3D2<const CCTK_REAL> gxz,
    const GF3D2<const CCTK_REAL> gyy, const GF3D2<const CCTK_REAL> gyz,
    const GF3D2<const CCTK_REAL> gzz, const GF3D2<const CCTK_REAL> kxx,
    const GF3D2<const CCTK_REAL> kxy, const GF3D2<const CCTK_REAL> kxz,
    const GF3D2<const CCTK_REAL> kyy, const GF3D2<const CCTK_REAL> kyz,
    const GF3D2<const CCTK_REAL> kzz) {
  const CCTK_REAL inv_dx[3] = {idx, idy, idz};
  const CCTK_REAL spacing[3] = {1.0 / idx, 1.0 / idy, 1.0 / idz};
  const CCTK_REAL base_x[3] = {p.x, p.y, p.z};
  CCTK_REAL divergence[3] = {0.0, 0.0, 0.0};

  for (int d = 0; d < 3; ++d) {
    CCTK_REAL fm2[3][3], fm1[3][3], fp1[3][3], fp2[3][3];
    CCTK_REAL xm2[3] = {base_x[0], base_x[1], base_x[2]};
    CCTK_REAL xm1[3] = {base_x[0], base_x[1], base_x[2]};
    CCTK_REAL xp1[3] = {base_x[0], base_x[1], base_x[2]};
    CCTK_REAL xp2[3] = {base_x[0], base_x[1], base_x[2]};
    xm2[d] -= 2.0 * spacing[d];
    xm1[d] -= spacing[d];
    xp1[d] += spacing[d];
    xp2[d] += 2.0 * spacing[d];

    VI_vacuumX_ADM_vector_flux_at(
        p, p.I - 2 * p.DI[d], xm2[0], xm2[1], xm2[2], center_x, center_y,
        center_z, angular_momentum, fm2, gxx, gxy, gxz, gyy, gyz, gzz, kxx,
        kxy, kxz, kyy, kyz, kzz);
    VI_vacuumX_ADM_vector_flux_at(
        p, p.I - p.DI[d], xm1[0], xm1[1], xm1[2], center_x, center_y,
        center_z, angular_momentum, fm1, gxx, gxy, gxz, gyy, gyz, gzz, kxx,
        kxy, kxz, kyy, kyz, kzz);
    VI_vacuumX_ADM_vector_flux_at(
        p, p.I + p.DI[d], xp1[0], xp1[1], xp1[2], center_x, center_y,
        center_z, angular_momentum, fp1, gxx, gxy, gxz, gyy, gyz, gzz, kxx,
        kxy, kxz, kyy, kyz, kzz);
    VI_vacuumX_ADM_vector_flux_at(
        p, p.I + 2 * p.DI[d], xp2[0], xp2[1], xp2[2], center_x, center_y,
        center_z, angular_momentum, fp2, gxx, gxy, gxz, gyy, gyz, gzz, kxx,
        kxy, kxz, kyy, kyz, kzz);

    for (int i = 0; i < 3; ++i)
      divergence[i] += inv_dx[d] *
                       (fm2[i][d] - 8.0 * fm1[i][d] +
                        8.0 * fp1[i][d] - fp2[i][d]) /
                       12.0;
  }

  output_x(p.I) = divergence[0];
  output_y(p.I) = divergence[1];
  output_z(p.I) = divergence[2];
  output_unused(p.I) = 0.0;
}

#endif // INTEGRANDS__VOLUMEINTEGRALS_VACUUMX__
