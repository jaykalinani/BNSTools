#ifndef VI_UTILS_HXX
#define VI_UTILS_HXX

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <mat.hxx>
#include <vec.hxx>
#include <sum.hxx>
#include <simd.hxx>

#include <algorithm>
#include <array>
#include <cmath>

using namespace std;
using namespace Loop;
using namespace Arith;

// Second-order average of vertex-centered grid functions to cell center
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_avg_v2c(const GF3D2<const T> &gf, const PointDesc &p) {
  T gf_avg = 0.0;

  for (int dk = 0; dk < 2; ++dk) {
    for (int dj = 0; dj < 2; ++dj) {
      for (int di = 0; di < 2; ++di) {
        gf_avg += gf(p.I + p.DI[0] * di + p.DI[1] * dj + p.DI[2] * dk);
      }
    }
  }
  return gf_avg / 8.0;
}

#endif
