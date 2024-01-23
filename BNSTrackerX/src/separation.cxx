#include "separation.hxx"
#include <cassert>
#include <cmath>
#include <vector>
#include <string>

namespace BNSTrackerX {

CCTK_DEVICE CCTK_HOST CCTK_REAL integrate_regspaced(const std::vector<CCTK_REAL>& ds, const CCTK_REAL dx)
{
  CCTK_REAL s(0.0);
  for (unsigned int k=0; k<ds.size()-1; k++) {
    s += 0.5*(ds[k]+ds[k+1]);
  }
  return dx*s;
}

CCTK_DEVICE CCTK_HOST  CCTK_REAL proper_distance(const cactus_glop& glop, 
                        const CCTK_REAL x0, const CCTK_REAL y0, 
                        const CCTK_REAL x1, const CCTK_REAL y1, int num_points) 
{
  enum {X=0, Y=1, Z=2, XX=0, XY=1, XZ=2, YY=3, YZ=4, ZZ=5, NUMG=6};

  assert(num_points>10 && num_points<10000);
  vec_u pos0, pos1;
  pos0(X)=x0; pos0(Y)=y0; pos0(Z)=0.0;
  pos1(X)=x1; pos1(Y)=y1; pos1(Z)=0.0;

  const CCTK_REAL dl = 1.0 / (num_points - 1);
  const vec_u   dp = (pos1-pos0) * dl;

  const std::string glo_names[NUMG] = {
    "ADMBase::gxx","ADMBase::gxy","ADMBase::gxz",
    "ADMBase::gyy","ADMBase::gyz","ADMBase::gzz"};
 
  std::vector<CCTK_REAL> glo[NUMG], ds(num_points);
  std::vector<vec_u> samp_coord(num_points);

  for (int k=0; k<num_points; k++) {
    CCTK_REAL l = dl * CCTK_REAL(k);
    samp_coord[k] = (1.0 - l) * pos0 + l * pos1;
  }

  for (int i=0; i<NUMG; ++i) {
    var_index gi(glo_names[i]);
    glop.interpolate(gi, samp_coord, glo[i]);
  }
  
  for (int k=0; k<num_points; k++) {
    ds[k] = sqrt( glo[XX][k] * dp(X) * dp(X) + 
                  glo[YY][k] * dp(Y) * dp(Y) +
                  glo[ZZ][k] * dp(Z) * dp(Z) +
                  2.0 * glo[XY][k] * dp(X) * dp(Y) +  
                  2.0 * glo[XZ][k] * dp(X) * dp(Z) + 
                  2.0 * glo[YZ][k] * dp(Y) * dp(Z) ) / dl;
  }

  return integrate_regspaced(ds, dl);
}

}

