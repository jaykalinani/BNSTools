#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>

namespace MHD_Diagnostics {

using namespace Loop;

extern "C" void MHD_Diagnostics_Sync_Recovered(CCTK_ARGUMENTS) {
  // The scheduled SYNC clause does the work.
}

extern "C" void MHD_Diagnostics_Reset_VelsOld_PostRegrid(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MHD_Diagnostics_Reset_VelsOld_PostRegrid;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        velxold(p.I) = velx(p.I);
        velyold(p.I) = vely(p.I);
        velzold(p.I) = velz(p.I);
        w_lorentzold(p.I) = w_lorentz(p.I);
      });
}

} // namespace MHD_Diagnostics
