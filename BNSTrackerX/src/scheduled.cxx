#include <cmath>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "bnstracker.hxx"
#include "separation.hxx"

namespace BNSTrackerX {

CCTK_DEVICE CCTK_HOST void bns_tracker::prep(const cGH *cctkGH) //TODO 
{
  cactus_grid::prep(cctkGH); //TODO
  DECLARE_CCTK_ARGUMENTS;

  const map_cart &m = mcart(); //TODO
  CCTK_REAL *p_glo[]    = {gxx, gxy, gyy, gxz, gyz, gzz};

  rmd_pc.init(m, rho);  //TODO
  alp_pc.init(m, alp);  //TODO
  lfac_pc.init(m, w_lorentz);  //TODO
  glo_pc.init(m, p_glo); //TODO
}

cactus_single<bns_tracker> tracker; //TODO

}// namespace BNSTrackerX



using namespace BNSTrackerX;
using namespace std;

extern "C" int BNSTrackerX_Startup(void)
{
  DECLARE_CCTK_PARAMETERS;
  CCTK_RegisterBanner("BNSTrackerX: Track binary neutron star positions");
  tracker = new bns_tracker(analysis_reflevel, dist_num_points); //TODO
  
  return 0;
}


extern "C" void BNSTrackerX_MoveSurfaces(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  double pi         = 3.141592653589793;
  const int* sidx[3]= {bns_surface_idx_1, bns_surface_idx_2,
                       bns_surface_idx_3};
  bool valid[3]     = {false,false,false};
  const double* rc[3] = {bns_surface_radius_1, bns_surface_radius_2, 
                       bns_surface_radius_3};
  double xc[3]      = {0.0, 0.0, 0.0};
  double yc[3]      = {0.0, 0.0, 0.0};
  double zc[3]      = {0.0, 0.0, 0.0};
  
  if (*bns_merger_stage == BNSTrackerX::INSPIRAL) {
    valid[0] = true;
    xc[0]    = *bns_x_1 + *bns_cms_x;
    yc[0]    = *bns_y_1 + *bns_cms_y;
    
    valid[1] = true;
    xc[1]    = *bns_x_2 + *bns_cms_x;
    yc[1]    = *bns_y_2 + *bns_cms_y;
  } 
  else if (*bns_merger_stage == BNSTrackerX::MERGING) {
    valid[2] = true;
    xc[2]    = *bns_cms_x;
    yc[2]    = *bns_cms_y;
    zc[2]    = move_z ? (*bns_cms_z) : 0.0;
  }
  else if (*bns_merger_stage == BNSTrackerX::COLLAPSING) {
    valid[2] = true;
    xc[2]    = *bh_x;
    yc[2]    = *bh_y;
    zc[2]    = move_z ? (*bh_z) : 0.0;
  }
  
  for (int n=0; n<3; ++n) {
    for (int sl=0; sl<num_surf; sl++) {
      int i = sidx[n][sl];
      if (i < 0) continue;
      double rci = rc[n][sl];
      sf_active[i]    = valid[n] ? 1 : 0;
      sf_valid[i]     = sf_active[i];
      sf_origin_x[i]  = xc[n];  
      sf_origin_y[i]  = yc[n];  
      sf_origin_z[i]  = zc[n];  
        
      sf_area[i]        = 4 * pi * std::pow(rci, 2);
      sf_mean_radius[i] = rci;
      sf_centroid_x[i]  = sf_origin_x[i];
      sf_centroid_y[i]  = sf_origin_y[i];
      sf_centroid_z[i]  = sf_origin_z[i];
        
      sf_quadrupole_xx[i] = 0.0;
      sf_quadrupole_xy[i] = 0.0;
      sf_quadrupole_xz[i] = 0.0;
      sf_quadrupole_yy[i] = 0.0;
      sf_quadrupole_yz[i] = 0.0;
      sf_quadrupole_zz[i] = 0.0;
      sf_min_radius[i]    = rci;
      sf_max_radius[i]    = rci;
      
      sf_min_x[i] = sf_origin_x[i] - rci;
      sf_min_y[i] = sf_origin_y[i] - rci;
      sf_min_z[i] = sf_origin_z[i] - rci;
      sf_max_x[i] = sf_origin_x[i] + rci;
      sf_max_y[i] = sf_origin_y[i] + rci;
      sf_max_z[i] = sf_origin_z[i] + rci;
      
      for (int j=0; j<sf_nphi[n]; ++j) {
        for (int k=0; k<sf_ntheta[n]; ++k) {
          int ind = k + maxntheta * (j + maxnphi * i);
          sf_radius[ind] = rci;
        }
      }
    }
  } 
}

extern "C" void BNSTrackerX_Move_Grids(CCTK_ARGUMENTS)
{  
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  if (*bns_merger_stage == BNSTrackerX::UNKNOWN) {
    CCTK_WARN(1, "BNSTrackerX_MoveGrids: BNS locations not computed yet, skipping.");
    return;
  }
  if (*bns_merger_stage == BNSTrackerX::INSPIRAL) {
    if (*bns_sep_tot < merge_separation) {
      if (! track_only) {
        active[index_region1]     = 0;
        active[index_region2]     = 0;
        active[index_region3]     = 1;
        num_levels[index_region3] += add_levels_post_merge;
        if (follow_hmns) {
          position_x[index_region3] = *bns_cms_x;
          position_y[index_region3] = *bns_cms_y;
          position_z[index_region3] = 0; 
        }        
      }
      *bns_merger_stage         = BNSTrackerX::MERGING;
    } else {
      if (! track_only) {
        position_x[index_region1] = *bns_x_1 + *bns_cms_x;
        position_y[index_region1] = *bns_y_1 + *bns_cms_y;
        position_z[index_region1] = 0;

        position_x[index_region2] = *bns_x_2 + *bns_cms_x;
        position_y[index_region2] = *bns_y_2 + *bns_cms_y;
        position_z[index_region2] = 0;
      }
    }
  }
  if (*bns_merger_stage == BNSTrackerX::MERGING) {
    if (*bns_sep_tot < collapse_separation) {
      if (! track_only) {
        CCTK_INFO("Activating additional refinement level post collapse");
        CCTK_VInfo(CCTK_THORNSTRING, "collapse_separation,add_levels_post_collapse = %26.16e, %d ", collapse_separation,add_levels_post_collapse);
        num_levels[index_region3] += add_levels_post_collapse;
      }
      
      *bns_merger_stage         = BNSTrackerX::COLLAPSING;

      *bns_r_1            = 0;
      *bns_x_1            = 0;
      *bns_y_1            = 0;
      *bns_r_2            = 0;
      *bns_x_2            = 0;
      *bns_y_2            = 0;
      *bns_x_md_1         = 0;
      *bns_y_md_1         = 0;
      *bns_r_md_1         = 0;
      *bns_x_md_2         = 0;
      *bns_y_md_2         = 0;
      *bns_r_md_2         = 0;
      *bns_sep_tot        = 0;
      *bns_coord_sep_bc   = 0;
      *bns_proper_sep_bc  = 0;
      *bns_coord_sep_md   = 0;
      *bns_proper_sep_md  = 0;
      *bns_cms_x          = 0;
      *bns_cms_y          = 0;
      *bns_cms_z          = 0;
      *bh_x               = 0;
      *bh_y               = 0;
      *bh_z               = 0;
      *bns_timestamp  = cctk_time;
    }
    else {
      if (follow_hmns && (! track_only)) {
        position_x[index_region3] = *bns_cms_x;
        position_y[index_region3] = *bns_cms_y;
        position_z[index_region3] = move_z ? (*bns_cms_z) : 0.0; 
      }
    }
  }
  if (*bns_merger_stage == BNSTrackerX::COLLAPSING) {
    if (follow_bh && (! track_only)) {
      position_x[index_region3] = *bh_x;
      position_y[index_region3] = *bh_y;
      position_z[index_region3] = move_z ? (*bh_z) : 0.0; 
    }
  }
  
}

extern "C" void BNSTrackerX_Track_Stars(CCTK_ARGUMENTS)
{  
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  if (cctk_levfac[0] != (1 << analysis_reflevel)) return;

  if (((cctk_iteration % analyze_every) != 0) &&
      (*bns_merger_stage != BNSTrackerX::UNKNOWN)) return;

  if (*bns_merger_stage == BNSTrackerX::UNKNOWN) {
    *bns_merger_stage = BNSTrackerX::INSPIRAL;
  }
  
  try { //TODO
    if (*bns_merger_stage == BNSTrackerX::COLLAPSING) {
      const bns_locations& locs = 
        tracker.prep(CCTK_PASS_CTOC).track_minalp();
      *bh_x = locs.minalp_x;
      *bh_y = locs.minalp_y;
      *bh_z = locs.minalp_z;
    }
    else {
      *bh_x = 0;
      *bh_y = 0;
      *bh_z = 0;
    
      const bns_locations& locs = 
        tracker.prep(CCTK_PASS_CTOC).track_stars(*bns_phase_tot);

      *bns_phase_1    = locs.star1.phase;
      *bns_r_1        = locs.star1.sep;
      *bns_x_1        = locs.star1.x;
      *bns_y_1        = locs.star1.y;
      *bns_phase_2    = locs.star2.phase;
      *bns_r_2        = locs.star2.sep;
      *bns_x_2        = locs.star2.x;
      *bns_y_2        = locs.star2.y;

      *bns_x_md_1     = locs.star1_maxloc.x;
      *bns_y_md_1     = locs.star1_maxloc.y;
      *bns_r_md_1     = locs.star1_maxloc.sep;
      *bns_phase_md_1 = locs.star1_maxloc.phase;

      *bns_x_md_2     = locs.star2_maxloc.x;
      *bns_y_md_2     = locs.star2_maxloc.y;
      *bns_r_md_2     = locs.star2_maxloc.sep;
      *bns_phase_md_2 = locs.star2_maxloc.phase;

      *bns_phase_tot      = locs.avg.phase;
      *bns_sep_tot        = locs.avg.sep;
      *bns_coord_sep_bc   = locs.coord_sep_bc;
      *bns_coord_sep_md   = locs.coord_sep_md;
      *bns_cms_x          = locs.cms_x;
      *bns_cms_y          = locs.cms_y;
      *bns_cms_z          = locs.cms_z;
      *bns_timestamp  = cctk_time;
    }
    
  }
  catch (exception &e) {
    CCTK_WARN(0, e.what());
  }
}

extern "C" void BNSTrackerX_Init_State(CCTK_ARGUMENTS)
{  
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;
  *bns_timestamp = -1;
  *bns_merger_stage = BNSTrackerX::UNKNOWN;
  *bns_phase_tot    = initial_phase;
}

extern "C" void BNSTrackerX_Separation(CCTK_ARGUMENTS)
{  
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  if ((cctk_iteration % analyze_every) != 0) return;

  if (*bns_merger_stage == BNSTrackerX::INSPIRAL) {
    cactus_glop glop(CCTK_PASS_CTOC); //TODO
    *bns_proper_sep_bc = proper_distance(glop, 
                *bns_x_1 + *bns_cms_x, *bns_y_1 + *bns_cms_y, 
                *bns_x_2 + *bns_cms_x, *bns_y_2 + *bns_cms_y, 
                dist_num_points);
    *bns_proper_sep_md = proper_distance(glop, 
                *bns_x_md_1 + *bns_cms_x, *bns_y_md_1 + *bns_cms_y, 
                *bns_x_md_2 + *bns_cms_x, *bns_y_md_2 + *bns_cms_y, 
                dist_num_points);

  } else {
    *bns_proper_sep_bc = 0;
    *bns_proper_sep_md = 0;
  }
}



