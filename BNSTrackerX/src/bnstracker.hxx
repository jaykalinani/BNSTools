#ifndef BNSTRACKERX_HXX
#define BNSTRACKERX_HXX

namespace BNSTrackerX {

enum merger_state_t {UNKNOWN=-1, INSPIRAL=0, MERGING=1, COLLAPSING=2};

struct orb_param {
  CCTK_REAL phase;
  CCTK_REAL sep;
  CCTK_REAL x;
  CCTK_REAL y;
  CCTK_DEVICE CCTK_HOST void set_phase_sep(CCTK_REAL phi, CCTK_REAL r);
};

struct bns_locations {
  orb_param star1;
  orb_param star1_maxloc;
  orb_param star2;
  orb_param star2_maxloc;
  orb_param avg;
  CCTK_REAL coord_sep_bc;
  CCTK_REAL coord_sep_md;
  CCTK_REAL cms_x;
  CCTK_REAL cms_y;
  CCTK_REAL cms_z;
  CCTK_REAL minalp_x;
  CCTK_REAL minalp_y;
  CCTK_REAL minalp_z;
};

class bns_tracker : bns_locations {
  //TODO: correct the defition type
  scalar_pc rmd_pc;
  scalar_pc alp_pc;
  scalar_pc lfac_pc;
  mats_l_pc glo_pc;

  int use_reflevel;
  int dist_num_points;

  CCTK_DEVICE CCTK_HOST static CCTK_REAL remove_phase_jump(CCTK_REAL old_phase, CCTK_REAL new_phase);
  CCTK_DEVICE CCTK_HOST void get_dens(const gpos& i, CCTK_REAL& rmd, CCTK_REAL& crmd);
  CCTK_DEVICE CCTK_HOST void track_individual(CCTK_REAL phase_est);
  CCTK_DEVICE CCTK_HOST void track_average(const double prev_phase);

  public:
  CCTK_DEVICE CCTK_HOST bns_tracker(int use_reflevel_, int dist_num_points_);
  CCTK_DEVICE CCTK_HOST void prep(const cGH *cgh);
  CCTK_DEVICE CCTK_HOST const bns_locations& track_stars(const double prev_phase);
  CCTK_DEVICE CCTK_HOST const bns_locations& track_minalp();

};

}


#endif

