#include "bnstracker.hxx"
#include <cassert>
#include <cmath>
#include <AMReX_GpuComplex.H>

using namespace BNSTrackerX;

CCTK_DEVICE CCTK_HOST bns_tracker::bns_tracker(int use_reflevel_, int dist_num_points_)
: use_reflevel(use_reflevel_), 
  dist_num_points(dist_num_points_)
{}

CCTK_DEVICE CCTK_HOST void bns_tracker::get_dens(const gpos& i, CCTK_REAL& rmd, CCTK_REAL& crmd) //TODO
{
  mats_l glo; //TODO
  glo_pc(i) >> glo; //TODO
  CCTK_REAL vc  = sqrt(det(glo));
  rmd         = rmd_pc(i);
  crmd        = vc * lfac_pc(i) * rmd;
}

CCTK_DEVICE CCTK_HOST void bns_tracker::track_individual(CCTK_REAL phase_est)
{
  amrex::GpuComplex rotate = exp(amrex::GpuComplex(0,-phase_est));
  amrex::GpuComplex avg_pos_1(0), avg_pos_2(0), pos_max_1(0), pos_max_2(0);
  CCTK_REAL wsum_1(0), wsum_2(0), maxrmd_1(-1), maxrmd_2(-1);
 
  for (region::iterator i(r_int());i;++i) { //TODO
    vec_u pos   = coord(i); //TODO
    CCTK_REAL weight, rmd; 
    get_dens(i, rmd, weight);
    amrex::GpuComplex cp(pos(0)-cms_x, pos(1)-cms_y);
    cp         *= rotate;
    if (cp.real() > 0) {
      avg_pos_1 += weight * cp;
      wsum_1    += weight;
      if (rmd > maxrmd_1) {
        maxrmd_1  = rmd;
        pos_max_1 = cp;
      }
    }
    else {
      avg_pos_2 += weight * cp;
      wsum_2    += weight;
      if (rmd > maxrmd_2) {
        maxrmd_2  = rmd;
        pos_max_2 = cp;
      }
    }
  }
  avg_pos_1      = global_sum(avg_pos_1); //TODO: check global sum
  wsum_1         = global_sum(wsum_1); //TODO
  avg_pos_2      = global_sum(avg_pos_2); //TODO
  wsum_2         = global_sum(wsum_2); //TODO
  avg_pos_1      = (avg_pos_1 / wsum_1) / rotate;
  avg_pos_2      = (avg_pos_2 / wsum_2) / rotate;

  CCTK_REAL mcnt_1 = (maxrmd_1 < global_max(maxrmd_1)) ? 0.0 : 1.0;
  pos_max_1      = pos_max_1 * mcnt_1;
  mcnt_1         = global_sum(mcnt_1); //TODO
  pos_max_1      = global_sum(pos_max_1); //TODO
  pos_max_1      = (pos_max_1 / mcnt_1) / rotate;

  CCTK_REAL mcnt_2 = (maxrmd_2 < global_max(maxrmd_2)) ? 0.0 : 1.0;
  pos_max_2      = pos_max_2 * mcnt_2;
  mcnt_2         = global_sum(mcnt_2); //TODO
  pos_max_2      = global_sum(pos_max_2); //TODO
  pos_max_2      = (pos_max_2 / mcnt_2) / rotate;

  star1.set_phase_sep(remove_phase_jump(phase_est, 
                        std::arg(avg_pos_1)), abs(avg_pos_1));
  star2.set_phase_sep(remove_phase_jump(phase_est + M_PI, 
                        std::arg(avg_pos_2)), abs(avg_pos_2));

  star1_maxloc.set_phase_sep(remove_phase_jump(phase_est, 
                               std::arg(pos_max_1)), abs(pos_max_1));
  star2_maxloc.set_phase_sep(remove_phase_jump(phase_est + M_PI, 
                               std::arg(pos_max_2)), abs(pos_max_2));
}

CCTK_DEVICE CCTK_HOST void bns_tracker::track_average(const double prev_phase)
{
  amrex::GpuComplex avg_pos(0), cms_pos(0);
  CCTK_REAL cms_zpos(0);
  CCTK_REAL avg_r(0), wsum(0);
  CCTK_REAL rmd(0), crmd(0);
  
  for (region::iterator i(r_int());i;++i) { //TODO
    vec_u pos         = coord(i); //TODO
    get_dens(i, rmd, crmd);
    amrex::GpuComplex cp(pos(0), pos(1));
    cms_pos           += crmd * cp;
    cms_zpos          += crmd * pos(2);
    wsum              += crmd;
  }
  cms_pos           = global_sum(cms_pos); //TODO
  cms_zpos          = global_sum(cms_zpos); //TODO
  wsum              = global_sum(wsum); //TODO
  cms_pos           /= wsum;
  cms_zpos          /= wsum;

  for (region::iterator i(r_int());i;++i) {
    vec_u pos         = coord(i);
    get_dens(i, rmd, crmd);
    amrex::GpuComplex cp       = amrex::GpuComplex(pos(0), pos(1)) - cms_pos;
    CCTK_REAL d         = abs(cp);
    amrex::GpuComplex cp2      = d > 0 ? (cp*(cp/d)) : 0.0;
    avg_r             += crmd * d;
    avg_pos           += crmd * cp2;
  }
  avg_r             = global_sum(avg_r); //TODO
  avg_pos           = global_sum(avg_pos); //TODO
  avg_r             /= wsum;
  avg_pos           /= wsum;
  CCTK_REAL new_phase = std::arg(avg_pos) / 2.0;
  CCTK_REAL new_sep   = avg_r;
 
  avg.set_phase_sep(remove_phase_jump(prev_phase*2.0, 
                        new_phase*2.0) / 2.0 , new_sep);  

  cms_x             = cms_pos.real();
  cms_y             = cms_pos.imag();
  cms_z             = cms_zpos;
}

CCTK_DEVICE CCTK_HOST const bns_locations& bns_tracker::track_minalp()
{
  CCTK_REAL mx(0),my(0),mz(0);
  CCTK_REAL minalp(0);
  bool first(true);
  for (region::iterator i(r_int());i;++i) { //TODO
    vec_u pos   = coord(i); //TODO
    CCTK_REAL alp = alp_pc(i);
    if (first || (alp < minalp)) {
      first = false;
      minalp  = alp;
      mx = pos(0);
      my = pos(1);
      mz = pos(2);
    }
  }
  
  CCTK_REAL mcnt   = (minalp > global_min(minalp)) ? 0.0 : 1.0;
  mx             = mx * mcnt;
  my             = my * mcnt;
  mz             = mz * mcnt;
  mcnt           = global_sum(mcnt); //TODO
  assert(mcnt>0);
  mx             = global_sum(mx); //TODO
  my             = global_sum(my); //TODO
  mz             = global_sum(mz); //TODO
  
  minalp_x       = mx / mcnt;
  minalp_y       = my / mcnt;
  minalp_z       = mz / mcnt;
  
  return *this;
}

CCTK_DEVICE CCTK_HOST CCTK_REAL bns_tracker::remove_phase_jump(CCTK_REAL old_phase, CCTK_REAL new_phase)
{
  while (new_phase > old_phase + M_PI) new_phase -= 2.0*M_PI;
  while (new_phase < old_phase - M_PI) new_phase += 2.0*M_PI;
  return new_phase;
}

CCTK_DEVICE CCTK_HOST void orb_param::set_phase_sep(CCTK_REAL phi, CCTK_REAL r)
{
  sep   = r;
  phase = phi;
  x     = r * cos(phi);
  y     = r * sin(phi);
}

CCTK_DEVICE CCTK_HOST const bns_locations& bns_tracker::track_stars(const double prev_phase) 
{
  track_average(prev_phase);
  track_individual(avg.phase);

  const double dx     = star1.x - star2.x;
  const double dy     = star1.y - star2.y;
  coord_sep_bc        = sqrt(dx * dx + dy * dy);

  const double dx_md  = star1_maxloc.x - star2_maxloc.x;
  const double dy_md  = star1_maxloc.y - star2_maxloc.y;
  coord_sep_md        = sqrt(dx_md * dx_md + dy_md * dy_md);

  return *this;
}


