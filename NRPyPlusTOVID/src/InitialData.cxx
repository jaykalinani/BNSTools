
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include <loop.hxx>
#include <loop_device.hxx>

#include "aster_utils.hxx"

#include "ADMQuantities.h"
#include "HydroQuantities.h"
#include "interpolate_TOV_solution_to_point.h"
#include "convert_TOV_spacetime_vars_to_ADM_vars.h"

namespace NRPyPlusTOVID {
using namespace std;
using namespace AsterUtils;
using namespace Loop;


// Alias for "vel" vector gridfunction:
// #define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
// #define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
// #define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

void read_TOV_input_data_from_file(ID_inputs *TOV_in) {

    DECLARE_CCTK_PARAMETERS;

    // Step 1: Set up TOV initial data
    // Step 1.a: Read TOV initial data from data file
    // Open the data file:
    char filename[100];
    sprintf(filename,"%s",TOV_filename); // TOV_filename is a CCTK_PARAMETER
    FILE *in1Dpolytrope = fopen(filename, "r");
    if (in1Dpolytrope == NULL) {
        fprintf(stderr,"ERROR: could not open file %s\n",filename);
        exit(1);
    }
    // Count the number of lines in the data file:
    int numlines_in_file = count_num_lines_in_file(in1Dpolytrope);
    // Allocate space for all data arrays:
    CCTK_REAL *r_Schw_arr     = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*numlines_in_file);
    CCTK_REAL *rho_arr        = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*numlines_in_file);
    CCTK_REAL *rho_baryon_arr = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*numlines_in_file);
    CCTK_REAL *P_arr          = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*numlines_in_file);
    CCTK_REAL *M_arr          = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*numlines_in_file);
    CCTK_REAL *expnu_arr      = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*numlines_in_file);
    CCTK_REAL *exp4phi_arr    = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*numlines_in_file);
    CCTK_REAL *rbar_arr       = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*numlines_in_file);

    // Read from the data file, filling in arrays.
    // read_datafile__set_arrays() may be found in TOV/tov_interp.h
    if(read_datafile__set_arrays(in1Dpolytrope, r_Schw_arr,rho_arr,rho_baryon_arr,P_arr,M_arr,expnu_arr,exp4phi_arr,rbar_arr) == 1) {
        fprintf(stderr,"ERROR WHEN READING FILE %s!\n",filename);
        exit(1);
    }
    fclose(in1Dpolytrope);
    REAL Rbar = -100;
    int Rbar_idx = -100;
    for(int i=1;i<numlines_in_file;i++) {
        if(rho_arr[i-1]>0 && rho_arr[i]==0) { Rbar = rbar_arr[i-1]; Rbar_idx = i-1; }
    }
    if(Rbar<0) {
        fprintf(stderr,"Error: could not find rbar=Rbar from data file.\n");
        exit(1);
    }

    TOV_in->Rbar = Rbar;
    TOV_in->Rbar_idx = Rbar_idx;

    const int interp_stencil_size = 12;
    TOV_in->interp_stencil_size = interp_stencil_size;
    TOV_in->numlines_in_file = numlines_in_file;

    TOV_in->r_Schw_arr     = r_Schw_arr;
    TOV_in->rho_arr        = rho_arr;
    TOV_in->rho_baryon_arr = rho_baryon_arr;
    TOV_in->P_arr          = P_arr;
    TOV_in->M_arr          = M_arr;
    TOV_in->expnu_arr      = expnu_arr;
    TOV_in->exp4phi_arr    = exp4phi_arr;
    TOV_in->rbar_arr       = rbar_arr;
    /* END TOV INPUT ROUTINE */
}

extern "C" void NRPyPlusTOVID_ET_InitialData(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_NRPyPlusTOVID_ET_InitialData;
  DECLARE_CCTK_PARAMETERS;

  ID_inputs TOV_in;
  read_TOV_input_data_from_file(&TOV_in);

  const array<CCTK_INT, dim> indextype = {1, 1, 1};
  const GF3D2layout CCC_layout(cctkGH, indextype);

  grid.loop_all<1, 1, 1>(grid.nghostzones,
                         [&](const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_INT idx = CCC_layout.linear(p.I);
        CCTK_REAL rr = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
        CCTK_REAL th = acos(p.z/rr);

        CCTK_REAL IDexp_4phi,IDnu,IDPressure,IDrho_baryonic,IDrho__total_energy_density;
        interpolate_TOV_solution_to_point(rr, TOV_in, &IDexp_4phi,&IDnu,
                                          &IDPressure,&IDrho_baryonic,&IDrho__total_energy_density);

        CCTK_REAL IDalpha,IDgammaDD00,IDgammaDD01,IDgammaDD02,IDgammaDD11,IDgammaDD12,IDgammaDD22;
        convert_TOV_spacetime_vars_to_ADM_vars(rr, th, IDexp_4phi,IDnu,
            &IDalpha,&IDgammaDD00,&IDgammaDD01,&IDgammaDD02,&IDgammaDD11,&IDgammaDD12,&IDgammaDD22);

        HydroQuantities(p,
                        IDPressure,IDrho_baryonic,IDrho__total_energy_density,
                        press,rho,eps,velx,vely,velz);

        ADMQuantities(p,
                      IDalpha,IDgammaDD00,IDgammaDD01,IDgammaDD02,IDgammaDD11,IDgammaDD12,IDgammaDD22,
                      alp_cc,betax_cc,betay_cc,betaz_cc,
                      gxx_cc,gxy_cc,gxz_cc,gyy_cc,gyz_cc,gzz_cc);
  });

  free(TOV_in.r_Schw_arr);
  free(TOV_in.rho_arr);
  free(TOV_in.rho_baryon_arr);
  free(TOV_in.P_arr);
  free(TOV_in.M_arr);
  free(TOV_in.expnu_arr);
  free(TOV_in.exp4phi_arr);
  free(TOV_in.rbar_arr);
}

extern "C" void NRPyPlusTOVID_Interpolation_C2V(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_NRPyPlusTOVID_Interpolation_C2V;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO("Starting interpolation for ADM variables.");
  grid.loop_int<0, 0, 0>(grid.nghostzones,
                         [=] CCTK_HOST(const Loop::PointDesc &p)
                             CCTK_ATTRIBUTE_ALWAYS_INLINE {
                               gxx(p.I) = calc_avg_c2v<4>(gxx_cc, p);
                               gxy(p.I) = calc_avg_c2v<4>(gxy_cc, p);
                               gxz(p.I) = calc_avg_c2v<4>(gxz_cc, p);
                               gyy(p.I) = calc_avg_c2v<4>(gyy_cc, p);
                               gyz(p.I) = calc_avg_c2v<4>(gyz_cc, p);
                               gzz(p.I) = calc_avg_c2v<4>(gzz_cc, p);

                               kxx(p.I) = 0.0;
                               kxy(p.I) = 0.0;
                               kxz(p.I) = 0.0;
                               kyy(p.I) = 0.0;
                               kyz(p.I) = 0.0;
                               kzz(p.I) = 0.0;

                               alp(p.I) = calc_avg_c2v<4>(alp_cc, p);
                               betax(p.I) = calc_avg_c2v<4>(betax_cc, p);
                               betay(p.I) = calc_avg_c2v<4>(betay_cc, p);
                               betaz(p.I) = calc_avg_c2v<4>(betaz_cc, p);

                               dtalp(p.I) = 0.0;
                               dtbetax(p.I) = 0.0;
                               dtbetay(p.I) = 0.0;
                               dtbetaz(p.I) = 0.0;
                             });

  CCTK_INFO("Done interpolation for ADM variables.");
}

extern "C" void NRPyPlusTOVID_ET_VelPert(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_NRPyPlusTOVID_ET_VelPert;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all<1, 1, 1>(grid.nghostzones,
                         [&](const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    if (rho(p.I) > rho_atmosphere) {
      velx(p.I) = pert_amp * p.x;
      vely(p.I) = pert_amp * p.y;
      velz(p.I) = pert_amp * p.z;
    }
    else {
      velx(p.I) = 0.0; 
      vely(p.I) = 0.0; 
      velz(p.I) = 0.0; 
    }
  });
}

} // namespace
