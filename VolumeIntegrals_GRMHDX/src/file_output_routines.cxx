//TODO: Fix memory allocation in the output routine

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "cctk.h"
#include "cctk_IOMethods.h"
#include "cctk_Parameters.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"

#include "number_of_reductions.cxx"

/********************************************************************
 ********************    Macro Definitions   ************************
 ********************************************************************/
CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void CREATE_OUTDIR(CCTK_ARGUMENTS,const char *input_dir,char *&actual_dir)
  {
    DECLARE_CCTK_PARAMETERS;
    const char *_dir = input_dir;


    /* check whether "dir" was set; if not default to "IO::out_dir" */
    if (*_dir == 0)
      {
        _dir = out_dir;
      }

    /* omit the directory name if it's the current working dir */
    if (strcmp (_dir, ".") == 0)
      {
        actual_dir = strdup ("");
      }
    else
      {
        actual_dir = (char *)malloc (strlen (_dir) + 2);
        sprintf (actual_dir, "%s/", _dir);
      }

    /* create the directory */
    int izzle = IOUtil_CreateDirectory (cctkGH, actual_dir, 0, 0);
    if (izzle < 0)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "VolumeIntegrals_GRMHDX:file_output_routines.C Problem creating directory '%s' "
                    "for output", actual_dir);
      }
    else if (izzle >= 0 && verbose==1)
      {
        CCTK_VInfo (CCTK_THORNSTRING, "Outputting to directory '%s'",
                    actual_dir);
      }
  }

CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void OUTDIR_NAME(CCTK_ARGUMENTS,const char *input_dir,char *&actual_dir)
{
    DECLARE_CCTK_PARAMETERS;
    const char *_dir = input_dir;


    /* check whether "dir" was set; if not default to "IO::out_dir" */
    if (*_dir == 0)
      {
        _dir = out_dir;
      }

    /* omit the directory name if it's the current working dir */
    if (strcmp (_dir, ".") == 0)
      {
        actual_dir = strdup ("");
      }
    else
      {
        actual_dir = (char *)malloc (strlen (_dir) + 2);
        sprintf (actual_dir, "%s/", _dir);
      }
}

extern "C" void VI_GRMHDX_file_output_routine_Startup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  if(enable_file_output==1) {
    char *actual_dir;
    CREATE_OUTDIR(CCTK_PASS_CTOC,outVolIntegral_dir,actual_dir);
  }
}

extern "C" void VI_GRMHDX_file_output(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_VI_GRMHDX_file_output;
  DECLARE_CCTK_PARAMETERS;
  if(CCTK_MyProc(cctkGH)==0 && cctk_iteration%VolIntegral_out_every==0 && enable_file_output==1) {
    char* actual_dir;
    OUTDIR_NAME(CCTK_PASS_CTOC,outVolIntegral_dir,actual_dir);

    /* Next loop over all integrals, output a file for each */

    /* allocate filename string buffer and build the filename */
    char *filename = (char *)malloc (strlen (actual_dir) +
                                     strlen ("volume_integrals-GRMHDX.asc") + 10); // <- Tiny amount of memory.
    sprintf (filename, "%svolume_integrals-GRMHDX.asc", actual_dir);
    FILE *file = fopen (filename,"a+");
    if (! file) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "VolumeIntegrals_GRMHDX: Cannot open output file '%s'", filename);
      exit(1);
    } else {
      fseek(file, 0, SEEK_END);
      long size = ftell(file);
      /* First print header if file size is zero. */
      if(size==0) {
        char *header_buffer = (char *) malloc(sizeof(char)*500000); // <- Tiny amount of memory.
        sprintf(header_buffer,"# Col. 1: Time (cctk_time)\n");

        int which_col=2;

        if(enable_time_reparameterization) { sprintf(header_buffer,"%s# Col. 2: Physical Time (physical_time)\n",header_buffer); which_col++; }

        for(int which_integral=1;which_integral<=NumIntegrals;which_integral++) {
          if(volintegral_inside_sphere__radius[which_integral]>0 && volintegral_outside_sphere__radius[which_integral]>0) {
            sprintf(header_buffer,
                    "%s# Col. %d: %s. Init. IN sphere @ (%e,%e,%e), r=%e. OUT sphere @ (%e,%e,%e), r=%e. Moves/Tracks AMR Centre %d/%d\n",
                    header_buffer,which_col,Integration_quantity_keyword[which_integral],
                    volintegral_sphere__center_x_initial[which_integral],
                    volintegral_sphere__center_y_initial[which_integral],
                    volintegral_sphere__center_z_initial[which_integral],
                    volintegral_inside_sphere__radius[which_integral],
                    volintegral_outside_sphere__center_x[which_integral],
                    volintegral_outside_sphere__center_y[which_integral],
                    volintegral_outside_sphere__center_z[which_integral],
                    volintegral_outside_sphere__radius[which_integral],
                    amr_centre__tracks__volintegral_inside_sphere[which_integral],
                    volintegral_sphere__tracks__amr_centre[which_integral]);
          } else if(volintegral_inside_sphere__radius[which_integral]>0) {
            sprintf(header_buffer,
                    "%s# Col. %d: %s. Init. IN sphere @ (%e,%e,%e), r=%e. Moves/Tracks AMR Centre %d/%d\n",
                    header_buffer,which_col,Integration_quantity_keyword[which_integral],
                    volintegral_sphere__center_x_initial[which_integral],
                    volintegral_sphere__center_y_initial[which_integral],
                    volintegral_sphere__center_z_initial[which_integral],
                    volintegral_inside_sphere__radius[which_integral],
                    amr_centre__tracks__volintegral_inside_sphere[which_integral],
                    volintegral_sphere__tracks__amr_centre[which_integral]);
          } else if(volintegral_outside_sphere__radius[which_integral]>0) {
            sprintf(header_buffer,
                    "%s# Col. %d: %s. OUT sphere @ (%e,%e,%e), r=%e.\n",
                    header_buffer,which_col,Integration_quantity_keyword[which_integral],
                    volintegral_outside_sphere__center_x[which_integral],
                    volintegral_outside_sphere__center_y[which_integral],
                    volintegral_outside_sphere__center_z[which_integral],
                    volintegral_outside_sphere__radius[which_integral]);
          } else {
            sprintf(header_buffer,
                    "%s# Col. %d: %s. (Integral over full grid)\n",
                    header_buffer,which_col,Integration_quantity_keyword[which_integral]);
          }
            which_col+=VI_GRMHDX_number_of_reductions(which_integral);
        }
        fprintf(file,"%s",header_buffer);
        free(header_buffer);
      }
      /* Next print one line of data. */
      char *buffer = (char *) malloc(sizeof(char)*500000); // <- Tiny amount of memory.
      sprintf(buffer,"%e ",cctk_time);
      if(enable_time_reparameterization) {
        double xx = cctk_time;
        //Mathematica input:  CForm[Integrate[((1/2)*(Erf[((t) - (10))/(5)] + 1)), {t, 0, xx}]]
        // Output: (((-Power(E,-Power(t0,2)/Power(w,2)) + Power(E,-Power(t0 - xx,2)/Power(w,2)))*w)/Sqrt(Pi) + xx - t0*Erf(t0/w) + (-t0 + xx)*Erf((-t0 + xx)/w))/2.
        // Parser: echo "(((-Power(E,-Power(t0,2)/Power(w,2)) + Power(E,-Power(t0 - xx,2)/Power(w,2)))*w)/Sqrt(Pi) + xx - t0*Erf(t0/w) + (-t0 + xx)*Erf((-t0 + xx)/w))/2."|sed "s/Power(E,/exp(/g;s/Power(/pow(/g;s/Sqrt(/sqrt(/g;s/Erf(/erf(/g;s/Pi/M_PI/g"
        // Parser output: (((-exp(-pow(t0,2)/pow(w,2)) + exp(-pow(t0 - xx,2)/pow(w,2)))*w)/sqrt(M_PI) + xx - t0*erf(t0/w) + (-t0 + xx)*erf((-t0 + xx)/w))/2.
        double t0 = VIv_time_reparam_t0;
        double w  = VIv_time_reparam_w;
        *physical_time = (((-exp(-pow(t0,2)/pow(w,2)) + exp(-pow(t0 - xx,2)/pow(w,2)))*w)/sqrt(M_PI) + xx - t0*erf(t0/w) + (-t0 + xx)*erf((-t0 + xx)/w))/2.;
        //0.5 * ( (xx-t0)*erf((xx-t0)/w) - t0*erf(2.0) + (w/sqrt(M_PI))*(exp(-(xx-t0)*(xx-t0)/(w*w)) - 1.0/exp(t0*t0/(w*w))) + xx);


        if(fabs(cctk_time) < 1e-8*cctk_delta_time) *physical_time=0.0;

        sprintf(buffer,"%s%e ",buffer,*physical_time);
      }
      for(int which_integral=1;which_integral<=NumIntegrals;which_integral++) {
        int num_reductions = VI_GRMHDX_number_of_reductions(which_integral);
        for(int reduction=0;reduction<num_reductions;reduction++) {
          if(which_integral==NumIntegrals && reduction==num_reductions-1) {
            sprintf(buffer, "%s%.15e\n",buffer,VolIntegral[4*(which_integral) + reduction]);
          } else {
            sprintf(buffer, "%s%.15e ", buffer,VolIntegral[4*(which_integral) + reduction]);
          }
        }
      }
      fprintf(file,"%s",buffer);
      fclose(file);
      free(buffer);
    }
    free(filename);
  }
}
