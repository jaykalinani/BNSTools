#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "util_Table.h"
#include "util_String.h"

/* definition of rest mass: Eq. (11) in arXiv:0812.2245v2 M = \int \rho W \sqrt\gamma d^3x
 * this meshes with the EOM for hydro which gives: D_{,t} + (D [\alpha v^i - \beta^i])_{,\i} = 0 
 * (eg. Font lrr-2008-7 Eq. (28) ff) where D=D_{whisky} = \sqrt\gamma*D_{lrr}
 *
 * with this we can easyly get M_{,t} = - \oint_{\partial V} D [\alpha v^i - \beta^i] d\Sigma_i
 * accessing only HydroBase and ADMBase variables this means we need: 
 *   D = \sqrt\gamma W \rho, \gamma = \sqrt{\det g_{ij}}, W^2 = 1/(1-v^2)
 * this calculation is performed in get_gab_ja_onto_detector below, although I think one should 
 * really consider using a generalized get_vars_onto_detector coupled to an expression evaluator 
 * later on (roland)
 */

/******************************
 ***** Scheduled routines ****
 ******************************/
void outflow (CCTK_ARGUMENTS);

/******************************
 ***** Hard coded limits ******
 ******************************/
#define DIM 3
#define MAX_NUMBER_DETECTORS 100
#define MAX_NUMBER_EXTRAS 20
#define MAX_NUMBER_TRESHOLDS 20
#define NUM_INPUT_ARRAYS 6+4+4
#define NUM_OUTPUT_ARRAYS 6+4+4
#define NGHOSTS 2
#define ARRAY_INIT_VALUE 0.

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279
#endif
static inline CCTK_REAL pow2(CCTK_REAL x) {return x*x;}

/******************************************************
 ***** records of which files need to be truncated ****
 ******************************************************/
static CCTK_INT files_created = -1;

/*************************
 ***** Local routines ****
 *************************/

/* parse string of variables into array of variable indices (copied from Multipole) */
static void fill_variable(int idx, const char *optstring, void *callback_arg);

/* compute dr/dtheta and dr/dphi */
static int drdth_drdph(int i, int j,
                int sn,
                CCTK_REAL dth, CCTK_REAL dph,
                CCTK_INT verbose,
                CCTK_INT maxntheta, CCTK_INT maxnphi,
                CCTK_REAL *sf_radius,
                CCTK_REAL *ht, CCTK_REAL *hp);

/* call interpolator and interpolate the current density and Lorentz factor onto
 * the detector surface */
static int get_ja_w_eninf_and_extras_onto_detector(CCTK_ARGUMENTS, CCTK_INT det,
        CCTK_INT num_extras, CCTK_INT extras_ind[MAX_NUMBER_EXTRAS], CCTK_REAL
        *jx, CCTK_REAL *jy, CCTK_REAL *jz, CCTK_REAL *w, CCTK_REAL * eninf,
        CCTK_REAL **extras);

/* utility routine to access an element on the interpolated surface */
static int get_j_w_and_eninf_local(int i, int j, int ntheta, CCTK_REAL
        *j1_det,CCTK_REAL *j2_det,CCTK_REAL *j3_det, CCTK_REAL *w_det,
        CCTK_REAL *eninf_det, CCTK_REAL jloc[3], CCTK_REAL *wloc,
        CCTK_REAL * eninf);
/* utility routine to get value of Cactus parameter surface_projection_<extra_num> */
static CCTK_REAL *get_surface_projection(CCTK_ARGUMENTS, int extra_num);

/* setup memory for interpolator call */
static CCTK_INT outflow_get_local_memory(CCTK_INT npoints);
static CCTK_REAL *outflow_allocate_array(CCTK_INT npoints, const char *name);

/* replace '/' by '\' in file name to have Util_TableSet* accept it as a key */
static char *sanitize_filename(const char *fn);
/* write results to disk */
static int Outflow_write_output(CCTK_ARGUMENTS, CCTK_INT det, CCTK_REAL flux,
        CCTK_REAL avg_w_lorentz, CCTK_REAL avg_eninf,
        const CCTK_REAL *threshold_fluxes);
static int Outflow_write_2d_output(CCTK_ARGUMENTS, const char *varname, CCTK_INT
        det, const CCTK_REAL *data_det, const CCTK_REAL *w_det,
        const CCTK_REAL *eninf_det, const CCTK_REAL *surfaceelement_det,
        int num_extras, const CCTK_INT *extras_ind, CCTK_REAL * const extras[]);


/**********************************************************************/
/*** Implementation                                                 ***/
/**********************************************************************/

/**********************************************************************/
/*** IO                                                             ***/
/**********************************************************************/
static char *sanitize_filename(const char *fn)
{
  char *sfn = Util_Strdup(fn);
  if (sfn == NULL) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Could not allocate memory for sanitized file name");
  }

  for (char *s = sfn ; *s != '\0' ; s++) {
    if(*s == '/')
      *s = '\\';
  }

  return sfn;
}

static int Outflow_write_2d_output(CCTK_ARGUMENTS, const char *varname, CCTK_INT
        det, const CCTK_REAL *data_det, const CCTK_REAL *w_det,
        const CCTK_REAL *eninf_det, const CCTK_REAL *surfaceelement_det,
        int num_extras, const CCTK_INT *extras_ind, CCTK_REAL * const extras[])
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  char const *fmode;
  int ierr;
  CCTK_INT file_created;
  char *filename, *key;
  char format_str_fixed[2048];
  char format_str_extras[128];
  size_t len_written;
  FILE *file;

  if (verbose>3) {
    CCTK_VInfo(CCTK_THORNSTRING, "writing output");
  }

  // check input data
  if (det>=MAX_NUMBER_DETECTORS) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "warn: det=%d, but MAX_NUMBER_DETECTORS=%d, increase",
               (int)det,MAX_NUMBER_DETECTORS);
  }
  
  // filename
  Util_asprintf (&filename, "%s/outflow_surface_det_%d_%s.asc", out_dir, det, varname);
  assert(filename);

  // file mode: append if already written
  file_created = 0;
  key = sanitize_filename(filename);
  ierr = Util_TableGetInt(files_created, &file_created, key);
  if (ierr < 0 && ierr != UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Internal error: could not check if file '%s' has already been created: %d",
               filename, ierr);
  }
  fmode = !file_created && IO_TruncateOutputFiles(cctkGH) ? "w" : "a";

  // open file
  file = fopen (filename, fmode);
  if (!file) {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "write_outflow: Could not open scalar output file '%s'",
                filename);
    free(key);
    free(filename);
    return -1;
  }

  // write header on startup
  if (!file_created) {
    const CCTK_INT sn = sf_IdFromName(surface_index[det], surface_name[det]);
    const CCTK_INT ntheta=sf_ntheta[sn]-2*nghoststheta[sn];
    const CCTK_INT nphi=sf_nphi[sn]-2*nghostsphi[sn];
    assert(sn >= 0);
    fprintf(file,"# 2d Outflow\n");
    fprintf(file,"# detector no.=%d ntheta=%d nphi=%d\n",(int)det,(int)ntheta,(int)nphi);
    fprintf(file,"# gnuplot column index:\n");
    fprintf(file,"# 1:it 2:t 3:x 4:y 5:z 6:%s 7:w_lorentz 8:eninf 9:surface_element", varname);
    for (int i = 0 ; i < num_extras ; i++) {
      fprintf(file, " %d:%s", 10+i, CCTK_VarName(extras_ind[i]));
    }
    fprintf(file, "\n");
  }

  // write data
  len_written = snprintf (format_str_fixed, sizeof(format_str_fixed)/sizeof(format_str_fixed[0]),
           "%%d\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s",
           out_format, out_format, out_format, out_format, out_format,
           out_format, out_format, out_format);
  assert(len_written < sizeof(format_str_fixed)/sizeof(format_str_fixed[0]));
  len_written = snprintf (format_str_extras, sizeof(format_str_extras)/sizeof(format_str_extras[0]),
           "\t%%%s", out_format);
  assert(len_written < sizeof(format_str_extras)/sizeof(format_str_extras[0]));

  const CCTK_INT sn = sf_IdFromName(surface_index[det], surface_name[det]);
  const CCTK_INT ntheta=sf_ntheta[sn]-2*nghoststheta[sn];
  const CCTK_INT imin=nghoststheta[sn], imax=sf_ntheta[sn]-nghoststheta[sn]-1;
  const CCTK_INT jmin=nghostsphi[sn], jmax=sf_nphi[sn]-nghostsphi[sn]-1;
  const CCTK_REAL oth=sf_origin_theta[sn];
  const CCTK_REAL oph=sf_origin_phi[sn];
  const CCTK_REAL dth=sf_delta_theta[sn];
  const CCTK_REAL dph=sf_delta_phi[sn];
  assert(sn >= 0);

  CCTK_REAL th, ph, ct,st, cp,sp,rp;
  CCTK_REAL det_x, det_y, det_z;
  for (int i=imin,n=0;i<=imax;i++,n++) { // theta in [0.5 delta_th, pi-0.5 delta_th]
    th=oth + i * dth;
    ct=cos(th);
    st=sin(th);

    for (int j=jmin,m=0;j<=jmax;j++,m++) { // phi in [0,2pi-delta_phi]
      int ind = i + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
      int ind2d = n + ntheta*m;

      ph=oph + j * dph;
      cp=cos(ph);
      sp=sin(ph);
      if (override_radius[det]) {
          rp = radius[det];
          assert(rp > 0.);
      } else {
        rp=rad_rescale[det]*sf_radius[ind];
      }

      if(output_relative_coordinates) {
        det_x=rp*cp*st;
        det_y=rp*sp*st;
        det_z=rp*ct;
      } else {
        det_x=sf_centroid_x[sn]+rp*cp*st;
        det_y=sf_centroid_y[sn]+rp*sp*st;
        det_z=sf_centroid_z[sn]+rp*ct;
      }

      fprintf(file, format_str_fixed, cctk_iteration, cctk_time,
              det_x,det_y,det_z, data_det[ind2d], w_det[ind2d],
              eninf_det[ind2d], surfaceelement_det[ind2d]);
      for (int e = 0 ; e < num_extras ; e++) {
        fprintf(file, format_str_extras, extras[e][ind2d]);
      }
      fprintf(file, "\n");
    }
    /* repeat first angle to ge a closed surface in gnuplot */
    int ind = i + maxntheta *(jmin+maxnphi*sn); // XXX not sf_ntheta!
    int ind2d = n + ntheta*0;
    ph=oph + jmin * dph;
    cp=cos(ph);
    sp=sin(ph);
    if (override_radius[det]) {
        rp = radius[det];
        assert(rp > 0.);
    } else {
      rp=rad_rescale[det]*sf_radius[ind];
    }
    if(output_relative_coordinates) {
      det_x=rp*cp*st;
      det_y=rp*sp*st;
      det_z=rp*ct;
    } else {
      det_x=sf_centroid_x[sn]+rp*cp*st;
      det_y=sf_centroid_y[sn]+rp*sp*st;
      det_z=sf_centroid_z[sn]+rp*ct;
    }

    fprintf(file, format_str_fixed, cctk_iteration, cctk_time, det_x,det_y,det_z,
            data_det[ind2d], w_det[ind2d], eninf_det[ind2d],
            surfaceelement_det[ind2d]);
    for (int e = 0 ; e < num_extras ; e++) {
      fprintf(file, format_str_extras, extras[e][ind2d]);
    }
    fprintf(file, "\n");

    fprintf(file, "\n"); /* create a grid edge for gnuplot */
  }
  fprintf(file, "\n"); /* create a block for gnuplot */

  fclose(file); 
  Util_TableSetInt(files_created, 1, key);
  free(key);
  free(filename);

  return 1;
}

static int Outflow_write_output(CCTK_ARGUMENTS, CCTK_INT det, CCTK_REAL flux,
        CCTK_REAL avg_w_lorentz, CCTK_REAL avg_eninf,
        const CCTK_REAL *threshold_fluxes)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  char const *fmode;
  int ierr;
  CCTK_INT file_created;
  char *filename, *key;
  char varname[1024];
  char file_extension[5]=".asc";
  char format_str_real[2048];
  int thresh, col;
  FILE *file;

  if (verbose>3) {
    CCTK_VInfo(CCTK_THORNSTRING, "writing output");
  }
  /* Threshold options */
  int threshold_on_w_lorentz = CCTK_Equals(threshold_on_var, "w_lorentz");
  int threshold_on_eninf = CCTK_Equals(threshold_on_var, "eninf");
  assert(threshold_on_w_lorentz || threshold_on_eninf);

  // filename
  sprintf(varname, "outflow_det_%d",(int)det);

  filename = (char *) malloc (strlen (out_dir) + strlen (varname) +
                              strlen (file_extension) +2);
  assert(filename);
  sprintf (filename, "%s/%s%s", out_dir, varname, file_extension);

  // file mode: append if already written
  file_created = 0;
  key = sanitize_filename(filename);
  ierr = Util_TableGetInt(files_created, &file_created, key);
  if (ierr < 0 && ierr != UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Internal error: could not check if file '%s' has already been created: %d",
               filename, ierr);
  }
  fmode = !file_created && IO_TruncateOutputFiles(cctkGH) ? "w" : "a";

  // open file
  file = fopen (filename, fmode);
  if (!file) {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "write_outflow: Could not open scalar output file '%s'",
                filename);
    free(key);
    free(filename);
    return -1;
  }

  // write header on startup
  if (!file_created) {
    fprintf(file,"# Outflow\n");
    fprintf(file,"# detector no.=%d\n",(int)det);
    fprintf(file,"# gnuplot column index:\n");
    fprintf(file,"# 1:it 2:t 3:flux 4:avg(w_lorentz) 5:avg(eninf)");
    if(num_thresholds > 0) {
      col = 6;
      for(thresh = 0 ; thresh < num_thresholds ; thresh++) {
        if(threshold_on_w_lorentz) {
          fprintf(file," %d:w>=%g",col++,threshold[thresh]);
        }
        if(threshold_on_eninf) {
          fprintf(file," %d:eninf>=%g",col++,threshold[thresh]);
        }
      }
    }
    fprintf(file,"\n");
  }

  // write data
  sprintf (format_str_real,
           "%%d\t%%%s\t%%%s\t%%%s\t%%%s",
           out_format,out_format,out_format,out_format);
  fprintf(file, format_str_real, cctk_iteration, cctk_time, flux, avg_w_lorentz,
          avg_eninf);
  sprintf (format_str_real, "\t%%%s", out_format);
  for(thresh = 0 ; thresh < num_thresholds ; thresh++) {
    fprintf(file,format_str_real,threshold_fluxes[thresh]);
  }
  fprintf(file,"\n");

  fclose(file); 
  Util_TableSetInt(files_created, 1, key);
  free(key);
  free(filename);

  return 1;
}

/**********************************************************************/
/* Utility                                                            */
/**********************************************************************/

/* return pointer to surface_projection_<extra_num> */
static CCTK_REAL *get_surface_projection(CCTK_ARGUMENTS, int extra_num)
{
  DECLARE_CCTK_ARGUMENTS;

  CCTK_REAL *retval = NULL;
  switch(extra_num)
  {
    case  0: retval = surface_projection_0; break;
    case  1: retval = surface_projection_1; break;
    case  2: retval = surface_projection_2; break;
    case  3: retval = surface_projection_3; break;
    case  4: retval = surface_projection_4; break;
    case  5: retval = surface_projection_5; break;
    case  6: retval = surface_projection_6; break;
    case  7: retval = surface_projection_7; break;
    case  8: retval = surface_projection_8; break;
    case  9: retval = surface_projection_9; break;
    case 10: retval = surface_projection_10; break;
    case 11: retval = surface_projection_11; break;
    case 12: retval = surface_projection_12; break;
    case 13: retval = surface_projection_13; break;
    case 14: retval = surface_projection_14; break;
    case 15: retval = surface_projection_15; break;
    case 16: retval = surface_projection_16; break;
    case 17: retval = surface_projection_17; break;
    case 18: retval = surface_projection_18; break;
    case 19: retval = surface_projection_19; break;
    default: CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "invalid extra variable number %d passed. Must be less than 20",
                extra_num);
             break;
  }

  return retval;
}

/**********************************************************************/
/* Interpolator call and surface access                               */
/**********************************************************************/

/* fills j1...j3,w,eninf and the extras with the interpolated numbers */
static int get_ja_w_eninf_and_extras_onto_detector(CCTK_ARGUMENTS, CCTK_INT det,
        CCTK_INT num_extras, CCTK_INT extras_ind[MAX_NUMBER_EXTRAS], CCTK_REAL
        *jx, CCTK_REAL *jy, CCTK_REAL *jz, CCTK_REAL *w, CCTK_REAL *eninf,
        CCTK_REAL **extras)
{
  DECLARE_CCTK_ARGUMENTS; 
  DECLARE_CCTK_PARAMETERS; 
  int ierr;
  CCTK_INT sn,ind,ind2d;
  CCTK_REAL th,ph,ct,st,cp,sp,rp;
  CCTK_INT ntheta,nphi,npoints;
  // auxilliary variables used in constructing j
  static CCTK_REAL *rho0 = NULL, *velx = NULL, *vely = NULL, *velz = NULL;
  static CCTK_REAL *beta1 = NULL, *beta2 = NULL, *beta3 = NULL, *alpha = NULL;
  static CCTK_REAL *g11 = NULL, *g12 = NULL, *g13 = NULL, *g22 = NULL;
  static CCTK_REAL *g23 =  NULL, *g33 = NULL;

  assert(det>=0);
  assert(det<num_detectors);
  assert(jx); assert(jy); assert(jz);
  assert(w); assert(eninf); assert(extras);
  for (int i=0; i<num_extras; i++)
      assert(extras[i]);

  sn = sf_IdFromName(surface_index[det], surface_name[det]);
  assert(sn>=0);

  ntheta=sf_ntheta[sn]-2*nghoststheta[sn];
  nphi=sf_nphi[sn]-2*nghostsphi[sn];
  npoints=ntheta*nphi;

  if (verbose>1) {
    CCTK_VInfo(CCTK_THORNSTRING,"surface %d (%g,%g,%g) nth,nph (%d,%d)",
	       (int)sn,sf_centroid_x[sn],sf_centroid_y[sn],sf_centroid_z[sn],
	       (int)ntheta,(int)nphi);
  }

  CCTK_INT interp_npoints=npoints;
  // uni-processor code - only work on CPU 0
  const CCTK_INT myproc= CCTK_MyProc(cctkGH);
  if ( myproc != 0 ) {
    interp_npoints=0;
  }

  // allocate memory for auxilliary arrays (of maximum possible size)
# define ALLOCATE_TEMP(name) \
  if(name == NULL) \
    name = outflow_allocate_array(maxntheta*maxnphi, #name); \
  assert(name)
  
  ALLOCATE_TEMP(g11);
  ALLOCATE_TEMP(g12);
  ALLOCATE_TEMP(g13);
  ALLOCATE_TEMP(g22);
  ALLOCATE_TEMP(g23);
  ALLOCATE_TEMP(g33);
  ALLOCATE_TEMP(rho0); 
  ALLOCATE_TEMP(velx);
  ALLOCATE_TEMP(vely);
  ALLOCATE_TEMP(velz);
  ALLOCATE_TEMP(beta1);
  ALLOCATE_TEMP(beta2);
  ALLOCATE_TEMP(beta3);
  ALLOCATE_TEMP(alpha);
# undef ALLOCATE_TEMP

  // coordinates setup
  const CCTK_INT imin=nghoststheta[sn], imax=sf_ntheta[sn]-nghoststheta[sn]-1;
  const CCTK_INT jmin=nghostsphi[sn], jmax=sf_nphi[sn]-nghostsphi[sn]-1;
  const CCTK_REAL oth=sf_origin_theta[sn];
  const CCTK_REAL oph=sf_origin_phi[sn];
  const CCTK_REAL dth=sf_delta_theta[sn];
  const CCTK_REAL dph=sf_delta_phi[sn];

  CCTK_REAL det_x[npoints], det_y[npoints], det_z[npoints];
  for (int i=imin,n=0;i<=imax;i++,n++) { // theta in [0.5 delta_th, pi-0.5 delta_th]
    th=oth + i * dth;
    ct=cos(th);
    st=sin(th);
    for (int j=jmin,m=0;j<=jmax;j++,m++) { // phi in [0,2pi-delta_phi]
      ind=i + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
      ph=oph + j * dph;
      cp=cos(ph);
      sp=sin(ph);
      ind2d=n + ntheta*m;
      if (override_radius[det]) {
        rp = radius[det];
        assert(rp > 0.);
      } else {
        rp = rad_rescale[det]*sf_radius[ind];
      }
      det_x[ind2d]=sf_centroid_x[sn]+rp*cp*st;
      det_y[ind2d]=sf_centroid_y[sn]+rp*sp*st;
      det_z[ind2d]=sf_centroid_z[sn]+rp*ct;
    }
  }

  const void* interp_coords[3] 
    = { (const void *) det_x,
        (const void *) det_y,
        (const void *) det_z };

  // 3d input arrays
  CCTK_STRING input_array_names[NUM_INPUT_ARRAYS]
    = { "ADMBase::gxx",
        "ADMBase::gxy",
        "ADMBase::gxz",
        "ADMBase::gyy",
        "ADMBase::gyz",
        "ADMBase::gzz",

        "HydroBase::vel[0]",
        "HydroBase::vel[1]",
        "HydroBase::vel[2]",
        "HydroBase::rho",

        "ADMBase::betax",
        "ADMBase::betay",
        "ADMBase::betaz",
        "ADMBase::alp",
      };
  CCTK_INT input_array_indices[NUM_INPUT_ARRAYS + MAX_NUMBER_EXTRAS];
  for(int i = 0 ; i < NUM_INPUT_ARRAYS ; i++) {
    input_array_indices[i] = CCTK_VarIndex(input_array_names[i]);
    if(input_array_indices[i] < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "couldn't find variable '%s'",
        input_array_names[i]);
        return -1; /*NOTREACHED*/
    }
  }
  for(int i = 0 ; i < num_extras ; i++) {
     input_array_indices[NUM_INPUT_ARRAYS + i] = extras_ind[i];
  }

  CCTK_INT output_array_types[NUM_OUTPUT_ARRAYS + MAX_NUMBER_EXTRAS];
  for(int i = 0 ; i < NUM_OUTPUT_ARRAYS + num_extras ; i++) {
    output_array_types[i] = CCTK_VARIABLE_REAL;
  }

  // 2d output arrays
  void * output_arrays[NUM_OUTPUT_ARRAYS + MAX_NUMBER_EXTRAS]
    = { (void *) g11, 
        (void *) g12,
        (void *) g13,
        (void *) g22,
        (void *) g23,
        (void *) g33,

        (void *) velx, 
        (void *) vely,
        (void *) velz,
        (void *) rho0, 

        (void *) beta1, 
        (void *) beta2,
        (void *) beta3,
        (void *) alpha, 
      };
  for(int i = 0 ; i < num_extras ; i++) {
     output_arrays[NUM_OUTPUT_ARRAYS + i] = extras[i];
  }

  CCTK_INT operand_indices[NUM_OUTPUT_ARRAYS + MAX_NUMBER_EXTRAS];
  for(int i = 0 ; i < NUM_OUTPUT_ARRAYS + num_extras  ; i++) {
    operand_indices[i] = i;
  }

  CCTK_INT operation_codes[NUM_OUTPUT_ARRAYS + MAX_NUMBER_EXTRAS];
  for(int i = 0 ; i < NUM_OUTPUT_ARRAYS + num_extras  ; i++) {
    operation_codes[i] = 0;
  }

  // handles setup
  const int operator_handle = CCTK_InterpHandle(interpolator_name);
  if (operator_handle < 0)
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "couldn't find interpolator \"%s\"!",
               interpolator_name);

  int param_table_handle = Util_TableCreateFromString(interpolator_pars);
  if (param_table_handle < 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "bad interpolator parameter(s) \"%s\"!",
               interpolator_pars);
  }
  
  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS + num_extras,
                        operand_indices, "operand_indices");
  
  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS + num_extras, 
                        operation_codes, "operation_codes");
  
  const int coord_system_handle = CCTK_CoordSystemHandle(coord_system);
  if (coord_system_handle < 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
        "can't get coordinate system handle for coordinate system \"%s\"!",
               coord_system);
  }

  // actual interpolation call
  ierr = CCTK_InterpGridArrays(cctkGH,
                               DIM, // number of dimensions 
                               operator_handle,
                               param_table_handle,
                               coord_system_handle,
                               interp_npoints,
                               CCTK_VARIABLE_REAL,
                               interp_coords,
                               NUM_INPUT_ARRAYS + num_extras, // Number of input arrays
                               input_array_indices,
                               NUM_OUTPUT_ARRAYS + num_extras, // Number of output arrays
                               output_array_types,
                               output_arrays);

  if (ierr<0) {
    CCTK_WARN(1,"interpolation screwed up");
    Util_TableDestroy(param_table_handle);
    return -1;
  }

  ierr = Util_TableDestroy(param_table_handle);
  if (ierr != 0) {
    CCTK_WARN(1,"Could not destroy table");
    return -1;
  }

  // compute current from primitive values
  for(int i = 0 ; i < interp_npoints ; i++) {
    CCTK_REAL detg, dens, vlowx, vlowy, vlowz, v2, my_w_lorentz, my_eninf;

    detg = 2*g12[i]*g13[i]*g23[i] + g33[i]*(g11[i]*g22[i] - pow2(g12[i])) -
        g22[i]*pow2(g13[i]) - g11[i]*pow2(g23[i]);
    if( detg < 0. ) 
    {
        static CCTK_INT last_warned = -1;

        if(verbose > 1 || (verbose > 0 && last_warned != cctk_iteration))
        {
          CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "%s: Metric determinant in iteration %6d at %15.6g,%15.6g,%15.6g is :" 
                "%15.6g from data g = [%15.6g,%15.6g,%15.6g,%15.6g,%15.6g,%15.6g]",
                __func__, cctk_iteration, det_x[i],det_y[i],det_z[i], detg, 
                g11[i],g12[i],g13[i],g22[i],g23[i],g33[i]);
          last_warned = cctk_iteration;
        }

        detg = 1.;
    }

    vlowx = g11[i]*velx[i] + g12[i]*vely[i] + g13[i]*velz[i];
    vlowy = g12[i]*velx[i] + g22[i]*vely[i] + g23[i]*velz[i];
    vlowz = g13[i]*velx[i] + g23[i]*vely[i] + g33[i]*velz[i];

    v2 = vlowx*velx[i] + vlowy*vely[i] + vlowz*velz[i];

    my_w_lorentz = sqrt(1. / (1. - v2));
    if( my_w_lorentz < 1. || v2 > 1 ) 
    {
        static CCTK_INT last_warned = -1;

        if(verbose > 1 || (verbose > 0 && last_warned != cctk_iteration))
        {
          CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                "%s: Unphysical Lorentz factor %15.6g, v2 = %15.6g for data "
                "g = [%15.6g,%15.6g,%15.6g,%15.6g,%15.6g,%15.6g] "
                "vel = [%15.6g,%15.6g,%15.6g] occured in iteration %d at location [%15.6g,%15.6g,%15.6g]",
                __func__, my_w_lorentz,v2, g11[i],g12[i],g13[i],g22[i],g23[i],g33[i],
                velx[i],vely[i],velz[i], cctk_iteration,
                det_x[i],det_y[i],det_z[i]);
          last_warned = cctk_iteration;
        }

        my_w_lorentz = 1.;
    }
    /* - u_t - 1 */
    my_eninf = - my_w_lorentz*(vlowx*beta1[i] + vlowy*beta2[i] + vlowz*beta3[i]
            - alpha[i]) - 1.0;
    dens = sqrt(detg)*rho0[i]*my_w_lorentz;

    jx[i] = dens * (alpha[i]*velx[i] - beta1[i]);
    jy[i] = dens * (alpha[i]*vely[i] - beta2[i]);
    jz[i] = dens * (alpha[i]*velz[i] - beta3[i]);
    w[i]  = my_w_lorentz;
    eninf[i] = my_eninf;
  }

  return interp_npoints;

}

static int get_j_w_and_eninf_local(int i, int j, int ntheta,
                  CCTK_REAL *j1_det,CCTK_REAL *j2_det,CCTK_REAL *j3_det,
                  CCTK_REAL *w_det, CCTK_REAL *eninf_det,
                  CCTK_REAL jloc[3], CCTK_REAL *wloc, CCTK_REAL *eninfloc)
{
  CCTK_INT ind2d=i + ntheta*j;
  /* jloc_i - upstairs index */
  jloc[0]=j1_det[ind2d];
  jloc[1]=j2_det[ind2d];
  jloc[2]=j3_det[ind2d];
  *wloc  =w_det[ind2d];
  *eninfloc=eninf_det[ind2d];

  return 1;
}

static int drdth_drdph(int i, int j, 
                int sn,
                CCTK_REAL dth, CCTK_REAL dph,
                CCTK_INT verbose,
                CCTK_INT maxntheta, CCTK_INT maxnphi,
                CCTK_REAL *sf_radius,
                CCTK_REAL *ht, CCTK_REAL *hp)
{
  CCTK_INT ind;
  CCTK_REAL htp1,htm1,hpp1,hpm1;
  CCTK_REAL htp2,htm2,hpp2,hpm2;

  /* dr/dth dr/dph */
  ind=(i+1) + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
  htp1=sf_radius[ind];
  ind=(i-1) + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
  htm1=sf_radius[ind];
  ind=(i+2) + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
  htp2=sf_radius[ind];
  ind=(i-2) + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
  htm2=sf_radius[ind];
  *ht = (1./12.*htm2-2./3.*htm1+2./3.*htp1-1./12.*htp2)/dth;
  if (verbose>5) {
    fprintf(stderr,"  normal : i=%d j=%d ht=%g\n",i,j,*ht);
  }

  ind=i + maxntheta *((j+1)+maxnphi*sn); // XXX not sf_ntheta!
  hpp1=sf_radius[ind];
  ind=i + maxntheta *((j-1)+maxnphi*sn); // XXX not sf_ntheta!
  hpm1=sf_radius[ind];
  ind=i + maxntheta *((j+2)+maxnphi*sn); // XXX not sf_ntheta!
  hpp2=sf_radius[ind];
  ind=i + maxntheta *((j-2)+maxnphi*sn); // XXX not sf_ntheta!
  hpm2=sf_radius[ind];
  *hp = (1./12.*hpm2-2./3.*hpm1+2./3.*hpp1-1./12.*hpp2)/dph;
  if (verbose>5) { 
    fprintf(stderr,"  normal : i=%d j=%d hp=%g\n",i,j,*hp);            
  }

  return 1;
}


static CCTK_REAL *outflow_allocate_array(CCTK_INT npoints, const char *name)
{
  DECLARE_CCTK_PARAMETERS;

  if (npoints<=0) {
    CCTK_WARN(0,"can't allocate array with npoints <=0");
  }
  if (name==NULL) {
    CCTK_WARN(0,"give a name");
  }

  CCTK_REAL *res;
  res=(CCTK_REAL *) malloc(sizeof(CCTK_REAL)*npoints);
  if (res==NULL) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "%s allocation for npoints=%d failed",name,(int)npoints);
  }
  for (int i=0;i<npoints;i++) res[i]=0;
  const CCTK_REAL mbyt=npoints*sizeof(CCTK_REAL)/(1024.0*1024.0);
  if (verbose>0) {
    CCTK_VInfo(CCTK_THORNSTRING,"allocated array %s with %d elements -> %g MB",name,(int)npoints,mbyt);
  }
  return res;
}


static CCTK_REAL *j1_det, *j2_det, *j3_det, *w_det, *eninf_det, *fluxdens_det, *surfaceelement_det;
static CCTK_INT outflow_get_local_memory(CCTK_INT npoints)
{
  DECLARE_CCTK_PARAMETERS;
  
  static CCTK_INT have_integrand_memory=0;

  if (verbose>1) CCTK_INFO("in allocate_memory");

  if (have_integrand_memory==0) {
    if (verbose>0) CCTK_INFO("allocating new memory");
    // current density on detector (vector)
    j1_det=outflow_allocate_array(npoints,"j1_det");
    j2_det=outflow_allocate_array(npoints,"j2_det");
    j3_det=outflow_allocate_array(npoints,"j3_det");
    w_det =outflow_allocate_array(npoints,"w_det");
    eninf_det=outflow_allocate_array(npoints,"eninf_det");
    fluxdens_det =outflow_allocate_array(npoints,"fluxdens_det");
    surfaceelement_det =outflow_allocate_array(npoints,"surfaceelement_det");
    // update memory allocation flag
    have_integrand_memory=1;
  }
  else {
    if (verbose>1) CCTK_INFO("already allocated memory");
    return 2;
  }

  return 1;
}

/* callback routine to add one extra variable to be output */
static void fill_variable(int idx, const char *optstring, void *callback_arg)
{
  assert(idx >= 0);
  assert(callback_arg);

  CCTK_INT *extras_ind = (CCTK_INT * ) callback_arg;

  if(extras_ind[MAX_NUMBER_EXTRAS-1] != -1) { /* no more free slots */
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "too many extra variables, ignoring variable '%s'.",
               CCTK_VarName(idx));
    return;
  }

  /* find the first free slot in extras_ind and store the new index in it */
  for(int i = 0 ; i < MAX_NUMBER_EXTRAS ; i++)
  {
    if(extras_ind[i] == -1) {
      extras_ind[i] = idx;
      break;
    }
  }
}

/**********************************************************************/
/* Scheduled routines                                                 */
/**********************************************************************/

void outflow (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_INT sn, ind, ind2d, ierr;
  CCTK_REAL ht,hp;
  CCTK_REAL sint,sinp,cost,cosp,rp;
  CCTK_INT interp_npoints;
  // variables related to the extra projected grid functions
  CCTK_INT num_extras;
  CCTK_INT extras_ind[MAX_NUMBER_EXTRAS];
  CCTK_REAL *extras[MAX_NUMBER_EXTRAS];

  /* Threshold options */
  int threshold_on_w_lorentz = CCTK_Equals(threshold_on_var, "w_lorentz");
  int threshold_on_eninf = CCTK_Equals(threshold_on_var, "eninf");
  assert(threshold_on_w_lorentz || threshold_on_eninf);

  /* local memory allocation */
  CCTK_INT maxnpoints=maxntheta*maxnphi;
  ierr=outflow_get_local_memory(maxnpoints);
  if (ierr<0) {
    CCTK_WARN(1,"failed to allocate memory");
    return;
  }

  /* clear the grid arrays */
  for(int i = 0 ; i < maxntheta*maxnphi*num_detectors ; i++)
    fluxdens_projected[i] = w_lorentz_projected[i] = eninf_projected[i] =
        ARRAY_INIT_VALUE;
  for(int e = 0 ; e < MAX_NUMBER_EXTRAS ; e++)
  {
    CCTK_REAL *surface_projection = get_surface_projection(CCTK_PASS_CTOC, e);
    for(int i = 0 ; i < maxntheta*maxnphi*num_detectors ; i++)
    {
      surface_projection[i] = ARRAY_INIT_VALUE;
    }
  }

  /* parse variables string and allocate temporary memory for them */
  if(!CCTK_Equals(extra_variables, "")) {
    for(int i = 0 ; i < MAX_NUMBER_EXTRAS ; i++) /* initialize so that we can count later */
    {
      extras_ind[i] = -1;
    }

    ierr = CCTK_TraverseString(extra_variables, fill_variable, extras_ind,
            CCTK_GROUP_OR_VAR);
    assert(ierr > 0);
  
    for(num_extras = 0 ; num_extras < MAX_NUMBER_EXTRAS ; num_extras++) /* count valid indices */
    {
      if(extras_ind[num_extras] == -1)
        break;
      extras[num_extras] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*maxnpoints);
      assert(extras[num_extras]);
    }
  } else {
    num_extras = 0;
  }

  /* loop over detectors */
  for (int det=0;det<num_detectors;det++)
  {
    CCTK_INT my_compute_every;

    /* check parameters and decide if we have to do anythin */
    if ( compute_every_det[det] >= 0 ) {
        my_compute_every = compute_every_det[det];
    } else {
        my_compute_every = compute_every;
    }
    if ( my_compute_every == 0 || cctk_iteration % my_compute_every != 0 ) {
      continue;
    }

    sn = sf_IdFromName(surface_index[det], surface_name[det]);
    if (sn < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "surface number sn=%d is invalid for detector %d", (int)sn,(int)det);
      continue;
    } else if (sn>=sphericalsurfaces_nsurfaces) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "surface number sn=%d too large, increase SphericalSurface::nsurfaces from its current value %d",
                 (int)sn,(int)sphericalsurfaces_nsurfaces);
      continue;
    }
    if (sf_valid[sn]<=0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "didn't find valid detector surface for sn=%d, det=%d",(int)sn,(int)det);
      continue;
    }

    if(nghoststheta[sn] < NGHOSTS || nghostsphi[sn] < NGHOSTS) { // we need at least NGHOSTS ghost zones 
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "number of ghost zones for spherical surface %d must be at least %d.",(int)sn,NGHOSTS);
      continue;
    }

    /* interpolate onto surface */
    const CCTK_INT imin=nghoststheta[sn], imax=sf_ntheta[sn]-nghoststheta[sn]-1;
    const CCTK_INT jmin=nghostsphi[sn], jmax=sf_nphi[sn]-nghostsphi[sn]-1;
    const CCTK_INT ntheta = imax-imin+1, nphi = jmax-jmin+1;
    const CCTK_REAL oth=sf_origin_theta[sn];
    const CCTK_REAL oph=sf_origin_phi[sn];
    const CCTK_REAL dth=sf_delta_theta[sn];
    const CCTK_REAL dph=sf_delta_phi[sn];
    const CCTK_REAL dtp=dth*dph;

    if (verbose>2) {
      CCTK_VInfo(CCTK_THORNSTRING,"ntheta=%d nphi=%d dth=%g dph=%g",
                 (int)ntheta,(int)nphi,dth,dph);
    }

    interp_npoints=get_ja_w_eninf_and_extras_onto_detector(CCTK_PASS_CTOC, det,
                    num_extras, extras_ind, j1_det, j2_det, j3_det, w_det,
                    eninf_det, extras);
    if (interp_npoints<0) {
      CCTK_WARN(1,"unable to get g_ab, j^a and the extra variables onto the detector. not doing anything.");
      continue;
    }
    if (interp_npoints==0) {
      /* nothing to do (we are not cpu 0) */
      if (verbose > 1) {
        CCTK_VInfo(CCTK_THORNSTRING, "I have nothing to do for detector %d", det);
      }
      continue;
    }

    /*****************************************/
    /* all code below executes only on cpu 0 */
    /*****************************************/

    if (verbose > 1) {
      CCTK_VInfo(CCTK_THORNSTRING, "integrating detector %d", det);
    }

    /* integrate on sphere */
    CCTK_REAL rdn[3], rhat[3], phihat[3], thetahat[3];
    CCTK_REAL jloc[3], wloc, eninfloc;
    CCTK_REAL th,ph;
    CCTK_REAL sum, sum_thresh[MAX_NUMBER_TRESHOLDS], sum_w_lorentz, area; // the value of the flux integral
    CCTK_REAL sum_eninf;

    CCTK_REAL iwtheta,iwphi,intweight;
    /* init integration vars */
    sum = sum_w_lorentz = sum_eninf = area = 0.;
    for (int t = 0 ; t < MAX_NUMBER_TRESHOLDS ; t++) {
      sum_thresh[t] = 0.;
    }

    /* loop over detector surface */
    for (int i=imin,n=0;i<=imax;i++,n++) // theta in [0.5 delta_th, pi-0.5 delta_th]
    {
      th=oth + i * dth;
      cost=cos(th);
      sint=sin(th);
      // weigths from NR(C++,2) 4.1.14 plus extrapolated 1/2 integral (see Maple worksheet)
      if      (i==imin+0 || i==imax-0) iwtheta=13.0/12.0;
      else if (i==imin+1 || i==imax-1) iwtheta= 7.0/ 8.0;
      else if (i==imin+2 || i==imax-2) iwtheta=25.0/24.0;
      else iwtheta=1;

      for (int j=jmin,m=0;j<=jmax;j++,m++) // phi in [0,2pi-delta_phi]
      {
        ph=oph + j * dph;
        cosp=cos(ph);
        sinp=sin(ph);
	iwphi=1; // trapezoid rule
	intweight=iwphi*iwtheta;

        ind=i + maxntheta *(j+maxnphi*sn); // XXX not sf_ntheta!
        if (override_radius[det]) {
            rp = radius[det];
            assert(rp > 0.);
        } else {
          rp=rad_rescale[det]*sf_radius[ind];
        }

        if (verbose>5) {
          fprintf(stderr,"r=%g theta=%g phi=%g\n",rp,th,ph);
        }

        // this computes Eq. (2), (7) of the documentation

        // operates on interpolated values
        get_j_w_and_eninf_local(n,m,ntheta,
                      j1_det,j2_det,j3_det,
                      w_det,eninf_det,jloc,
                      &wloc,&eninfloc);

        // the flat space-like 3d unit vectors
        rhat    [0] =  cosp*sint;rhat    [1] =  sinp*sint;rhat    [2] =  cost;
        thetahat[0] =  cosp*cost;thetahat[1] =  sinp*cost;thetahat[2] = -sint;
        phihat  [0] = -sinp     ;phihat  [1] =  cosp     ;phihat  [2] =     0;

        /* get derivatives of r in theta and phi direction */
        if (override_radius[det]) {
          ht = hp = 0.; /* spherical */
        } else {
          ierr=drdth_drdph(i, j, sn, dth,dph, verbose, maxntheta, maxnphi,
                  sf_radius, &ht, &hp);
	  ht=rad_rescale[det]*ht;
	  hp=rad_rescale[det]*hp;
          if (ierr<0) {
            CCTK_WARN(1,"derivative computation failed");
            continue;
          }
        }

        // the vector surface element
        CCTK_REAL mag_rdn = 0;
        for(int idir = 0 ; idir < 3 ; idir++)
        {
          rdn[idir] = pow2(rp)*sint*rhat[idir] - ht*rp*sint*thetahat[idir] -
                      hp*rp*phihat[idir];
          mag_rdn += pow2(rdn[idir]);
        }
        mag_rdn = sqrt(mag_rdn);

        // sum the integral
        CCTK_REAL df = 0, fluxdens_temp = 0;
        for (int a=0;a<3;a++) {
          df += jloc[a] * rdn[a] * intweight * dtp;
          fluxdens_temp += jloc[a] * rdn[a]/mag_rdn;
        }
        fluxdens_det[n + ntheta*m] = fluxdens_temp; /* store flux density for output */
        surfaceelement_det[n + ntheta*m] = mag_rdn * dtp; /* area element */

        // flux
        sum += df;
        if (verbose>4) {
          fprintf(stderr,"sum=%g\n",sum);
        }

        // Lorentz factor
        sum_w_lorentz += wloc * intweight * mag_rdn * dtp;

        // Specific energy at infinity
        sum_eninf += eninfloc * intweight * mag_rdn * dtp;

        // area of detector
        area += intweight * mag_rdn * dtp;

        for(int t = 0 ; t < num_thresholds ; t++)
        {
          if(threshold_on_w_lorentz) {
            if(wloc >= threshold[t]) {
              sum_thresh[t] += df;
            }
          }
          else if(threshold_on_eninf) {
            if(eninfloc >= threshold[t]) {
              sum_thresh[t] += df;
            }
          }
          if (verbose>4) {
            fprintf(stderr,"sum_thresh[%d]=%g\n",t,sum_thresh[t]);
          }
        }

      } // j : phi
    } // i : theta
    sum_w_lorentz /= area; // average w_lorentz
    sum_eninf /= area; // average eninf

    if (verbose>0) {
      CCTK_VInfo(CCTK_THORNSTRING,"flux value=%g on detector %d", sum,det);
    }

    // main output
    outflow_flux[det]=sum;

    /* store results in grid arrays, translating indices on the way */
    /* we fill in only the upper left corner of the grid array */
    for (int i=imin,n=0;i<=imax;i++,n++) // theta in [0.5 delta_th, pi-0.5 delta_th]
    {
      for (int j=jmin,m=0;j<=jmax;j++,m++) // phi in [0,2pi-delta_phi]
      {
        ind = i + maxntheta * (j + maxnphi*det);
        ind2d = n + ntheta * m;
        fluxdens_projected[ind] = fluxdens_det[ind2d];
        w_lorentz_projected[ind] = w_det[ind2d];
        eninf_projected[ind] = eninf_det[ind2d];
      }
    }
    for(int e = 0 ; e < num_extras ; e++)
    {
      CCTK_REAL *surface_projection = get_surface_projection(CCTK_PASS_CTOC, e);
      for (int i=imin,n=0;i<=imax;i++,n++) // theta in [0.5 delta_th, pi-0.5 delta_th]
      {
        for (int j=jmin,m=0;j<=jmax;j++,m++) // phi in [0,2pi-delta_phi]
        {
          ind = i + maxntheta * (j + maxnphi*det);
          ind2d = n + ntheta * m;
          surface_projection[ind] = extras[e][ind2d];
        }
      }
    }

    /* IO (we only get here if we are CPU #0) */
    ierr=Outflow_write_output(CCTK_PASS_CTOC,det, sum, sum_w_lorentz, sum_eninf, sum_thresh);
    if (ierr<0) {
      CCTK_WARN(1,"writing of information to files failed");
    }
    if (output_2d_data) {
      ierr=Outflow_write_2d_output(CCTK_PASS_CTOC, "fluxdens", det,
              fluxdens_det, w_det, eninf_det, surfaceelement_det,
              num_extras, extras_ind, extras);
      if (ierr<0) {
        CCTK_WARN(1,"writing of fluxdens information to files failed");
      }
    }
  } // det loop over detector number

  /* free temporary memory */
  for(int i = 0 ; i < num_extras ; i++)
  {
    free(extras[i]);
  }
}

int outflow_setup (void)
{
  files_created = Util_TableCreate(0);
  if (files_created < 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Internal error: could not create files_created table: %d",
               (int)files_created);
  }

  return 0;
}
