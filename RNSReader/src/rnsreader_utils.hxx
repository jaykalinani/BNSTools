#include <loop.hxx>
#include "consts.h"

namespace RNSReader {

using namespace Loop;
using namespace amrex;

class Omega_th_reader{
private:
  double *v_th_arr;
  double *th_arr;
  FILE *in1D;
  int interp_stencil_size;
  // int nrho;

  // Interpolation Helpers
  double *th_sample;
  double *l_i_of_r;
public:
  bool interp_err{false};

  CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  init(FILE *init_in1D, const int stencil_size) {
    // nrho = init_nrho;
    in1D = init_in1D;
    interp_stencil_size = stencil_size;

    // Init Y_e array
    read_1dfile__set_array();

    // Init Interpolation Helpers
    th_sample =
      (double *)The_Managed_Arena()->alloc(stencil_size * sizeof(double));
    l_i_of_r =
      (double *)The_Managed_Arena()->alloc(stencil_size * sizeof(double));

    return;
  }

  CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  read_1dfile__set_array(const int num_header_lines = 0) {
    char *line = NULL;
    char *pEnd;
    size_t len = 0;
    ssize_t read;
    int which_line = 0;
    if (!(v_th_arr =
        (double *)The_Managed_Arena()->alloc((MDIV - 1) * sizeof(double)))) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Cannot allocate memory for v(theta) table");
    }
    // Init Rho Array
    if (!(th_arr =
        (double *)The_Managed_Arena()->alloc((MDIV - 1) * sizeof(double)))) {
        CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                    "Cannot allocate memory for theta table");
    }

    // Skip header line
    while ((read = getline(&line, &len, in1D)) != -1) {
      if (which_line >= num_header_lines) {
        th_arr[which_line] = strtod(line, &pEnd);
        v_th_arr[which_line] = strtod(pEnd, NULL);
      }
      which_line++;
    }
    free(line);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  interpolate_1d_quantity_as_function_of_th(const int numlines_in_file,
                                             const CCTK_REAL th,
                                             CCTK_REAL *restrict f_of_th) {

    // First find the central interpolation stencil index:
    int idx = bisection_idx_finder(th, numlines_in_file, th_arr);
    
#ifdef MAX
#undef MAX
#endif
#define MAX(A, B) (((A) > (B)) ? (A) : (B))

#ifdef MIN
#undef MIN
#endif
#define MIN(A, B) (((A) < (B)) ? (A) : (B))
    
    int idxmin = MAX(0, idx - interp_stencil_size / 2 - 1);
    idxmin = MIN(idxmin, numlines_in_file - interp_stencil_size);
    
    // Now perform the Lagrange polynomial interpolation:
    
    // First set the interpolation coefficients:
    for (int i = idxmin; i < idxmin + interp_stencil_size; i++) {
      th_sample[i - idxmin] = th_arr[i];
    }
    for (int i = 0; i < interp_stencil_size; i++) {
      CCTK_REAL numer = 1.0;
      CCTK_REAL denom = 1.0;
      for (int j = 0; j < i; j++) {
        numer *= th - th_sample[j];
        denom *= th_sample[i] - th_sample[j];
      }
      for (int j = i + 1; j < interp_stencil_size; j++) {
        numer *= th - th_sample[j];
        denom *= th_sample[i] - th_sample[j];
      }
      l_i_of_r[i] = numer / denom;
    }

    // Then perform the interpolation:
    *f_of_th = 0.0;
    for (int i = idxmin; i < idxmin + interp_stencil_size; i++) {
      *f_of_th += l_i_of_r[i - idxmin] * v_th_arr[i];
    }
  }
    // Find interpolation index using Bisection root-finding algorithm:
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline int
  bisection_idx_finder(const CCTK_REAL rrbar, const int numlines_in_file,
                        const CCTK_REAL *restrict rbar_arr) {
    int x1 = 0;
    int x2 = numlines_in_file - 1;
    CCTK_REAL y1 = rrbar - rbar_arr[x1];
    CCTK_REAL y2 = rrbar - rbar_arr[x2];
    if (y1 * y2 >= 0) {
      // Cannot print on GPU
      // fprintf(stderr,"INTERPOLATION BRACKETING ERROR %e | %e
      // %e\n",rrbar,y1,y2); exit(1);

      // Return poison value instead
      interp_err = true;
      return 2555;
    }
    for (int i = 0; i < numlines_in_file; i++) {
      int x_midpoint = (x1 + x2) / 2;
      CCTK_REAL y_midpoint = rrbar - rbar_arr[x_midpoint];
      if (y_midpoint * y1 < 0) {
        x2 = x_midpoint;
        y2 = y_midpoint;
      } else {
        x1 = x_midpoint;
        y1 = y_midpoint;
      }
      if (std::abs(x2 - x1) == 1) {
        // If rbar_arr[x1] is closer to rrbar than rbar_arr[x2] then return x1:
        // if(fabs(rrbar-rbar_arr[x1]) < fabs(rrbar-rbar_arr[x2])) return x1;
        // Otherwiser return x2:
        // return x2;
        // Always return the left value
        return x1;
      }
    }
    // Cannot print on GPU
    // fprintf(stderr,"INTERPOLATION BRACKETING ERROR: DID NOT CONVERGE.\n");
    // exit(1);
    
    // Return poison value instead
    interp_err = true;
    return 2555;
  }
};

template <typename T>
CCTK_DEVICE CCTK_HOST T calc_avg_c2v(const GF3D2<T> &gf, const PointDesc &p) {
  //  constexpr auto DI = PointDesc::DI;
  CCTK_REAL gf_avg = 0.0;
  for (int dk = 0; dk < 2; ++dk) {
    for (int dj = 0; dj < 2; ++dj) {
      for (int di = 0; di < 2; ++di) {
        gf_avg += gf(p.I - p.DI[0] * di - p.DI[1] * dj - p.DI[2] * dk);
      }
    }
  }
  return gf_avg / 8.0;
}

} // namespace RNSReader
