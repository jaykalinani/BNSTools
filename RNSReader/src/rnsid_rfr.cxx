/* Changed by N. Stergioulas, 24/10/2001: include latest source files from
 * cactus 3 */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

// #include "SpaceMask.h"

#include "cctk_DefineThorn.h"
#include "rnsid.h"
#include "consts.h"
#include "equil.h"
#include "equil_util.h"
#include "rnsid_util.h"
#include "hdf5_save.h"

#include <loop.hxx>
#include <loop_device.hxx>
#include <rnsreader_utils.hxx>

namespace RNSReader {
using namespace Loop;

static const char *rcsid = "$Header$";
CCTK_FILEVERSION(Hydro_RNSID_rnsid_rfr_c)

static inline int exists_file_name(const char *fname) {
  FILE *file;
  if (file = fopen(fname, "r")) {
    fclose(file);
    return 1;
  }
  return 0;
}

// void Hydro_rnsid_init(CCTK_ARGUMENTS);

extern "C" void Hydro_rnsid_init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_Hydro_rnsid_init;
  DECLARE_CCTK_PARAMETERS;

  /* EQUILIBRIUM VARIABLES */

  int n_tab,       /* Number of points in EOS file */
      print_dif;   /* if =1, monitor convergence */

  double log_e_tab[MAX_NTAB], /* energy dens./c^2 in tab. EOS */
      log_p_tab[MAX_NTAB],    /* pressure in tabulated EOS */
      log_h_tab[MAX_NTAB],    /* enthalpy in EOS file */
      log_n0_tab[MAX_NTAB],   /* number density in EOS file */
      e_center,               /* central energy density */
      p_center,               /* central pressure */
      h_center,               /* central enthalpy */
      e_surface,              /* surface en. density */
      *s_gp,                  /* s grid points */
      *mu,                    /* \mu grid points */
      **rho_potential,        /* potential \rho_potential */
      **gama,                 /* potential \gamma */
      **omega,                /* potential \omega */
      **alpha,                /* potential \alpha */
      **energy,               /* energy density \epsilon */
      **pressure,             /* pressure */
      **enthalpy,             /* enthalpy */
      **velocity_sq,          /* square of velocity */
      R_e,                    /* Circumferential radius */
      Mass,                   /* Gravitational mass */
      Mass_0,                 /* Baryon Mass */
      T,                      /* Rotational kinetic energy */
      W,                      /* Gravitational binding energy */
      Omega,                  /* Angular velocity */
      Omega_K,                /* Ang. vel. of part. in orbit at eq.*/
      r_e,                    /* coord. radius at equator */
      Omega_e,                /* Ang. vel. at equator, when difrot. */
      **Omega_diff,           /* Diff. ang. vel. */
      J;

  CCTK_REAL eos_k, eos_ideal_fluid_gamma, rnsid_rho_min;

  if ((RNS_K < 0.0) || (RNS_Gamma < 0.0)) {
    CCTK_WARN(0,
              "RNS_K and RNS_Gamma must be greater than 0: using 100.0 and 2!");
    eos_k = 100.0;
    eos_ideal_fluid_gamma = 2.0;
  } else {
    eos_k = RNS_K;
    eos_ideal_fluid_gamma = RNS_Gamma;
  }

  rnsid_rho_min = RNS_rho_min;

  char *message;

  /* These used to be params, but are now only to be read from precomputed RNS */

  char *eos_type, *eos_file, *rotation_type;
  eos_type = (char *)malloc(200 * sizeof(char));
  eos_file = (char *)malloc(200 * sizeof(char));
  rotation_type = (char *)malloc(200 * sizeof(char));
  double A_diff, axes_ratio;

  /* INITIAL DATA VARIABLES */

  CCTK_INT m,  /* counter */
      s;       /* counter */

  CCTK_REAL
      **nu,             /* potential nu */
      **B,              /* potential B */
      **rho_0,          /* rest mass density */
      **nu_dr,          /* r-der. in s-coord. of nu */
      **B_dr,           /* r-der. in s-coord. of B */
      **alpha_dr,       /* r-der. in s-coord. of alpha */
      **omega_dr,       /* r-der. in s-coord. of omega */
      **nu_dth,         /* theta-der. in mu-coord. of nu */
      **B_dth,          /* theta-der. in mu-coord. of B */
      **alpha_dth,      /* theta-der. in mu-coord. of alpha */
      **omega_dth,      /* theta-der. in mu-coord. of omega */
      x_i,              /* x at i */
      y_j,              /* y at j */
      z_k,              /* z at k */
      nu_ijk,           /* nu at ijk point */
      exp_nu_ijk,       /* exp(nu) at ijk point */
      B_ijk,            /* B at ijk point */
      omega_ijk,        /* omega at ijk point */
      alpha_ijk,        /* alpha at ijk point */
      exp_alpha_ijk,    /* exp(alpha) at ijk point */
      rho_0_ijk,        /* rho_0 at ijk point */
      energy_ijk,       /* energy at ijk point */
      pressure_ijk,     /* pressure at ijk point */
      nu_dx,            /* derivative of nu w.r.t. x */
      nu_dy,            /* derivative of nu w.r.t. y */
      B_dx,             /* derivative of B w.r.t. x */
      B_dy,             /* derivative of B w.r.t. y */
      omega_dx,         /* derivative of omega w.r.t. x */
      omega_dy,         /* derivative of omega w.r.t. y */
      omega_dz,         /* derivative of omega w.r.t. z */
      alpha_dx,         /* derivative of alpha w.r.t. x */
      alpha_dy,         /* derivative of alpha w.r.t. y */
      r_ijk,            /* r at ijk point */
      r_bar_ijk,        /* sqrt(x^2+y^2) at ijk point */
      dr_dx,            /* dr/dx */
      dr_dy,            /* dr/dy */
      dr_dz,            /* dr/dz */
      dtheta_dx,        /* dtheta/dx */
      dtheta_dy,        /* dtheta/dy */
      dtheta_dz,        /* dtheta/dz */
      nu_dr_ijk,        /* dnu/dr at ijk */
      B_dr_ijk,         /* dB/dr at ijk */
      alpha_dr_ijk,     /* dalpha/dr at ijk */
      omega_dr_ijk,     /* domega/dr at ijk */
      nu_dtheta_ijk,    /* dnu/dtheta at ijk */
      B_dtheta_ijk,     /* dB/dtheta at ijk */
      alpha_dtheta_ijk, /* dalpha/dtheta at ijk */
      omega_dtheta_ijk, /* domega/dtheta at ijk */
      gamma_ijk,        /* gamma = det(3g) */
      h_ijk,            /* h = 1 + eps + P/rho_potential */
      distance_ijk = 0, /* Signed distance to surface */
      e_atm,            /* energy density of atmosphere */
      p_atm,            /* pressure of atmosphere */
      rho_0_atm,        /* rest mass density of atmosphere */
      dens_atm,         /* D of atmosphere */
      tau_atm,          /* tau of atmosphere */
      Omega_ijk;

  int nx = cctkGH->cctk_lsh[0];
  int ny = cctkGH->cctk_lsh[1];
  int nz = cctkGH->cctk_lsh[2];

  CCTK_REAL rho0_center;

 /*
 *
 *  HISTORICAL NOTE ON NAMES OF VARIABLES:
 *
 *   old name               new name (from version 1.25 on)
 *
 *    pert_amp               pert_amplitude (now a parameter)
 *     Gamma_P                eos_ideal_fluid_gamma (now a parameter)
 *      r_ratio                axes_ratio (now a parameter)
 *       rho                    rho_potential
 *
 *        */


  /* COMPUTE POLYTROPIC INDEX AND CENTRAL ENERGY DENSITY */
  e_center = (eos_k*pow(rho_central,eos_ideal_fluid_gamma)/(eos_ideal_fluid_gamma-1.0)+rho_central);

  /* TABULATED EOS OPTION */
  if(strcmp(eos_type,"tab")==0) {
    /* --V0-- load_eos( eos_file, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, Gamma_tab, &n_tab ); */
    int n_nearest;
    double n0;
    /* ==================================================== */
    /* printf(" TAB eos from file: %s\n",eos_file); */
    message = (char *)malloc(200*sizeof(char));
    sprintf(message," TAB eos from file: %s",eos_file);
    CCTK_INFO(message);
    free(message);
    /* ==================================================== */
    load_eos( eos_file, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, &n_tab );
    n_nearest = 50;
    n0 = rho_central/(MB*cactusM);
    e_center = pow(10.0,interp(log_n0_tab, log_e_tab, n_tab,log10(n0), &n_nearest));
  }

  /* SET UP GRID */
  s_gp=(double *)malloc((SDIV+1)*sizeof(double));
  mu=(double *)malloc((MDIV+1)*sizeof(double));
  make_grid(s_gp, mu);

  /* ALLLOCATE MEMORY */

  rho_potential = array_allocate(1,SDIV,1,MDIV);
  gama = array_allocate(1,SDIV,1,MDIV);
  alpha = array_allocate(1,SDIV,1,MDIV);
  omega = array_allocate(1,SDIV,1,MDIV);
  energy = array_allocate(1,SDIV,1,MDIV);
  pressure = array_allocate(1,SDIV,1,MDIV);
  enthalpy = array_allocate(1,SDIV,1,MDIV);
  velocity_sq = array_allocate(1,SDIV,1,MDIV);
  Omega_diff = array_allocate(1,SDIV,1,MDIV);

  /* INITIALIZE VARIABLES WITH ZERO */

  #pragma omp parallel for
  for(s=1;s<=SDIV;s++)
    for(m=1;m<=MDIV;m++) {
      rho_potential[s][m] = 0.0e0;
      gama[s][m] = 0.0e0;
      alpha[s][m] = 0.0e0;
      omega[s][m] = 0.0e0;
      energy[s][m] = 0.0e0;
      pressure[s][m] = 0.0e0;
      enthalpy[s][m] = 0.0e0;
      velocity_sq[s][m] = 0.0e0;
      Omega_diff[s][m] = 0.0e0;
    }

  /* SET DEFAULT EQUILIBRIUM PARAMETERS */

  if(strcmp(eos_type,"tab")==0) {
    e_surface=7.8e-15;
  }
  else {
        e_surface=0.0;
  }

  Omega_e=0.0; /* initialize ang. vel. at equator for diff. rot. */

  print_dif=1;

  /* CHANGE e_center TO POLYTROPIC DIMENSIONLESS UNITS */
  /*      e_center /= ( 1.0/pow(eos_k, n_P) ); */

  /* MAKE e_center DIMENSIONLESS FOR TAB. EOS */

  // if(strcmp(eos_type,"tab")==0) 
  //   e_center *= (C*C*KSCALE);


  /* COMPUTE DIMENSIONLESS CENTRAL PRESSURE AND ENTHALPY */

  /*-V0-- make_center( e_center, log_e_tab, log_p_tab, log_h_tab, n_tab,      */
  /*-V0--             eos_type, eos_ideal_fluid_gamma, &p_center, &h_center); */
  make_center( e_center, log_e_tab, log_p_tab, log_h_tab, n_tab,
  eos_type, eos_k,eos_ideal_fluid_gamma, &p_center, &h_center);

  rho0_center =  (e_center+p_center)*exp(-h_center);

  /* SET-UP INITIAL DATA */

  /* I will only port over the data reading portion. Keepin it simple for now */

  if (print_dif == 1) {
    CCTK_INFO(" ****************************************************");
    CCTK_INFO(" ****************************************************");
    CCTK_INFO(" **              HYDRO - RNSID                     **");
    CCTK_INFO(" **      ROTATING NEUTRON STAR INITIAL DATA        **");
    CCTK_INFO(" ****************************************************");
    CCTK_INFO(" ****************************************************");
  }

  if (strcmp(recover_2Dmodel, "yes") == 0 &&
      exists_file_name(model2D_file) == 1) {

    /* RECOVER FROM 2D FILE */
    /* ==================================================== */
    message = (char *)malloc(200 * sizeof(char));
    sprintf(message, " Recovering 2D model form file %s", model2D_file);
    CCTK_INFO(message);
    free(message);
    /* ==================================================== */

    int sdiv, mdiv;
    hdf5_read_var(&sdiv, &mdiv, model2D_file, eos_type, eos_file, &eos_k,
                  &eos_ideal_fluid_gamma, rotation_type, &A_diff, &axes_ratio,
                  &rho0_center, &r_e, s_gp, mu, rho_potential, gama, alpha,
                  omega, energy, pressure, enthalpy, velocity_sq, &Omega,
                  &Omega_e, Omega_diff);
  } /* END RECOVER FROM 2D FILE, so also end of finding the equil state */

  /* COMPUTE EQUILIBRIUM QUANTITIES (Mass, Radius, T/W etc.) */

  comp_values(s_gp, mu, axes_ratio, e_surface, r_e, eos_type, log_e_tab,
              log_n0_tab, n_tab, Omega, rho_potential, gama, alpha, omega,
              energy, pressure, enthalpy, velocity_sq, &Mass, &Mass_0, &T, &W,
              &Omega_K, &R_e, rotation_type, Omega_diff, &J);

  /* TRANSFORM UNITS TO c=G=M_sun=1 */
  /*
  transform_units( eos_type, n_P, eos_k, &rho0_center, &e_center,
                                                                   &p_center,
  &r_e, omega, energy, pressure, &Mass, &Mass_0, &T, &W, &Omega, &Omega_K, &R_e,
                                                                   &Omega_e,
  Omega_diff, &J);
  */

  /* RETURN OMEGA AND R_E VALUES */

  /* PRINT-OUT SOME EQUILIBRIUM QUANTITIES */
  if (print_dif == 1) {
    message = (char *)malloc(200 * sizeof(char));
    sprintf(message, " %5.4e %5.4e %5.4e %5.4e %5.4e", rho0_center, e_center,
            Mass, Mass_0, R_e);
    CCTK_INFO("Equilibrium model done in c=G=M_sun=1 dimensionless form");
    CCTK_INFO(" rho_center   e_center    Mass      Mass_0      R_e");
    CCTK_INFO(message);
    if (strcmp(rotation_type, "uniform") == 0) {
      CCTK_INFO(
          "     J         T/W       Omega   Omega_Kepler axes_ratio  J/M^2 ");
      sprintf(message, " %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e",
              ((Omega > 0.0) ? J : 0.0), T / W, Omega, Omega_K, axes_ratio,
              J / (Mass * Mass));
      CCTK_INFO(message);
    } else {
      CCTK_INFO("     J         T/W       Omega_c Omega_e   Omega_Kepler "
                "axes_ratio  J/M^2 ");
      sprintf(message, " %5.4e %5.4e %5.4e %5.4e %5.4e  %5.4e %5.4e",
              ((Omega > 0.0) ? J : 0.0), T / W, Omega, Omega_e, Omega_K,
              axes_ratio, J / (Mass * Mass));
      CCTK_INFO(message);
    }
    free(message);
  }

  /* CONSTRUCT ARRAYS WITH NEEDED POLAR QUANTITIES */

  nu = array_allocate(1, SDIV, 1, MDIV);
  B = array_allocate(1, SDIV, 1, MDIV);
  rho_0 = array_allocate(1, SDIV, 1, MDIV);

  for (m = 1; m <= MDIV; m++)
    for (s = 1; s <= SDIV; s++) {
      nu[s][m] = (gama[s][m] + rho_potential[s][m]) / 2.0;
      B[s][m] = exp(gama[s][m]);
      rho_0[s][m] = (energy[s][m] + pressure[s][m]) * exp(-enthalpy[s][m]);
    }

  array_free(rho_potential, 1, SDIV, 1, MDIV);
  array_free(gama, 1, SDIV, 1, MDIV);
  array_free(enthalpy, 1, SDIV, 1, MDIV);
  array_free(velocity_sq, 1, SDIV, 1, MDIV);

  nu_dr = array_allocate(1, SDIV, 1, MDIV);
  B_dr = array_allocate(1, SDIV, 1, MDIV);
  alpha_dr = array_allocate(1, SDIV, 1, MDIV);
  omega_dr = array_allocate(1, SDIV, 1, MDIV);
  nu_dth = array_allocate(1, SDIV, 1, MDIV);
  B_dth = array_allocate(1, SDIV, 1, MDIV);
  alpha_dth = array_allocate(1, SDIV, 1, MDIV);
  omega_dth = array_allocate(1, SDIV, 1, MDIV);

  for (m = 1; m <= MDIV; m++)
    for (s = 1; s <= SDIV; s++) {
      nu_dr[s][m] = deriv_s(nu, s, m) * SQ(1.0 - s_gp[s]) / r_e;
      B_dr[s][m] = deriv_s(B, s, m) * SQ(1.0 - s_gp[s]) / r_e;
      alpha_dr[s][m] = deriv_s(alpha, s, m) * SQ(1.0 - s_gp[s]) / r_e;
      omega_dr[s][m] = deriv_s(omega, s, m) * SQ(1.0 - s_gp[s]) / r_e;
      nu_dth[s][m] = deriv_m(nu, s, m) * (-sqrt(1.0 - SQ(mu[m])));
      B_dth[s][m] = deriv_m(B, s, m) * (-sqrt(1.0 - SQ(mu[m])));
      alpha_dth[s][m] = deriv_m(alpha, s, m) * (-sqrt(1.0 - SQ(mu[m])));
      omega_dth[s][m] = deriv_m(omega, s, m) * (-sqrt(1.0 - SQ(mu[m])));
    }

  /* COMPUTE INITIAL DATA */

  rho_0_atm = rnsid_rho_min; /* rename the constant for historical reasons */
  e_atm = rho_0_atm;
  p_atm = eos_k * pow(rho_0_atm, eos_ideal_fluid_gamma);

  grid.loop_all<1, 1, 1>(
      grid.nghostzones, [&](const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        x_i = p.x;
        y_j = p.y;
        z_k = p.z;

        grid_interp_all(s_gp, mu, r_e, nx, ny, nz,
                        // x, y, z,
                        // i, j, k,
                        p, nu, B, alpha, omega, nu_dr, B_dr, alpha_dr, omega_dr,
                        nu_dth, B_dth, alpha_dth, omega_dth, rho_0, energy,
                        pressure, &nu_ijk, &B_ijk, &alpha_ijk, &omega_ijk,
                        &nu_dr_ijk, &B_dr_ijk, &alpha_dr_ijk, &omega_dr_ijk,
                        &nu_dtheta_ijk, &B_dtheta_ijk, &alpha_dtheta_ijk,
                        &omega_dtheta_ijk, &rho_0_ijk, &energy_ijk,
                        &pressure_ijk, &distance_ijk, Omega_diff, &Omega_ijk);

        /* *************************************** */
        /* DETECT if it is in the ATMOSPHERE       */
        /* *************************************** */
        if ((rho_0_ijk <= 0.0) || (energy_ijk <= 0.0) ||
            (pressure_ijk <= 0.0)) {

          rho_0_ijk = rho_0_atm;
          energy_ijk = e_atm + 1.e-20;
          pressure_ijk = p_atm;
        }
        /* *************************************** */
        /* END ATMOSPERE SETTINGs                      */
        /* *************************************** */

        exp_nu_ijk = exp(nu_ijk);
        exp_alpha_ijk = exp(alpha_ijk);

        r_ijk = sqrt(SQ(x_i) + SQ(y_j) + SQ(z_k));
        r_bar_ijk = sqrt(SQ(x_i) + SQ(y_j));

        alp_cc(p.I) = exp_nu_ijk;

        if (x_i == 0.0 && y_j == 0.0) {

          gxx_cc(p.I) = SQ(exp_alpha_ijk);
          gyy_cc(p.I) = SQ(exp_alpha_ijk);
          gzz_cc(p.I) = SQ(exp_alpha_ijk);

          gxy_cc(p.I) = 0.0;
          gxz_cc(p.I) = 0.0;
          gyz_cc(p.I) = 0.0;

          kxx_cc(p.I) = 0.0;
          kyy_cc(p.I) = 0.0;
          kzz_cc(p.I) = 0.0;
          kxy_cc(p.I) = 0.0;
          kxz_cc(p.I) = 0.0;
          kyz_cc(p.I) = 0.0;

        } else {
          dr_dx = x_i / r_ijk;
          dr_dy = y_j / r_ijk;
          dr_dz = z_k / r_ijk;

          dtheta_dx = x_i * z_k / (SQ(r_ijk) * r_bar_ijk);
          dtheta_dy = y_j * z_k / (SQ(r_ijk) * r_bar_ijk);
          dtheta_dz = -r_bar_ijk / SQ(r_ijk);

          nu_dx = dr_dx * nu_dr_ijk + dtheta_dx * nu_dtheta_ijk;
          nu_dy = dr_dy * nu_dr_ijk + dtheta_dy * nu_dtheta_ijk;

          B_dx = dr_dx * B_dr_ijk + dtheta_dx * B_dtheta_ijk;
          B_dy = dr_dy * B_dr_ijk + dtheta_dy * B_dtheta_ijk;

          alpha_dx = dr_dx * alpha_dr_ijk + dtheta_dx * alpha_dtheta_ijk;
          alpha_dy = dr_dy * alpha_dr_ijk + dtheta_dy * alpha_dtheta_ijk;

          omega_dx = dr_dx * omega_dr_ijk + dtheta_dx * omega_dtheta_ijk;
          omega_dy = dr_dy * omega_dr_ijk + dtheta_dy * omega_dtheta_ijk;

          /* enforce omega_dz=0 at z=0 (it is slightly nonzero due
                   to O(h) forwards formula in computing derivative) */
          if (z_k == 0.0)
            omega_dz = 0.0;
          else
            omega_dz = dr_dz * omega_dr_ijk + dtheta_dz * omega_dtheta_ijk;

          gxx_cc(p.I) =
              (SQ(B_ijk * y_j / exp_nu_ijk) + SQ(exp_alpha_ijk * x_i)) /
              (SQ(x_i) + SQ(y_j));

          gxy_cc(p.I) = (SQ(exp_alpha_ijk) - SQ(B_ijk / exp_nu_ijk)) * x_i *
                        y_j / (SQ(x_i) + SQ(y_j));

          gxz_cc(p.I) = 0.0;

          gyy_cc(p.I) =
              (SQ(B_ijk * x_i / exp_nu_ijk) + SQ(exp_alpha_ijk * y_j)) /
              (SQ(x_i) + SQ(y_j));

          gyz_cc(p.I) = 0.0;

          gzz_cc(p.I) = SQ(exp_alpha_ijk);

          kxx_cc(p.I) =
              ((SQ(r_bar_ijk) * y_j * omega_dx +
                (x_i * nu_dy - y_j * nu_dx) * SQ(y_j) * omega_ijk) *
                   SQ(B_ijk) +
               (y_j * B_dx - x_i * B_dy) * omega_ijk * SQ(y_j) * B_ijk +
               (y_j * alpha_dx - x_i * alpha_dy) * omega_ijk *
                   SQ(x_i * exp_alpha_ijk * exp_nu_ijk)) /
              (SQ(r_bar_ijk * exp_nu_ijk) * exp_nu_ijk);

          kxy_cc(p.I) =
              ((0.5 * SQ(r_bar_ijk) * (y_j * omega_dy - x_i * omega_dx) +
                (y_j * nu_dx - x_i * nu_dy) * x_i * y_j * omega_ijk) *
                   SQ(B_ijk) +
               (-y_j * B_dx + x_i * B_dy) * omega_ijk * x_i * y_j * B_ijk +
               (y_j * alpha_dx - x_i * alpha_dy) * omega_ijk * x_i * y_j *
                   SQ(exp_alpha_ijk * exp_nu_ijk)) /
              (SQ(r_bar_ijk * exp_nu_ijk) * exp_nu_ijk);

          kxz_cc(p.I) =
              0.5 * SQ(B_ijk) * y_j * omega_dz / (SQ(exp_nu_ijk) * exp_nu_ijk);

          kyy_cc(p.I) =
              ((-SQ(r_bar_ijk) * x_i * omega_dy +
                (x_i * nu_dy - y_j * nu_dx) * SQ(x_i) * omega_ijk) *
                   SQ(B_ijk) +
               (y_j * B_dx - x_i * B_dy) * omega_ijk * SQ(x_i) * B_ijk +
               (y_j * alpha_dx - x_i * alpha_dy) * omega_ijk *
                   SQ(y_j * exp_alpha_ijk * exp_nu_ijk)) /
              (SQ(r_bar_ijk * exp_nu_ijk) * exp_nu_ijk);

          kyz_cc(p.I) =
              -0.5 * SQ(B_ijk) * x_i * omega_dz / (SQ(exp_nu_ijk) * exp_nu_ijk);

          kzz_cc(p.I) = (y_j * alpha_dx - x_i * alpha_dy) * omega_ijk *
                        SQ(exp_alpha_ijk) / exp_nu_ijk;
        }

        betax_cc(p.I) = omega_ijk * y_j;
        betay_cc(p.I) = -omega_ijk * x_i;

        betaz_cc(p.I) = 0.0;

        rho(p.I) = rho_0_ijk;

        eps(p.I) = energy_ijk / rho_0_ijk - 1.0;
        h_ijk = (energy_ijk + pressure_ijk) / rho_0_ijk;

        gamma_ijk = SQ(exp_alpha_ijk) * SQ(B_ijk * exp_alpha_ijk / exp_nu_ijk);

        // W_ijk = 1.0/sqrt(1.0-SQ((omega_ijk-Omega_ijk)*B_ijk*
        // 				r_bar_ijk/SQ(exp_nu_ijk)));

        // w_lorentz(p.I) = W_ijk;

        velx(p.I) = (omega_ijk - Omega_ijk) * y_j / exp_nu_ijk;

        vely(p.I) = -(omega_ijk - Omega_ijk) * x_i / exp_nu_ijk;

        velz(p.I) = 0.0;

        press(p.I) = pressure_ijk;

        dens_atm = sqrt(gamma_ijk) * rho_0_atm;
        tau_atm = sqrt(gamma_ijk) * eos_k *
                  pow(rho_0_atm, eos_ideal_fluid_gamma) /
                  (eos_ideal_fluid_gamma - 1.0);

        /* *************************************** */
        /* ATMOSPERE SETTINGs                      */
        /* *************************************** */
        if ((rho(p.I) < (1.0 + RNS_atmo_tolerance) * rho_0_atm)) {
          const CCTK_REAL radial_distance = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
          CCTK_REAL rho_atm = (radial_distance > r_atmo)
            ? (rho_abs_min * pow((r_atmo / radial_distance), n_rho_atmo))
            : rho_abs_min;
          rho(p.I) = rho_atm; // This will be overwritten by C2P anyway

          velx(p.I) = 0.0;
          vely(p.I) = 0.0;
          velz(p.I) = 0.0;
        }
        /* *************************************** */
        /* ATMOSPERE SETTINGs                      */
        /* *************************************** */
      }); /* END GRID LOOP */

  if (strcmp(zero_shift, "yes") == 0) {

    /* SET SHIFT TO ZERO */
    grid.loop_all<1, 1, 1>(grid.nghostzones, [&](const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   betax_cc(p.I) = 0.0;
                                                   betay_cc(p.I) = 0.0;
                                                   betaz_cc(p.I) = 0.0;
                                                 });
  }

  /* FREE MEMORY */

  array_free(alpha, 1, SDIV, 1, MDIV);
  array_free(omega, 1, SDIV, 1, MDIV);
  array_free(rho_0, 1, SDIV, 1, MDIV);
  array_free(energy, 1, SDIV, 1, MDIV);
  array_free(pressure, 1, SDIV, 1, MDIV);

  array_free(nu, 1, SDIV, 1, MDIV);
  array_free(B, 1, SDIV, 1, MDIV);

  array_free(nu_dr, 1, SDIV, 1, MDIV);
  array_free(B_dr, 1, SDIV, 1, MDIV);
  array_free(alpha_dr, 1, SDIV, 1, MDIV);
  array_free(omega_dr, 1, SDIV, 1, MDIV);
  array_free(nu_dth, 1, SDIV, 1, MDIV);
  array_free(B_dth, 1, SDIV, 1, MDIV);
  array_free(alpha_dth, 1, SDIV, 1, MDIV);
  array_free(omega_dth, 1, SDIV, 1, MDIV);

  free(s_gp);
  free(mu);

  /* Add code for filling time levels here if needed in the future. */

  return;
}

/*
extern "C" void Hydro_RNSID_CheckParameters(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (timelevels < 2)
  {
      CCTK_PARAMWARN("You have to set 'HydroBase::timelevels to at least 2");
  }
}
*/

extern "C" void RNSReader_Interpolation_C2V(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_RNSReader_Interpolation_C2V;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO("Starting interpolation for ADM variables.");
  grid.loop_int<0, 0, 0>(grid.nghostzones,
                         [=] CCTK_HOST(const Loop::PointDesc &p)
                             CCTK_ATTRIBUTE_ALWAYS_INLINE {
                               gxx(p.I) = calc_avg_c2v(gxx_cc, p);
                               gxy(p.I) = calc_avg_c2v(gxy_cc, p);
                               gxz(p.I) = calc_avg_c2v(gxz_cc, p);
                               gyy(p.I) = calc_avg_c2v(gyy_cc, p);
                               gyz(p.I) = calc_avg_c2v(gyz_cc, p);
                               gzz(p.I) = calc_avg_c2v(gzz_cc, p);

                               kxx(p.I) = calc_avg_c2v(kxx_cc, p);
                               kxy(p.I) = calc_avg_c2v(kxy_cc, p);
                               kxz(p.I) = calc_avg_c2v(kxz_cc, p);
                               kyy(p.I) = calc_avg_c2v(kyy_cc, p);
                               kyz(p.I) = calc_avg_c2v(kyz_cc, p);
                               kzz(p.I) = calc_avg_c2v(kzz_cc, p);

                               alp(p.I) = calc_avg_c2v(alp_cc, p);
                               betax(p.I) = calc_avg_c2v(betax_cc, p);
                               betay(p.I) = calc_avg_c2v(betay_cc, p);
                               betaz(p.I) = calc_avg_c2v(betaz_cc, p);

                               dtalp(p.I) = calc_avg_c2v(dtalp_cc, p);
                               dtbetax(p.I) = calc_avg_c2v(dtbetax_cc, p);
                               dtbetay(p.I) = calc_avg_c2v(dtbetay_cc, p);
                               dtbetaz(p.I) = calc_avg_c2v(dtbetaz_cc, p);
                             });

  CCTK_INFO("Done interpolation for ADM variables.");
}

} // namespace RNSReader
