

namespace NRPyPlusTOVID{
using namespace std;
using namespace Loop;

static inline
void HydroQuantities(const PointDesc &p,
                   const CCTK_REAL IDPressure, const CCTK_REAL IDrho_baryonic,
                   const CCTK_REAL IDrho__total_energy_density,
                   const GF3D2<CCTK_REAL> &PressureGF,
                   const GF3D2<CCTK_REAL> &rho_baryonicGF,
                   const GF3D2<CCTK_REAL> &epsilonGF,
                   const GF3D2<CCTK_REAL> &Valencia3velocityU0GF,
                   const GF3D2<CCTK_REAL> &Valencia3velocityU1GF,
                   const GF3D2<CCTK_REAL> &Valencia3velocityU2GF) {
    DECLARE_CCTK_PARAMETERS;
    if(IDrho__total_energy_density <= 0 || IDrho_baryonic <= 0 || IDPressure <= 0) {
        rho_baryonicGF(p.I) = rho_atmosphere;
        PressureGF(p.I)     = K_atmosphere*pow(rho_atmosphere,Gamma_atmosphere);
        epsilonGF(p.I)      = 0;
        Valencia3velocityU0GF(p.I) = 0;
        Valencia3velocityU1GF(p.I) = 0;
        Valencia3velocityU2GF(p.I) = 0;
    } else {
      /*
       * NRPy+ Finite Difference Code Generation, Step 1 of 1: Evaluate SymPy expressions and write to main memory:
       */
      PressureGF(p.I) = IDPressure;
      rho_baryonicGF(p.I) = IDrho_baryonic;
      epsilonGF(p.I) = IDrho__total_energy_density/IDrho_baryonic - 1;
      Valencia3velocityU0GF(p.I) = 0;
      Valencia3velocityU1GF(p.I) = 0;
      Valencia3velocityU2GF(p.I) = 0;

        // Apply pressure depletion.
        PressureGF(p.I) *= (1.0 - Pressure_depletion_factor);
    }
}
    
} // namespace
