
/* This function converts TOV quantities into
 * ADM quantities in Spherical coordinates.
 */
void convert_TOV_spacetime_vars_to_ADM_vars( const CCTK_REAL rr, const CCTK_REAL th,
                           const CCTK_REAL IDexp_4phi, const CCTK_REAL IDexpnu,
                           CCTK_REAL *IDalpha,
                           CCTK_REAL *IDgammaDD00, CCTK_REAL *IDgammaDD01, CCTK_REAL *IDgammaDD02,
                           CCTK_REAL *IDgammaDD11, CCTK_REAL *IDgammaDD12, CCTK_REAL *IDgammaDD22) {

  /***************************************************************
   * Convert TOV quantities to ADM quantities in Spherical basis *
   ***************************************************************
   *
   * First we convert the lapse function:
   * .------------------.
   * | alpha = e^(nu/2) |
   * .------------------.
   */
   *IDalpha = sqrt(IDexpnu);

  /* Next we convert the metric function:
   * .----------------------------------------.
   * | gamma_{00} = e^{4phi}                  |
   * .----------------------------------------.
   * | gamma_{11} = e^{4phi} r^2              |
   * .----------------------------------------.
   * | gamma_{22} = e^{4phi} r^2 sin^2(theta) |
   * .----------------------------------------.
   * | All other components are zero.         |
   * .----------------------------------------.
   */
   *IDgammaDD00 = IDexp_4phi;
   *IDgammaDD11 = IDexp_4phi * rr * rr;
   *IDgammaDD22 = IDexp_4phi * rr * rr * sin(th) * sin(th);
   *IDgammaDD01 = 0.0;
   *IDgammaDD02 = 0.0;
   *IDgammaDD12 = 0.0;

}
