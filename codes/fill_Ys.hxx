/* fill_Ys: implement stage 2-3 of section C.1 of
 * https://github.com/lwJi/Docs-Subcycling/blob/main/notes/MC/MC.pdf
 * please take a look at AMReX-Codes/amrex/Src/AmrCore/AMReX_FillPatcher.H
 * for more information */
void fill_Ys(
        /* output: fine Y_i, depend on which stage (1,2,3,4))*/
        CCTK_REAL Yi,
        /* input: which RK stage (1,2,3,4) */
        CCTK_INT stage
        /* input: coarse k_i^(c) */
        CCTK_REAL k1c,
        CCTK_REAL k2c,
        CCTK_REAL k3c,
        CCTK_REAL k4c,
        /* input: coarse y_n^(c) */
        CCTK_REAL yn,
        /* input: coarse dt */
        CCTK_REAL dtc,
        /* input: which iteration for fine grid, 1-first, 2-second */
        CCTK_INT iteration
) {
  int rk_order = 4;
  CCTK_INT ncycle = 2; /* 2-to-1 time refinement */

  CCTK_REAL r = CCTK_REAL(1) / CCTK_REAL(ncycle);
  CCTK_REAL theta = CCTK_REAL(iteration - 1) / CCTK_REAL(ncycle);

  if (rk_order == 4) {
    CCTK_REAL theta2 = theta * theta;
    CCTK_REAL theta3 = theta2 * theta;
    // coefficients for U
    CCTK_REAL b1 = theta - CCTK_REAL(1.5) * theta2 + CCTK_REAL(2. / 3.) * theta3;
    CCTK_REAL b2 = theta2 - CCTK_REAL(2. / 3.) * theta3;
    CCTK_REAL b3 = b2;
    CCTK_REAL b4 = CCTK_REAL(-0.5) * theta2 + CCTK_REAL(2. / 3.) * theta3;
    // coefficients for Ut
    CCTK_REAL c1 = CCTK_REAL(1.) - CCTK_REAL(3.) * theta + CCTK_REAL(2.) * theta2;
    CCTK_REAL c2 = CCTK_REAL(2.) * theta - CCTK_REAL(2.) * theta2;
    CCTK_REAL c3 = c2;
    CCTK_REAL c4 = -theta + CCTK_REAL(2.) * theta2;
    // coefficients for Utt
    CCTK_REAL d1 = CCTK_REAL(-3.) + CCTK_REAL(4.) * theta;
    CCTK_REAL d2 = CCTK_REAL(2.) - CCTK_REAL(4.) * theta;
    CCTK_REAL d3 = d2;
    CCTK_REAL d4 = CCTK_REAL(-1.) + CCTK_REAL(4.) * theta;
    // coefficients for Uttt
    constexpr CCTK_REAL e1 = CCTK_REAL(4.);
    constexpr CCTK_REAL e2 = CCTK_REAL(-4.);
    constexpr CCTK_REAL e3 = CCTK_REAL(-4.);
    constexpr CCTK_REAL e4 = CCTK_REAL(4.);

    if (stage == 1) {
      CCTK_REAL uu = b1 * k1c + b2 * k2c + b3 * k3c + b4 * k4c;
      Yi = yn + dtc * uu;
    } else if (stage == 2) {
      CCTK_REAL uu = b1 * k1c + b2 * k2c + b3 * k3c + b4 * k4c;
      CCTK_REAL ut = c1 * k1c + c2 * k2c + c3 * k3c + c4 * k4c;
      Yi = yn + dtc * (uu + CCTK_REAL(0.5) * r * ut);
    } else if (stage == 3 || stage == 4) {
      CCTK_REAL r2 = r * r;
      CCTK_REAL r3 = r2 * r;
      CCTK_REAL at = (stage == 3) ? CCTK_REAL(0.5) * r : r;
      CCTK_REAL att = (stage == 3) ? CCTK_REAL(0.25) * r2 : CCTK_REAL(0.5) * r2;
      CCTK_REAL attt = (stage == 3) ? CCTK_REAL(0.0625) * r3 : CCTK_REAL(0.125) * r3;
      CCTK_REAL akk = (stage == 3) ? CCTK_REAL(-4.) : CCTK_REAL(4.);
      CCTK_REAL uu = b1 * k1c + b2 * k2c + b3 * k3c + b4 * k4c;
      CCTK_REAL ut = c1 * k1c + c2 * k2c + c3 * k3c + c4 * k4c;
      CCTK_REAL utt = d1 * k1c + d2 * k2c + d3 * k3c + d4 * k4c;
      CCTK_REAL uttt = e1 * k1c + e2 * k2c + e3 * k3c + e4 * k4c;
      Yi = yn + dtc * (uu + at * ut + att * utt + attt * (uttt + akk * (k3c - k2c)));
    }
  }

}
