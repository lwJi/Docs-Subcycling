void fill_Ys(int stage, int iteration, int ncycle, MF &mf, Real time, BC &cbc,
             BC &fbc, Vector<BCRec> const &bcs) {
  int rk_order = m_cf_crse_data.size() - 1;
  if (rk_order != 3 && rk_order != 4) {
    amrex::Abort("FillPatcher: unsupported RK order " +
                 std::to_string(rk_order));
    return;
  }
  AMREX_ASSERT(stage > 0 && stage <= rk_order);

  auto const &fpc = getFPinfo();
  if (m_cf_crse_data_tmp == nullptr) {
    m_cf_crse_data_tmp =
        std::make_unique<MF>(make_mf_crse_patch<MF>(fpc, m_ncomp));
  }

  auto const &u = m_cf_crse_data_tmp->arrays();
  auto const &u0 = m_cf_crse_data[0].second->const_arrays();
  auto const &k1 = m_cf_crse_data[1].second->const_arrays();
  auto const &k2 = m_cf_crse_data[2].second->const_arrays();
  auto const &k3 = m_cf_crse_data[3].second->const_arrays();

  Real dtc = m_dt_coarse;
  Real r = Real(1) / Real(ncycle);
  Real xsi = Real(iteration - 1) / Real(ncycle);

  int const ng_space_interp = 8; // Need to be big enough
  Box cdomain = m_cgeom.growPeriodicDomain(ng_space_interp);
  cdomain.convert(m_cf_crse_data_tmp->ixType());

  if (rk_order == 3) {
    // coefficients for U
    Real b1 = xsi - Real(5. / 6.) * xsi * xsi;
    Real b2 = Real(1. / 6.) * xsi * xsi;
    Real b3 = Real(2. / 3) * xsi * xsi;
    // coefficients for Ut
    Real c1 = Real(1.) - Real(5. / 3.) * xsi;
    Real c2 = Real(1. / 3.) * xsi;
    Real c3 = Real(4. / 3.) * xsi;
    // coefficients for Utt
    constexpr Real d1 = Real(-5. / 3.);
    constexpr Real d2 = Real(1. / 3.);
    constexpr Real d3 = Real(4. / 3.);
    if (stage == 1) {
      amrex::ParallelFor(
          *m_cf_crse_data_tmp, IntVect(0), m_ncomp,
          [=] AMREX_GPU_DEVICE(int bi, int i, int j, int k, int n) noexcept {
            if (cdomain.contains(i, j, k)) {
              Real kk1 = k1[bi](i, j, k, n);
              Real kk2 = k2[bi](i, j, k, n);
              Real kk3 = k3[bi](i, j, k, n);
              Real uu = b1 * kk1 + b2 * kk2 + b3 * kk3;
              u[bi](i, j, k, n) = u0[bi](i, j, k, n) + dtc * uu;
            }
          });
    } else if (stage == 2) {
      amrex::ParallelFor(
          *m_cf_crse_data_tmp, IntVect(0), m_ncomp,
          [=] AMREX_GPU_DEVICE(int bi, int i, int j, int k, int n) noexcept {
            if (cdomain.contains(i, j, k)) {
              Real kk1 = k1[bi](i, j, k, n);
              Real kk2 = k2[bi](i, j, k, n);
              Real kk3 = k3[bi](i, j, k, n);
              Real uu = b1 * kk1 + b2 * kk2 + b3 * kk3;
              Real ut = c1 * kk1 + c2 * kk2 + c3 * kk3;
              u[bi](i, j, k, n) = u0[bi](i, j, k, n) + dtc * (uu + r * ut);
            }
          });
    } else if (stage == 3) {
      amrex::ParallelFor(
          *m_cf_crse_data_tmp, IntVect(0), m_ncomp,
          [=] AMREX_GPU_DEVICE(int bi, int i, int j, int k, int n) noexcept {
            if (cdomain.contains(i, j, k)) {
              Real kk1 = k1[bi](i, j, k, n);
              Real kk2 = k2[bi](i, j, k, n);
              Real kk3 = k3[bi](i, j, k, n);
              Real uu = b1 * kk1 + b2 * kk2 + b3 * kk3;
              Real ut = c1 * kk1 + c2 * kk2 + c3 * kk3;
              Real utt = d1 * kk1 + d2 * kk2 + d3 * kk3;
              u[bi](i, j, k, n) =
                  u0[bi](i, j, k, n) +
                  dtc * (uu + Real(0.5) * r * ut + Real(0.25) * r * r * utt);
            }
          });
    }
  } else if (rk_order == 4) {
    auto const &k4 = m_cf_crse_data[4].second->const_arrays();
    Real xsi2 = xsi * xsi;
    Real xsi3 = xsi2 * xsi;
    // coefficients for U
    Real b1 = xsi - Real(1.5) * xsi2 + Real(2. / 3.) * xsi3;
    Real b2 = xsi2 - Real(2. / 3.) * xsi3;
    Real b3 = b2;
    Real b4 = Real(-0.5) * xsi2 + Real(2. / 3.) * xsi3;
    // coefficients for Ut
    Real c1 = Real(1.) - Real(3.) * xsi + Real(2.) * xsi2;
    Real c2 = Real(2.) * xsi - Real(2.) * xsi2;
    Real c3 = c2;
    Real c4 = -xsi + Real(2.) * xsi2;
    // coefficients for Utt
    Real d1 = Real(-3.) + Real(4.) * xsi;
    Real d2 = Real(2.) - Real(4.) * xsi;
    Real d3 = d2;
    Real d4 = Real(-1.) + Real(4.) * xsi;
    // coefficients for Uttt
    constexpr Real e1 = Real(4.);
    constexpr Real e2 = Real(-4.);
    constexpr Real e3 = Real(-4.);
    constexpr Real e4 = Real(4.);
    if (stage == 1) {
      amrex::ParallelFor(
          *m_cf_crse_data_tmp, IntVect(0), m_ncomp,
          [=] AMREX_GPU_DEVICE(int bi, int i, int j, int k, int n) noexcept {
            if (cdomain.contains(i, j, k)) {
              Real kk1 = k1[bi](i, j, k, n);
              Real kk2 = k2[bi](i, j, k, n);
              Real kk3 = k3[bi](i, j, k, n);
              Real kk4 = k4[bi](i, j, k, n);
              Real uu = b1 * kk1 + b2 * kk2 + b3 * kk3 + b4 * kk4;
              u[bi](i, j, k, n) = u0[bi](i, j, k, n) + dtc * uu;
            }
          });
    } else if (stage == 2) {
      amrex::ParallelFor(
          *m_cf_crse_data_tmp, IntVect(0), m_ncomp,
          [=] AMREX_GPU_DEVICE(int bi, int i, int j, int k, int n) noexcept {
            if (cdomain.contains(i, j, k)) {
              Real kk1 = k1[bi](i, j, k, n);
              Real kk2 = k2[bi](i, j, k, n);
              Real kk3 = k3[bi](i, j, k, n);
              Real kk4 = k4[bi](i, j, k, n);
              Real uu = b1 * kk1 + b2 * kk2 + b3 * kk3 + b4 * kk4;
              Real ut = c1 * kk1 + c2 * kk2 + c3 * kk3 + c4 * kk4;
              u[bi](i, j, k, n) =
                  u0[bi](i, j, k, n) + dtc * (uu + Real(0.5) * r * ut);
            }
          });
    } else if (stage == 3 || stage == 4) {
      Real r2 = r * r;
      Real r3 = r2 * r;
      Real at = (stage == 3) ? Real(0.5) * r : r;
      Real att = (stage == 3) ? Real(0.25) * r2 : Real(0.5) * r2;
      Real attt = (stage == 3) ? Real(0.0625) * r3 : Real(0.125) * r3;
      Real akk = (stage == 3) ? Real(-4.) : Real(4.);
      amrex::ParallelFor(
          *m_cf_crse_data_tmp, IntVect(0), m_ncomp,
          [=] AMREX_GPU_DEVICE(int bi, int i, int j, int k, int n) noexcept {
            if (cdomain.contains(i, j, k)) {
              Real kk1 = k1[bi](i, j, k, n);
              Real kk2 = k2[bi](i, j, k, n);
              Real kk3 = k3[bi](i, j, k, n);
              Real kk4 = k4[bi](i, j, k, n);
              Real uu = b1 * kk1 + b2 * kk2 + b3 * kk3 + b4 * kk4;
              Real ut = c1 * kk1 + c2 * kk2 + c3 * kk3 + c4 * kk4;
              Real utt = d1 * kk1 + d2 * kk2 + d3 * kk3 + d4 * kk4;
              Real uttt = e1 * kk1 + e2 * kk2 + e3 * kk3 + e4 * kk4;
              u[bi](i, j, k, n) = u0[bi](i, j, k, n) +
                                  dtc * (uu + at * ut + att * utt +
                                         attt * (uttt + akk * (kk3 - kk2)));
            }
          });
    }
  }
  Gpu::streamSynchronize();

  cbc(*m_cf_crse_data_tmp, 0, m_ncomp, m_nghost, time, 0);

  if (m_cf_fine_data == nullptr) {
    m_cf_fine_data = std::make_unique<MF>(make_mf_fine_patch<MF>(fpc, m_ncomp));
  }

  FillPatchInterp(
      *m_cf_fine_data, 0, *m_cf_crse_data_tmp, 0, m_ncomp, IntVect(0), m_cgeom,
      m_fgeom,
      amrex::grow(amrex::convert(m_fgeom.Domain(), mf.ixType()), m_nghost),
      m_ratio, m_interp, bcs, 0);

  // xxxxx We can optimize away this ParallelCopy by making a special fpinfo.
  mf.ParallelCopy(*m_cf_fine_data, 0, 0, m_ncomp, IntVect(0), m_nghost);

  mf.FillBoundary(m_fgeom.periodicity());
  fbc(mf, 0, m_ncomp, m_nghost, time, 0);
}
