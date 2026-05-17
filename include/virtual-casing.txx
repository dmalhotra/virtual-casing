
template <class Real> VirtualCasing<Real>::VirtualCasing() : comm_(sctl::Comm::Self()), BiotSavartFxU(comm_), LaplaceFxdU(comm_), BiotSavartFxdU(comm_), LaplaceFxd2U(comm_), Svec(1), NFP_(0), digits_(10), dosetup(true), dosetup_grad(true) {
}

template <class Real> void VirtualCasing<Real>::Setup(const sctl::Integer digits, const sctl::Integer NFP, const bool half_period, const sctl::Long surf_Nt, const sctl::Long surf_Np, const std::vector<Real>& X, const sctl::Long src_Nt, const sctl::Long src_Np, const sctl::Long trg_Nt, const sctl::Long trg_Np) {
  dosetup = true;
  dosetup_grad = true;
  SCTL_ASSERT(surf_Nt*surf_Np*COORD_DIM == (sctl::Long)X.size());
  if (half_period) { // upsample surf_Nt by 1
    sctl::Vector<Real> X0, X1;
    Svec[0] = biest::Surface<Real>(NFP*2*(surf_Nt+1), surf_Np, biest::SurfType::None);
    biest::SurfaceOp<Real>::CompleteVecField(X0, true, half_period, NFP, surf_Nt, surf_Np, sctl::Vector<Real>(X), -sctl::const_pi<Real>()/(NFP*surf_Nt*2));
    biest::SurfaceOp<Real>::Resample(X1, NFP*2*(surf_Nt+1), surf_Np, X0, NFP*2*surf_Nt, surf_Np);
    biest::SurfaceOp<Real>::RotateToroidal(Svec[0].Coord(), X1, NFP*2*(surf_Nt+1), surf_Np, sctl::const_pi<Real>()/(NFP*trg_Nt*2));
  } else {
    Svec[0] = biest::Surface<Real>(NFP*surf_Nt, surf_Np, biest::SurfType::None);
    biest::SurfaceOp<Real>::CompleteVecField(Svec[0].Coord(), true, half_period, NFP, surf_Nt, surf_Np, sctl::Vector<Real>(X), (Real)0);
  }

  NFP_ = NFP;
  half_period_ = half_period;
  src_Nt_ = src_Nt;
  src_Np_ = src_Np;
  trg_Nt_ = trg_Nt;
  trg_Np_ = trg_Np;
  digits_ = digits;
}


template <class Real> std::vector<Real> VirtualCasing<Real>::ComputeBext(const std::vector<Real>& B0) const {
  return ComputeB(B0, true);
}

template <class Real> std::vector<Real> VirtualCasing<Real>::ComputeBextOffSurf(const std::vector<Real>& B0, const std::vector<Real>& Xt, const sctl::Long max_Nt, const sctl::Long max_Np) const {
  return ComputeBOffSurf(B0, Xt, max_Nt, max_Np, true);
}

template <class Real> std::vector<Real> VirtualCasing<Real>::ComputeGradBext(const std::vector<Real>& B0) const {
  return ComputeGradB(B0, true);
}


template <class Real> std::vector<Real> VirtualCasing<Real>::ComputeBint(const std::vector<Real>& B0) const {
  return ComputeB(B0, false);
}

template <class Real> std::vector<Real> VirtualCasing<Real>::ComputeBintOffSurf(const std::vector<Real>& B0, const std::vector<Real>& Xt, const sctl::Long max_Nt, const sctl::Long max_Np) const {
  return ComputeBOffSurf(B0, Xt, max_Nt, max_Np, false);
}

//template <class Real> std::vector<Real> VirtualCasing<Real>::ComputeGradBint(const std::vector<Real>& B0) const {
//  return ComputeGradB(B0, false);
//}


template <class Real> std::vector<Real> VirtualCasing<Real>::ComputeB(const std::vector<Real>& B0, bool ext) const {
  SCTL_ASSERT((sctl::Long)B0.size() == COORD_DIM * src_Nt_ * src_Np_);
  SCTL_ASSERT(trg_Nt_ > 0);
  SCTL_ASSERT(trg_Np_ > 0);

  if (dosetup) {
    //BiotSavartFxU.SetupSingular(Svec, biest::BiotSavart3D<Real>::FxU(), digits_, NFP_*(half_period_?2:1), NFP_*(half_period_?2:1)*src_Nt_, src_Np_, trg_Nt_, trg_Np_);
    LaplaceFxdU.SetupSingular(Svec, biest::Laplace3D<Real>::FxdU(), digits_, NFP_*(half_period_?2:1), NFP_*(half_period_?2:1)*src_Nt_, src_Np_, trg_Nt_, trg_Np_);
    quad_Nt_ = LaplaceFxdU.QuadNt();
    quad_Np_ = LaplaceFxdU.QuadNp();

    sctl::Vector<Real> XX;
    biest::SurfaceOp<Real>::Resample(XX, quad_Nt_, quad_Np_, Svec[0].Coord(), Svec[0].NTor(), Svec[0].NPol());

    biest::SurfaceOp<Real> SurfOp(comm_, quad_Nt_, quad_Np_);
    SurfOp.Grad2D(dX, XX);
    SurfOp.SurfNormalAreaElem(&normal, nullptr, dX, &XX);

    dosetup = false;
  }

  sctl::Vector<Real> B0_, B;
  biest::SurfaceOp<Real>::CompleteVecField(B0_, false, half_period_, NFP_, src_Nt_, src_Np_, sctl::Vector<Real>(B0), (half_period_?sctl::const_pi<Real>()*(1/(Real)(NFP_*trg_Nt_*2)-1/(Real)(NFP_*src_Nt_*2)):0));
  biest::SurfaceOp<Real>::Resample(B, quad_Nt_, quad_Np_, B0_, NFP_*(half_period_?2:1)*src_Nt_, src_Np_);

  std::vector<Real> Bvc;
  Real sign = (ext ? 1 : -1); // flip sign for internal currents
  { // Bvc = BiotSavartFxU.Eval(normal x B) * sign;
    sctl::Vector<Real> J;
    CrossProd(J, normal, B);
    if (0) {
      sctl::Vector<Real> Bvc_;
      BiotSavartFxU.Eval(Bvc_, -J * sign);
      Bvc.assign(Bvc_.begin(), Bvc_.end());
    } else {
      const sctl::Long N = trg_Nt_ * trg_Np_;
      sctl::Vector<Real> gradG_J(N * COORD_DIM * COORD_DIM); gradG_J = 0;
      sctl::Vector<Real> gradG_J0(N*COORD_DIM, gradG_J.begin() + N*COORD_DIM*0, false);
      sctl::Vector<Real> gradG_J1(N*COORD_DIM, gradG_J.begin() + N*COORD_DIM*1, false);
      sctl::Vector<Real> gradG_J2(N*COORD_DIM, gradG_J.begin() + N*COORD_DIM*2, false);
      LaplaceFxdU.Eval(gradG_J0, sctl::Vector<Real>(J.Dim()/3, J.begin() + (J.Dim()/3)*0, false));
      LaplaceFxdU.Eval(gradG_J1, sctl::Vector<Real>(J.Dim()/3, J.begin() + (J.Dim()/3)*1, false));
      LaplaceFxdU.Eval(gradG_J2, sctl::Vector<Real>(J.Dim()/3, J.begin() + (J.Dim()/3)*2, false));

      if ((sctl::Long)Bvc.size() != N * COORD_DIM) Bvc.resize(N * COORD_DIM);
      for (sctl::Long i = 0; i < N; i++) {
        for (sctl::Integer k = 0; k < COORD_DIM; k++) {
          const sctl::Integer k1 = (k+1)%COORD_DIM;
          const sctl::Integer k2 = (k+2)%COORD_DIM;
          Bvc[k*N+i] = (gradG_J[(k1*COORD_DIM+k2)*N+i] - gradG_J[(k2*COORD_DIM+k1)*N+i]) * sign;
        }
      }
    }
  }
  { // Bvc += gradG[B . normal] * sign + B/2
    sctl::Vector<Real> B_;
    const sctl::Long trg_Nt__ = (half_period_?2:1)*trg_Nt_;
    biest::SurfaceOp<Real>::Resample(B_, NFP_*trg_Nt__, trg_Np_, B0_, NFP_*(half_period_?2:1)*src_Nt_, src_Np_);

    sctl::Vector<Real> BdotN, Bvc_;
    DotProd(BdotN, B, normal);
    LaplaceFxdU.Eval(Bvc_, BdotN);
    for (sctl::Long k = 0; k < COORD_DIM; k++) {
      for (sctl::Long i = 0; i < trg_Nt_; i++) {
        for (sctl::Long j = 0; j < trg_Np_; j++) {
          Bvc[(k*trg_Nt_+i)*trg_Np_+j] += Bvc_[(k*trg_Nt_+i)*trg_Np_+j] * sign + (Real)0.5 * B_[(k*NFP_*trg_Nt__+i)*trg_Np_+j];
        }
      }
    }
  }
  return Bvc;
}

template <class Real> std::vector<Real> VirtualCasing<Real>::ComputeBOffSurf(const std::vector<Real>& B0, const std::vector<Real>& Xt, const sctl::Long max_Nt, const sctl::Long max_Np, bool ext) const {
  SCTL_ASSERT((sctl::Long)B0.size() == COORD_DIM * src_Nt_ * src_Np_);
  const sctl::Long quad_Nt_ = std::max(NFP_*(half_period_?2:1)*src_Nt_, Svec[0].NTor());
  const sctl::Long quad_Np_ = std::max(src_Np_, Svec[0].NPol());

  sctl::Vector<Real> XX;
  biest::SurfaceOp<Real>::Resample(XX, quad_Nt_, quad_Np_, Svec[0].Coord(), Svec[0].NTor(), Svec[0].NPol());
  sctl::Vector<biest::Surface<Real>> Svec_(1);
  Svec_[0] = biest::Surface<Real>(quad_Nt_, quad_Np_);
  Svec_[0].Coord() = XX;

  biest::SurfaceOp<Real> SurfOp(comm_, quad_Nt_, quad_Np_);
  sctl::Vector<Real> dX, normal;
  SurfOp.Grad2D(dX, XX);
  SurfOp.SurfNormalAreaElem(&normal, nullptr, dX, &XX);

  sctl::Vector<Real> B0_, B;
  biest::SurfaceOp<Real>::CompleteVecField(B0_, false, half_period_, NFP_, src_Nt_, src_Np_, sctl::Vector<Real>(B0), (half_period_?sctl::const_pi<Real>()*(1/(Real)(NFP_*trg_Nt_*2)-1/(Real)(NFP_*src_Nt_*2)):0));
  biest::SurfaceOp<Real>::Resample(B, quad_Nt_, quad_Np_, B0_, NFP_*(half_period_?2:1)*src_Nt_, src_Np_);

  sctl::Vector<Real> BdotN, J;
  DotProd(BdotN, B, normal);
  CrossProd(J, normal, B);

  std::vector<Real> Bvc;
  Real sign = (ext ? 1 : -1); // flip sign for internal currents
  if (0) { // Bvc = BiotSavartFxU.Eval(normal x B) * sign;
    SCTL_ASSERT(quad_Nt_ >= 13);
    SCTL_ASSERT(quad_Np_ >= 13);
    biest::BoundaryIntegralOp<Real, 3, 3, 1, 6, 1> BiotSavartFxU;
    BiotSavartFxU.SetupSingular(Svec_, biest::BiotSavart3D<Real>::FxU());

    sctl::Vector<Real> Bvc_;
    BiotSavartFxU.EvalOffSurface(Bvc_, sctl::Vector<Real>(Xt), J * (-sign));
    Bvc.assign(Bvc_.begin(), Bvc_.end());
  }
  if (0) { // Bvc += gradG[B . normal]
    sctl::Vector<Real> Bvc_;
    biest::BoundaryIntegralOp<Real, 1, 3, 1, 6, 1> LaplaceFxdU;
    LaplaceFxdU.SetupSingular(Svec_, biest::Laplace3D<Real>::FxdU());
    LaplaceFxdU.EvalOffSurface(Bvc_, sctl::Vector<Real>(Xt), BdotN);
    for (sctl::Long i = 0; i < Bvc_.Dim(); i++) Bvc[i] += Bvc_[i] * sign;
  }
  { // Adaptive evaluation
    const auto& X = Svec[0].Coord();
    std::vector<Real> XX(X.begin(), X.end());
    std::vector<Real> J_(J.begin(), J.end());
    for (auto& j : J_) j = j * (-sign);

    std::vector<Real> BdotN_(BdotN.Dim());
    for (sctl::Long i = 0; i < BdotN.Dim(); i++) BdotN_[i] = BdotN[i] * (-sign);

    biest::ExtVacuumField<Real> ext_vacuum;
    ext_vacuum.Setup(digits_, 1, Svec[0].NTor(), Svec[0].NPol(), XX, quad_Nt_, quad_Np_);
    Bvc = ext_vacuum.EvalOffSurface(Xt, BdotN_, J_, (half_period_?2:1)*max_Nt, max_Np);
  }
  return Bvc;
}

template <class Real> std::vector<Real> VirtualCasing<Real>::ComputeGradB(const std::vector<Real>& B0, bool ext) const {
  SCTL_ASSERT((sctl::Long)B0.size() == COORD_DIM * src_Nt_ * src_Np_);

  if (dosetup_grad) {
    //BiotSavartFxdU.SetupSingular(Svec, biest::BiotSavart3D<Real>::FxdU(), digits_, NFP_*(half_period_?2:1), NFP_*(half_period_?2:1)*src_Nt_, src_Np_, (half_period_?2:1)*trg_Nt_, trg_Np_);
    LaplaceFxd2U.SetupSingular(Svec, biest::Laplace3D<Real>::Fxd2U(), digits_, NFP_*(half_period_?2:1), NFP_*(half_period_?2:1)*src_Nt_, src_Np_, trg_Nt_, trg_Np_);
    grad_quad_Nt_ = LaplaceFxd2U.QuadNt();
    grad_quad_Np_ = LaplaceFxd2U.QuadNp();

    sctl::Vector<Real> XX;
    biest::SurfaceOp<Real>::Resample(XX, grad_quad_Nt_, grad_quad_Np_, Svec[0].Coord(), Svec[0].NTor(), Svec[0].NPol());

    biest::SurfaceOp<Real> SurfOp(comm_, grad_quad_Nt_, grad_quad_Np_);
    SurfOp.Grad2D(dX, XX);
    SurfOp.SurfNormalAreaElem(&normal, nullptr, dX, &XX);

    dosetup_grad = false;
  }

  sctl::Vector<Real> B0_, B;
  biest::SurfaceOp<Real>::CompleteVecField(B0_, false, half_period_, NFP_, src_Nt_, src_Np_, sctl::Vector<Real>(B0), (half_period_?sctl::const_pi<Real>()*(1/(Real)(NFP_*trg_Nt_*2)-1/(Real)(NFP_*src_Nt_*2)):0));
  biest::SurfaceOp<Real>::Resample(B, grad_quad_Nt_, grad_quad_Np_, B0_, NFP_*(half_period_?2:1)*src_Nt_, src_Np_);

  std::vector<Real> gradBvc;
  Real sign = (ext ? 1 : -1); // flip sign for internal currents
  SCTL_ASSERT(ext = true); // Hedgehog quadrature currently does not allow flipping the normal vector
  { // gradBvc = gradBiotSavartFxU.Eval(normal x B) * sign;
    sctl::Vector<Real> J;
    CrossProd(J, normal, B);
    if (0) {
      sctl::Vector<Real> gradBvc_;
      BiotSavartFxdU.Eval(gradBvc_, -J * sign);
      gradBvc.assign(gradBvc_.begin(), gradBvc_.end());
    } else {
      const sctl::Long N = trg_Nt_ * trg_Np_;
      sctl::Vector<Real> gradG_J(N * COORD_DIM * COORD_DIM * COORD_DIM); gradG_J = 0;
      sctl::Vector<Real> gradG_J0(N*COORD_DIM*COORD_DIM, gradG_J.begin() + N*COORD_DIM*COORD_DIM*0, false);
      sctl::Vector<Real> gradG_J1(N*COORD_DIM*COORD_DIM, gradG_J.begin() + N*COORD_DIM*COORD_DIM*1, false);
      sctl::Vector<Real> gradG_J2(N*COORD_DIM*COORD_DIM, gradG_J.begin() + N*COORD_DIM*COORD_DIM*2, false);
      LaplaceFxd2U.Eval(gradG_J0, sctl::Vector<Real>(J.Dim()/3, J.begin() + (J.Dim()/3)*0, false));
      LaplaceFxd2U.Eval(gradG_J1, sctl::Vector<Real>(J.Dim()/3, J.begin() + (J.Dim()/3)*1, false));
      LaplaceFxd2U.Eval(gradG_J2, sctl::Vector<Real>(J.Dim()/3, J.begin() + (J.Dim()/3)*2, false));

      if ((sctl::Long)gradBvc.size() != COORD_DIM * COORD_DIM * N) gradBvc.resize(COORD_DIM * COORD_DIM * N);
      for (sctl::Long i = 0; i < COORD_DIM*N; i++) {
        for (sctl::Integer k = 0; k < COORD_DIM; k++) {
          const sctl::Integer k1 = (k+1)%COORD_DIM;
          const sctl::Integer k2 = (k+2)%COORD_DIM;
          gradBvc[k*COORD_DIM*N+i] = (gradG_J[(k1*COORD_DIM+k2)*COORD_DIM*N+i] - gradG_J[(k2*COORD_DIM+k1)*COORD_DIM*N+i]) * sign;
        }
      }
    }
  }
  { // gradBvc += gradgradG[B . normal] * sign
    sctl::Vector<Real> B_;
    const sctl::Long trg_Nt__ = (half_period_?2:1)*trg_Nt_;
    biest::SurfaceOp<Real>::Resample(B_, NFP_*trg_Nt__, trg_Np_, B0_, NFP_*(half_period_?2:1)*src_Nt_, src_Np_);

    sctl::Vector<Real> BdotN, gradBvc_;
    DotProd(BdotN, B, normal);
    LaplaceFxd2U.Eval(gradBvc_, BdotN);
    for (sctl::Long k = 0; k < COORD_DIM*COORD_DIM; k++) {
      for (sctl::Long i = 0; i < trg_Nt_; i++) {
        for (sctl::Long j = 0; j < trg_Np_; j++) {
          gradBvc[(k*trg_Nt_+i)*trg_Np_+j] += gradBvc_[(k*trg_Nt_+i)*trg_Np_+j] * sign;
        }
      }
    }
  }
  return gradBvc;
}


template <class Real> std::vector<Real> VirtualCasing<Real>::GetNormal(const sctl::Integer NFP, const bool half_period, const sctl::Long Nt, const sctl::Long Np) const {
  SCTL_ASSERT(Svec[0].NTor() && Svec[0].NPol());
  const sctl::Long Nt_ = NFP*(half_period?2:1)*Nt;
  const sctl::Long skip_tor = (sctl::Long)sctl::ceil(Svec[0].NTor()/(Real)Nt_);
  const sctl::Long skip_pol = (sctl::Long)sctl::ceil(Svec[0].NPol()/(Real)Np);
  const sctl::Long Nt0 = skip_tor*Nt_;
  const sctl::Long Np0 = skip_pol*Np;

  sctl::Vector<Real> X1, X2, dX, normal;
  biest::SurfaceOp<Real>::RotateToroidal(X1, Svec[0].Coord(), Svec[0].NTor(), Svec[0].NPol(), (half_period?sctl::const_pi<Real>()/(NFP*Nt*2):0)-(half_period_?sctl::const_pi<Real>()/(NFP_*trg_Nt_*2):0));
  biest::SurfaceOp<Real>::Resample(X2, Nt0, Np0, X1, Svec[0].NTor(), Svec[0].NPol());
  biest::SurfaceOp<Real> SurfOp(comm_, Nt0, Np0);
  SurfOp.Grad2D(dX, X2);
  SurfOp.SurfNormalAreaElem(&normal, nullptr, dX, &X2);

  std::vector<Real> Xn(COORD_DIM*Nt*Np);
  for (sctl::Integer k = 0; k < COORD_DIM; k++) {
    for (sctl::Long t = 0; t < Nt; t++) {
      for (sctl::Long p = 0; p < Np; p++) {
        Xn[(k*Nt+t)*Np+p] = normal[(k*Nt0+t*skip_tor)*Np0+p*skip_pol];
      }
    }
  }
  return Xn;
}

template <class Real> void VirtualCasing<Real>::DotProd(sctl::Vector<Real>& AdotB, const sctl::Vector<Real>& A, const sctl::Vector<Real>& B) {
  sctl::Long N = A.Dim() / COORD_DIM;
  SCTL_ASSERT(A.Dim() == COORD_DIM * N);
  SCTL_ASSERT(B.Dim() == COORD_DIM * N);
  if (AdotB.Dim() != N) AdotB.ReInit(N);
  for (sctl::Long i = 0; i < N; i++) {
    Real AdotB_ = 0;
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      AdotB_ += A[k*N+i] * B[k*N+i];
    }
    AdotB[i] = AdotB_;
  }
};

template <class Real> void VirtualCasing<Real>::CrossProd(sctl::Vector<Real>& AcrossB, const sctl::Vector<Real>& A, const sctl::Vector<Real>& B) {
  sctl::Long N = A.Dim() / COORD_DIM;
  SCTL_ASSERT(A.Dim() == COORD_DIM * N);
  SCTL_ASSERT(B.Dim() == COORD_DIM * N);
  if (AcrossB.Dim() != COORD_DIM * N) AcrossB.ReInit(COORD_DIM * N);
  for (sctl::Long i = 0; i < N; i++) {
    sctl::StaticArray<Real,COORD_DIM> A_, B_, AcrossB_;
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      A_[k] = A[k*N+i];
      B_[k] = B[k*N+i];
    }
    AcrossB_[0] = A_[1] * B_[2] - B_[1] * A_[2];
    AcrossB_[1] = A_[2] * B_[0] - B_[2] * A_[0];
    AcrossB_[2] = A_[0] * B_[1] - B_[0] * A_[1];
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      AcrossB[k*N+i] = AcrossB_[k];
    }
  }
};

template <class Real> std::vector<Real> VirtualCasingTestData<Real>::SurfaceCoordinates(const sctl::Integer NFP, const bool half_period, const sctl::Long Nt, const sctl::Long Np, const biest::SurfType surf_type) {
  sctl::Vector<Real> X_;
  const sctl::Long Nt_ = (half_period?2:1)*Nt;
  biest::Surface<Real> S(NFP*Nt_, Np, surf_type);
  biest::SurfaceOp<Real>::RotateToroidal(X_, S.Coord(), NFP*Nt_, Np, (half_period?sctl::const_pi<Real>()/(NFP*Nt*2):0));

  std::vector<Real> X(COORD_DIM*Nt*Np);
  for (sctl::Long k = 0; k < COORD_DIM; k++) {
    for (sctl::Long i = 0; i < Nt*Np; i++) {
      X[k*Nt*Np+i] = X_[k*NFP*Nt_*Np+i];
    }
  }
  return X;
}

template <class Real> std::tuple<std::vector<Real>, std::vector<Real>> VirtualCasingTestData<Real>::BFieldDataOffSurf(const sctl::Integer NFP, const bool half_period, const sctl::Long surf_Nt, const sctl::Long surf_Np, const std::vector<Real>& X, const std::vector<Real>& X_trg) {
  auto WriteVTK_ = [](const std::string& fname, const sctl::Vector<sctl::Vector<Real>>& coords, const sctl::Vector<sctl::Vector<Real>>& values) {
    biest::VTKData data;
    typedef biest::VTKData::VTKReal VTKReal;
    auto& point_coord =data.point_coord ;
    auto& point_value =data.point_value ;
    auto& line_connect=data.line_connect;
    auto& line_offset =data.line_offset ;
    constexpr sctl::Integer COORD_DIM = biest::VTKData::COORD_DIM;

    SCTL_ASSERT(coords.Dim() == values.Dim());
    for (sctl::Long j = 0; j < coords.Dim(); j++) { // set point_coord, line_connect
      const auto& coord = coords[j];
      const auto& value = values[j];
      sctl::Long N = coord.Dim() / COORD_DIM;
      sctl::Long dof = value.Dim() / N;
      SCTL_ASSERT(value.Dim() == dof * N);
      for (sctl::Long i = 0; i < N; i++) {
        line_connect.push_back(point_coord.size()/COORD_DIM);
        point_coord.push_back((VTKReal)coord[0*N+i]);
        point_coord.push_back((VTKReal)coord[1*N+i]);
        point_coord.push_back((VTKReal)coord[2*N+i]);
        for (sctl::Long k = 0; k < dof; k++) {
          point_value.push_back((VTKReal)value[k*N+i]);
        }
      }
      line_offset.push_back(line_connect.size());
    }
    data.WriteVTK(fname.c_str(), sctl::Comm::Self());
  };
  auto eval_BiotSavart = [](const sctl::Vector<Real>& Xt, const sctl::Vector<sctl::Vector<Real>>& source, const sctl::Vector<sctl::Vector<Real>>& density) {
    const auto& kernel = biest::BiotSavart3D<Real>::FxU();
    sctl::Long Nt = Xt.Dim() / COORD_DIM;
    SCTL_ASSERT(Xt.Dim() == COORD_DIM * Nt);
    SCTL_ASSERT(source.Dim() == density.Dim());

    sctl::Vector<Real> B(COORD_DIM*Nt);
    B = 0;
    for (sctl::Long i = 0; i < source.Dim(); i++) {
      const auto& Xs = source[i];
      const auto& Fs = density[i];
      sctl::Long Ns = Xs.Dim() / COORD_DIM;
      SCTL_ASSERT(Xs.Dim() == COORD_DIM * Ns);
      SCTL_ASSERT(Fs.Dim() == COORD_DIM * Ns);
      sctl::Vector<Real> SrcNormal(COORD_DIM*Ns);
      kernel(Xs,SrcNormal,Fs, Xt,B);
    }
    return B;
  };
  auto add_source_loop = [](sctl::Vector<sctl::Vector<Real>>& source, sctl::Vector<sctl::Vector<Real>>& density, const std::initializer_list<Real> coord, const std::initializer_list<Real> normal, const Real radius) {
    auto cross_norm = [](const sctl::Vector<Real>& A, const sctl::Vector<Real>& B) {
      sctl::Vector<Real> C(COORD_DIM);
      C[0] = A[1]*B[2] - B[1]*A[2];
      C[1] = A[2]*B[0] - B[2]*A[0];
      C[2] = A[0]*B[1] - B[0]*A[1];
      Real r = sctl::sqrt<Real>(C[0]*C[0]+C[1]*C[1]+C[2]*C[2]);
      return C*(1/r);
    };
    sctl::Vector<Real> coord_(COORD_DIM), normal_(COORD_DIM), e0(COORD_DIM), e1(COORD_DIM);
    coord_[0] = coord.begin()[0];
    coord_[1] = coord.begin()[1];
    coord_[2] = coord.begin()[2];
    normal_[0] = normal.begin()[0];
    normal_[1] = normal.begin()[1];
    normal_[2] = normal.begin()[2];
    Real normal_scal = 1/sctl::sqrt<Real>(normal_[0]*normal_[0]+normal_[1]*normal_[1]+normal_[2]*normal_[2]);
    normal_ *= normal_scal;

    e0[0] = (Real)drand48();
    e0[1] = (Real)drand48();
    e0[2] = (Real)drand48();
    e0 = cross_norm(e0,normal_)*radius;
    e1 = cross_norm(e0,normal_)*radius;

    sctl::Long N = 10000;
    sctl::Vector<Real> X(COORD_DIM * N);
    sctl::Vector<Real> dX(COORD_DIM * N);
    for (sctl::Long i = 0; i < N; i++) {
      Real t = 2*sctl::const_pi<Real>()*i/N;
      sctl::Vector<Real> r = coord_ + e0*sctl::sin<Real>(t) + e1*sctl::cos<Real>(t);
      sctl::Vector<Real> dr = e0*sctl::cos<Real>(t) - e1*sctl::sin<Real>(t);
      X[0*N+i] = r[0];
      X[1*N+i] = r[1];
      X[2*N+i] = r[2];
      dX[0*N+i] = dr[0];
      dX[1*N+i] = dr[1];
      dX[2*N+i] = dr[2];
    }
    source.PushBack(X);
    density.PushBack(dX);
  };

  sctl::Comm comm = sctl::Comm::Self();
  sctl::Vector<Real> X_surf;
  { // Set X_surf
    sctl::Vector<Real> XX;
    biest::SurfaceOp<Real>::CompleteVecField(XX, true, half_period, NFP, surf_Nt, surf_Np, sctl::Vector<Real>(X), (half_period?-sctl::const_pi<Real>()/(NFP*surf_Nt*2):0));
    biest::SurfaceOp<Real>::Resample(X_surf, NFP*(half_period?2:1)*(surf_Nt+1), surf_Np, XX, NFP*(half_period?2:1)*surf_Nt, surf_Np);
  }

  sctl::Vector<sctl::Vector<Real>> source0, density0;
  sctl::Vector<sctl::Vector<Real>> source1, density1;
  { // Set inside sources (source0, density0)
    sctl::Long N = 20000;
    sctl::Vector<Real> X(COORD_DIM*N), dX(COORD_DIM*N);
    { // Set X, dX
      sctl::Long Nt = 100, Np = 100;
      sctl::Vector<Real> coord(COORD_DIM*Nt);
      { // Set coord
        auto S = biest::Surface<Real>(Nt,Np, biest::SurfType::None);
        biest::SurfaceOp<Real>::Upsample(X_surf, NFP*(half_period?2:1)*(surf_Nt+1), surf_Np, S.Coord(), Nt, Np);

        sctl::Vector<Real> normal, dX;
        biest::SurfaceOp<Real> SurfOp(comm, Nt, Np);
        SurfOp.Grad2D(dX, S.Coord());
        SurfOp.SurfNormalAreaElem(&normal, nullptr, dX, &S.Coord());
        S.Coord() += (Real)(-2.17)*normal;

        coord = 0;
        for (sctl::Long t = 0; t < Nt; t++) {
          for (sctl::Long p = 0; p < Np; p++) {
            coord[0*Nt+t] += S.Coord()[(0*Nt+t)*Np+p]/Np;
            coord[1*Nt+t] += S.Coord()[(1*Nt+t)*Np+p]/Np;
            coord[2*Nt+t] += S.Coord()[(2*Nt+t)*Np+p]/Np;
          }
        }
      }

      sctl::Vector<Real> dX_;
      biest::SurfaceOp<Real>::Upsample(coord,Nt,1, X,N,1);
      biest::SurfaceOp<Real> SurfOp(comm,N,1);
      SurfOp.Grad2D(dX_, X);
      for (sctl::Long i = 0; i < N; i++) {
        for (sctl::Integer k = 0; k < COORD_DIM; k++) {
          dX[k*N+i] = dX_[(2*k+0)*N+i];
        }
      }
    }
    source0.PushBack(X);
    density0.PushBack(dX*(Real)0.05);
  }
  for (int i = 0; i < NFP; i++) { // Set outside sources (source1, density1)
    const sctl::Long N = source0[0].Dim()/COORD_DIM;
    const sctl::Long Nskip = i * N / NFP;

    const sctl::StaticArray<Real,COORD_DIM> X{source0[0][0*N+Nskip],source0[0][1*N+Nskip],source0[0][2*N+Nskip]};
    const sctl::StaticArray<Real,COORD_DIM> Xn{density0[0][0*N+Nskip],density0[0][1*N+Nskip],density0[0][2*N+Nskip]};
    const Real R = sctl::sqrt<Real>(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
    add_source_loop(source1, density1, {X[0],X[1],X[2]}, {Xn[0],Xn[1],Xn[2]}, R);
  }

  const auto Bint_ = eval_BiotSavart(sctl::Vector<Real>(X_trg), source0, density0);
  const auto Bext_ = eval_BiotSavart(sctl::Vector<Real>(X_trg), source1, density1);

  const sctl::Long Ntrg = X_trg.size() / COORD_DIM;
  std::vector<Real> Bint(3*Ntrg);
  std::vector<Real> Bext(3*Ntrg);
  for (sctl::Integer k = 0; k < COORD_DIM; k++) { // Set Bint, Bext
    for (sctl::Long i = 0; i < Ntrg; i++) {
      Bint[k*Ntrg+i] = Bint_[k*Ntrg+i];
      Bext[k*Ntrg+i] = Bext_[k*Ntrg+i];
    }
  }

  if (0) { // Visualization
    WriteVTK_("loop0", source0, density0);
    WriteVTK_("loop1", source1, density1);
  }

  SCTL_UNUSED(WriteVTK_);

  return std::make_tuple(std::move(Bext), std::move(Bint));
}

template <class Real> std::tuple<std::vector<Real>, std::vector<Real>> VirtualCasingTestData<Real>::BFieldData(const sctl::Integer NFP, const bool half_period, const sctl::Long surf_Nt, const sctl::Long surf_Np, const std::vector<Real>& X, const sctl::Long trg_Nt, const sctl::Long trg_Np) {
  std::vector<Real> X_trg;
  { // Set X_trg
    sctl::Vector<Real> XX, X_surf;
    biest::SurfaceOp<Real>::CompleteVecField(XX, true, half_period, NFP, surf_Nt, surf_Np, sctl::Vector<Real>(X), (half_period?-sctl::const_pi<Real>()/(NFP*surf_Nt*2):0));
    biest::SurfaceOp<Real>::Resample(X_surf, NFP*(half_period?2:1)*(surf_Nt+1), surf_Np, XX, NFP*(half_period?2:1)*surf_Nt, surf_Np);

    sctl::Vector<Real> X_surf_shifted, X_trg_;
    const sctl::Long trg_Nt_ = (half_period?2:1)*trg_Nt;
    biest::SurfaceOp<Real>::RotateToroidal(X_surf_shifted, X_surf, NFP*(half_period?2:1)*(surf_Nt+1), surf_Np, (half_period?sctl::const_pi<Real>()/(NFP*trg_Nt*2):0));
    biest::SurfaceOp<Real>::Resample(X_trg_, NFP*trg_Nt_, trg_Np, X_surf_shifted, NFP*(half_period?2:1)*(surf_Nt+1), surf_Np);
    X_trg.resize(COORD_DIM*trg_Nt_*trg_Np);
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      for (sctl::Long i = 0; i < trg_Nt_*trg_Np; i++) {
        X_trg[k*trg_Nt_*trg_Np+i] = X_trg_[k*NFP*trg_Nt_*trg_Np+i];
      }
    }
  }

  std::vector<Real> Bint_, Bext_;
  std::tie(Bext_, Bint_) = BFieldDataOffSurf(NFP, half_period, surf_Nt, surf_Np, X, X_trg);

  std::vector<Real> Bint(COORD_DIM*trg_Nt*trg_Np);
  std::vector<Real> Bext(COORD_DIM*trg_Nt*trg_Np);
  for (sctl::Integer k = 0; k < COORD_DIM; k++) { // Set Bint, Bext
    const sctl::Long trg_Nt_ = (half_period?2:1)*trg_Nt;
    for (sctl::Long i = 0; i < trg_Nt*trg_Np; i++) {
      Bint[k*trg_Nt*trg_Np+i] = Bint_[k*trg_Nt_*trg_Np+i];
      Bext[k*trg_Nt*trg_Np+i] = Bext_[k*trg_Nt_*trg_Np+i];
    }
  }

  if (0) { // Visualization
    biest::WriteVTK("B", NFP, half_period, surf_Nt, surf_Np, sctl::Vector<Real>(X), trg_Nt, trg_Np, sctl::Vector<Real>(Bext)+sctl::Vector<Real>(Bint));
  }

  return std::make_tuple(std::move(Bext), std::move(Bint));
}

template <class Real> std::tuple<std::vector<Real>, std::vector<Real>> VirtualCasingTestData<Real>::GradBFieldData(const sctl::Integer NFP, const bool half_period, const sctl::Long surf_Nt, const sctl::Long surf_Np, const std::vector<Real>& X, const sctl::Long trg_Nt, const sctl::Long trg_Np) {
  auto WriteVTK_ = [](const std::string& fname, const sctl::Vector<sctl::Vector<Real>>& coords, const sctl::Vector<sctl::Vector<Real>>& values) {
    biest::VTKData data;
    typedef biest::VTKData::VTKReal VTKReal;
    auto& point_coord =data.point_coord ;
    auto& point_value =data.point_value ;
    auto& line_connect=data.line_connect;
    auto& line_offset =data.line_offset ;
    constexpr sctl::Integer COORD_DIM = biest::VTKData::COORD_DIM;

    SCTL_ASSERT(coords.Dim() == values.Dim());
    for (sctl::Long j = 0; j < coords.Dim(); j++) { // set point_coord, line_connect
      const auto& coord = coords[j];
      const auto& value = values[j];
      sctl::Long N = coord.Dim() / COORD_DIM;
      sctl::Long dof = value.Dim() / N;
      SCTL_ASSERT(value.Dim() == dof * N);
      for (sctl::Long i = 0; i < N; i++) {
        line_connect.push_back(point_coord.size()/COORD_DIM);
        point_coord.push_back((VTKReal)coord[0*N+i]);
        point_coord.push_back((VTKReal)coord[1*N+i]);
        point_coord.push_back((VTKReal)coord[2*N+i]);
        for (sctl::Long k = 0; k < dof; k++) {
          point_value.push_back((VTKReal)value[k*N+i]);
        }
      }
      line_offset.push_back(line_connect.size());
    }
    data.WriteVTK(fname.c_str(), sctl::Comm::Self());
  };
  auto eval_BiotSavartGrad = [](const sctl::Vector<Real>& Xt, const sctl::Vector<sctl::Vector<Real>>& source, const sctl::Vector<sctl::Vector<Real>>& density) {
    const auto& kernel = biest::BiotSavart3D<Real>::FxdU();
    sctl::Long Nt = Xt.Dim() / COORD_DIM;
    SCTL_ASSERT(Xt.Dim() == COORD_DIM * Nt);
    SCTL_ASSERT(source.Dim() == density.Dim());

    sctl::Vector<Real> B(COORD_DIM*COORD_DIM*Nt);
    B = 0;
    for (sctl::Long i = 0; i < source.Dim(); i++) {
      const auto& Xs = source[i];
      const auto& Fs = density[i];
      sctl::Long Ns = Xs.Dim() / COORD_DIM;
      SCTL_ASSERT(Xs.Dim() == COORD_DIM * Ns);
      SCTL_ASSERT(Fs.Dim() == COORD_DIM * Ns);
      sctl::Vector<Real> SrcNormal(COORD_DIM*Ns);
      kernel(Xs,SrcNormal,Fs, Xt,B);
    }
    return B;
  };
  auto add_source_loop = [](sctl::Vector<sctl::Vector<Real>>& source, sctl::Vector<sctl::Vector<Real>>& density, const std::initializer_list<Real> coord, const std::initializer_list<Real> normal, const Real radius) {
    auto cross_norm = [](const sctl::Vector<Real>& A, const sctl::Vector<Real>& B) {
      sctl::Vector<Real> C(COORD_DIM);
      C[0] = A[1]*B[2] - B[1]*A[2];
      C[1] = A[2]*B[0] - B[2]*A[0];
      C[2] = A[0]*B[1] - B[0]*A[1];
      Real r = sctl::sqrt<Real>(C[0]*C[0]+C[1]*C[1]+C[2]*C[2]);
      return C*(1/r);
    };
    sctl::Vector<Real> coord_(COORD_DIM), normal_(COORD_DIM), e0(COORD_DIM), e1(COORD_DIM);
    coord_[0] = coord.begin()[0];
    coord_[1] = coord.begin()[1];
    coord_[2] = coord.begin()[2];
    normal_[0] = normal.begin()[0];
    normal_[1] = normal.begin()[1];
    normal_[2] = normal.begin()[2];
    Real normal_scal = 1/sctl::sqrt<Real>(normal_[0]*normal_[0]+normal_[1]*normal_[1]+normal_[2]*normal_[2]);
    normal_ *= normal_scal;

    e0[0] = drand48();
    e0[1] = drand48();
    e0[2] = drand48();
    e0 = cross_norm(e0,normal_)*radius;
    e1 = cross_norm(e0,normal_)*radius;

    sctl::Long N = 10000;
    sctl::Vector<Real> X(COORD_DIM * N);
    sctl::Vector<Real> dX(COORD_DIM * N);
    for (sctl::Long i = 0; i < N; i++) {
      Real t = 2*sctl::const_pi<Real>()*i/N;
      sctl::Vector<Real> r = coord_ + e0*sctl::sin<Real>(t) + e1*sctl::cos<Real>(t);
      sctl::Vector<Real> dr = e0*sctl::cos<Real>(t) - e1*sctl::sin<Real>(t);
      X[0*N+i] = r[0];
      X[1*N+i] = r[1];
      X[2*N+i] = r[2];
      dX[0*N+i] = dr[0];
      dX[1*N+i] = dr[1];
      dX[2*N+i] = dr[2];
    }
    source.PushBack(X);
    density.PushBack(dX);
  };
  auto DotProd = [](sctl::Vector<Real>& AdotB, const sctl::Vector<Real>& A, const sctl::Vector<Real>& B) {
    sctl::Long N = A.Dim() / COORD_DIM;
    SCTL_ASSERT(A.Dim() == COORD_DIM * N);
    SCTL_ASSERT(B.Dim() == COORD_DIM * N);
    if (AdotB.Dim() != N) AdotB.ReInit(N);
    for (sctl::Long i = 0; i < N; i++) {
      Real AdotB_ = 0;
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        AdotB_ += A[k*N+i] * B[k*N+i];
      }
      AdotB[i] = AdotB_;
    }
  };

  sctl::Comm comm = sctl::Comm::Self();
  sctl::Vector<Real> X_surf, X_trg;
  { // Set X_surf, X_trg
    sctl::Vector<Real> XX;
    biest::SurfaceOp<Real>::CompleteVecField(XX, true, half_period, NFP, surf_Nt, surf_Np, sctl::Vector<Real>(X), (half_period?-sctl::const_pi<Real>()/(NFP*surf_Nt*2):0));
    biest::SurfaceOp<Real>::Resample(X_surf, NFP*(half_period?2:1)*(surf_Nt+1), surf_Np, XX, NFP*(half_period?2:1)*surf_Nt, surf_Np);

    sctl::Vector<Real> X_surf_shifted, X_trg_;
    const sctl::Long trg_Nt_ = (half_period?2:1)*trg_Nt;
    biest::SurfaceOp<Real>::RotateToroidal(X_surf_shifted, X_surf, NFP*(half_period?2:1)*(surf_Nt+1), surf_Np, (half_period?sctl::const_pi<Real>()/(NFP*trg_Nt*2):0));
    biest::SurfaceOp<Real>::Resample(X_trg_, NFP*trg_Nt_, trg_Np, X_surf_shifted, NFP*(half_period?2:1)*(surf_Nt+1), surf_Np);
    X_trg.ReInit(COORD_DIM*trg_Nt_*trg_Np);
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      for (sctl::Long i = 0; i < trg_Nt_*trg_Np; i++) {
        X_trg[k*trg_Nt_*trg_Np+i] = X_trg_[k*NFP*trg_Nt_*trg_Np+i];
      }
    }
  }

  sctl::Vector<sctl::Vector<Real>> source0, density0;
  sctl::Vector<sctl::Vector<Real>> source1, density1;
  { // Set inside sources (source0, density0)
    sctl::Long N = 20000;
    sctl::Vector<Real> X(COORD_DIM*N), dX(COORD_DIM*N);
    { // Set X, dX
      sctl::Long Nt = 100, Np = 100;
      sctl::Vector<Real> coord(COORD_DIM*Nt);
      { // Set coord
        auto S = biest::Surface<Real>(Nt,Np, biest::SurfType::None);
        biest::SurfaceOp<Real>::Upsample(X_surf, NFP*(half_period?2:1)*(surf_Nt+1), surf_Np, S.Coord(), Nt, Np);

        sctl::Vector<Real> normal, dX;
        biest::SurfaceOp<Real> SurfOp(comm, Nt, Np);
        SurfOp.Grad2D(dX, S.Coord());
        SurfOp.SurfNormalAreaElem(&normal, nullptr, dX, &S.Coord());
        S.Coord() += -2.17*normal;

        coord = 0;
        for (sctl::Long t = 0; t < Nt; t++) {
          for (sctl::Long p = 0; p < Np; p++) {
            coord[0*Nt+t] += S.Coord()[(0*Nt+t)*Np+p]/Np;
            coord[1*Nt+t] += S.Coord()[(1*Nt+t)*Np+p]/Np;
            coord[2*Nt+t] += S.Coord()[(2*Nt+t)*Np+p]/Np;
          }
        }
      }

      sctl::Vector<Real> dX_;
      biest::SurfaceOp<Real>::Upsample(coord,Nt,1, X,N,1);
      biest::SurfaceOp<Real> SurfOp(comm,N,1);
      SurfOp.Grad2D(dX_, X);
      for (sctl::Long i = 0; i < N; i++) {
        for (sctl::Integer k = 0; k < COORD_DIM; k++) {
          dX[k*N+i] = dX_[(2*k+0)*N+i];
        }
      }
    }
    source0.PushBack(X);
    density0.PushBack(dX*0.05);
  }
  for (int i = 0; i < NFP; i++) { // Set outside sources (source1, density1)
    const sctl::Long N = source0[0].Dim()/COORD_DIM;
    const sctl::Long Nskip = i * N / NFP;

    const sctl::StaticArray<Real,COORD_DIM> X{source0[0][0*N+Nskip],source0[0][1*N+Nskip],source0[0][2*N+Nskip]};
    const sctl::StaticArray<Real,COORD_DIM> Xn{density0[0][0*N+Nskip],density0[0][1*N+Nskip],density0[0][2*N+Nskip]};
    const Real R = sctl::sqrt<Real>(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
    add_source_loop(source1, density1, {X[0],X[1],X[2]}, {Xn[0],Xn[1],Xn[2]}, R);
  }

  const auto GradBint_ = eval_BiotSavartGrad(X_trg, source0, density0);
  const auto GradBext_ = eval_BiotSavartGrad(X_trg, source1, density1);

  std::vector<Real> GradBint(COORD_DIM*COORD_DIM*trg_Nt*trg_Np);
  std::vector<Real> GradBext(COORD_DIM*COORD_DIM*trg_Nt*trg_Np);
  for (sctl::Integer k = 0; k < COORD_DIM*COORD_DIM; k++) { // Set GradBint, GradBext
    const sctl::Long trg_Nt_ = (half_period?2:1)*trg_Nt;
    for (sctl::Long i = 0; i < trg_Nt*trg_Np; i++) {
      GradBint[k*trg_Nt*trg_Np+i] = GradBint_[k*trg_Nt_*trg_Np+i];
      GradBext[k*trg_Nt*trg_Np+i] = GradBext_[k*trg_Nt_*trg_Np+i];
    }
  }

  if (0) { // Visualization
    biest::WriteVTK("GradB", NFP, half_period, surf_Nt, surf_Np, sctl::Vector<Real>(X), trg_Nt, trg_Np, sctl::Vector<Real>(GradBext)+sctl::Vector<Real>(GradBint));
    WriteVTK_("loop0", source0, density0);
    WriteVTK_("loop1", source1, density1);
  }

  SCTL_UNUSED(WriteVTK_);
  SCTL_UNUSED(DotProd);

  return std::make_tuple(std::move(GradBext), std::move(GradBint));
}

