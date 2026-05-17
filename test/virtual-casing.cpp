#include <virtual-casing.hpp>

template <class Real> void test(int digits, int NFP, bool half_period, long Nt, long Np, biest::SurfType surf_type, long src_Nt, long src_Np, long trg_Nt, long trg_Np) {
  // Construct the surface
  std::vector<Real> X(3*Nt*Np);
  X = VirtualCasingTestData<Real>::SurfaceCoordinates(NFP, half_period, Nt, Np, surf_type);
  //for (long t = 0; t < Nt; t++) { // toroidal direction
  //  for (long p = 0; p < Np; p++) { // poloidal direction
  //    Real x = (2 + 0.5*cos(2*M_PI*p/Np)) * cos(2*M_PI*t/(NFP*Nt));
  //    Real y = (2 + 0.5*cos(2*M_PI*p/Np)) * sin(2*M_PI*t/(NFP*Nt));
  //    Real z = 0.5*sin(2*M_PI*p/Np);
  //    X[(0*Nt+t)*Np+p] = x;
  //    X[(1*Nt+t)*Np+p] = y;
  //    X[(2*Nt+t)*Np+p] = z;
  //  }
  //}

  const auto print_relative_err = [NFP,half_period,Nt,Np,trg_Nt,trg_Np,&X](const std::string& name, const std::vector<Real>& B, const std::vector<Real>& Bext_ref, const std::vector<Real>& Bext) {
    auto err = Bext;
    Real max_err = 0, max_val = 0;
    for (long i = 0; i < (long)err.size(); i++) err[i] -= Bext_ref[i];
    for (const auto& x:B  ) max_val = std::max<Real>(max_val,fabs(x));
    for (const auto& x:err) max_err = std::max<Real>(max_err,fabs(x));
    std::cout<<"Maximum relative error ("<<name<<"): "<<max_err / max_val<<'\n';
    //if (err.size()%(trg_Nt*trg_Np)==0) biest::WriteVTK(("err_"+name).c_str(), NFP, half_period, Nt, Np, sctl::Vector<Real>(X), trg_Nt, trg_Np, sctl::Vector<Real>(err));
  };

  // Setup
  VirtualCasing<Real> virtual_casing;
  virtual_casing.Setup(digits, NFP, half_period, Nt, Np, X, src_Nt, src_Np, trg_Nt, trg_Np);

  // Get off-surface target points
  std::vector<Real> Xtrg_offsurf_int, Xtrg_offsurf_ext, Xn = virtual_casing.GetNormal(NFP, half_period, Nt, Np);
  for (long i = 0; i < (long)X.size(); i++) Xtrg_offsurf_int.push_back(X[i] - 0.08*Xn[i]);
  for (long i = 0; i < (long)X.size(); i++) Xtrg_offsurf_ext.push_back(X[i] + 0.08*Xn[i]);

  // Generate B fields for testing virtual-casing principle
  std::vector<Real> B, Bext, Bint;
  std::vector<Real> GradB, GradBext, GradBint;
  std::vector<Real> Bext_offsurf, Bint_offsurf;
  { // Set B, Bext, Bint
    std::vector<Real> Bint_, Bext_;
    std::tie(Bext, Bint) = VirtualCasingTestData<Real>::BFieldData(NFP, half_period, Nt, Np, X, trg_Nt, trg_Np);
    std::tie(Bext_, Bint_) = VirtualCasingTestData<Real>::BFieldData(NFP, half_period, Nt, Np, X, src_Nt, src_Np);
    const auto B_ = sctl::Vector<Real>(Bint_) + sctl::Vector<Real>(Bext_);
    B.assign(B_.begin(), B_.end());
  }
  { // Set GradB, GradBext
    std::tie(GradBext, GradBint) = VirtualCasingTestData<Real>::GradBFieldData(NFP, half_period, Nt, Np, X, trg_Nt, trg_Np);
    const auto GradB_ = sctl::Vector<Real>(GradBint) + sctl::Vector<Real>(GradBext);
    GradB.assign(GradB_.begin(), GradB_.end());
  }
  { // Set Bint_offsurf, Bext_offsurf
    std::tie(Bext_offsurf, std::ignore) = VirtualCasingTestData<Real>::BFieldDataOffSurf(NFP, half_period, Nt, Np, X, Xtrg_offsurf_int);
    std::tie(std::ignore, Bint_offsurf) = VirtualCasingTestData<Real>::BFieldDataOffSurf(NFP, half_period, Nt, Np, X, Xtrg_offsurf_ext);
  }

  // Compute Bext field
  auto Bext_ = virtual_casing.ComputeBext(B);
  auto Bint_ = virtual_casing.ComputeBint(B);
  print_relative_err("Bext", B, Bext, Bext_);
  print_relative_err("Bint", B, Bint, Bint_);

  // Compute GradBext field
  auto GradBext_ = virtual_casing.ComputeGradBext(B);
  print_relative_err("GradBext", GradB, GradBext, GradBext_);

  // Compute off-surface error
  const auto Bext_offsurf_ = virtual_casing.ComputeBextOffSurf(B, Xtrg_offsurf_int, (half_period?1:2)*800, 800); // max_Nt, max_Np = 800
  const auto Bint_offsurf_ = virtual_casing.ComputeBintOffSurf(B, Xtrg_offsurf_ext, (half_period?1:2)*800, 800); // max_Nt, max_Np = 800
  print_relative_err("Bext_offsurf", B, Bext_offsurf, Bext_offsurf_);
  print_relative_err("Bint_offsurf", B, Bint_offsurf, Bint_offsurf_);
  std::cout<<'\n';

  if (0) { // Visualize off-surface evaluation error
    const long N = 60;
    const Real L = 1.6, offset[]={5.8,0,0};

    std::vector<Real> Xtrg(3*N*N*N); // N x N cube of target points
    for (long i = 0; i < N; i++) {
      for (long j = 0; j < N; j++) {
        for (long k = 0; k < N; k++) {
          const long N3 = N*N*N;
          const long idx = (i*N+j)*N+k;
          Xtrg[0*N3+idx] = ((i/(Real)(N-1))-0.5)*L + offset[0];
          Xtrg[1*N3+idx] = ((j/(Real)(N-1))-0.5)*L + offset[1];
          Xtrg[2*N3+idx] = ((k/(Real)(N-1))-0.5)*L + offset[2];
        }
      }
    }

    std::vector<Real> Bext_offsurf; // reference data
    std::tie(Bext_offsurf, std::ignore) = VirtualCasingTestData<Real>::BFieldDataOffSurf(NFP, half_period, Nt, Np, X, Xtrg);
    const auto Bext_offsurf_ = virtual_casing.ComputeBextOffSurf(B, Xtrg, (half_period?1:2)*400, 400); // max_Nt, max_Np = 400

    { // Write VTK
      sctl::VTUData vtu_data;
      for (long i = 0; i < N; i++) {
        for (long j = 0; j < N; j++) {
          for (long k = 0; k < N; k++) {
            const long N3 = N*N*N;
            const long idx = (i*N+j)*N+k;
            vtu_data.coord.PushBack((float)Xtrg[0*N3+idx]);
            vtu_data.coord.PushBack((float)Xtrg[1*N3+idx]);
            vtu_data.coord.PushBack((float)Xtrg[2*N3+idx]);

            vtu_data.value.PushBack((float)(Bext_offsurf_[0*N3+idx]-Bext_offsurf[0*N3+idx]));
            vtu_data.value.PushBack((float)(Bext_offsurf_[1*N3+idx]-Bext_offsurf[1*N3+idx]));
            vtu_data.value.PushBack((float)(Bext_offsurf_[2*N3+idx]-Bext_offsurf[2*N3+idx]));
          }
        }
      }
      for (long i = 0; i < N-1; i++) {
        for (long j = 0; j < N-1; j++) {
          for (long k = 0; k < N-1; k++) {
            auto idx = [N](long i, long j, long k) {
              return (i*N+j)*N+k;
            };
            vtu_data.connect.PushBack(idx(i+0,j+0,k+0));
            vtu_data.connect.PushBack(idx(i+0,j+0,k+1));
            vtu_data.connect.PushBack(idx(i+0,j+1,k+1));
            vtu_data.connect.PushBack(idx(i+0,j+1,k+0));
            vtu_data.connect.PushBack(idx(i+1,j+0,k+0));
            vtu_data.connect.PushBack(idx(i+1,j+0,k+1));
            vtu_data.connect.PushBack(idx(i+1,j+1,k+1));
            vtu_data.connect.PushBack(idx(i+1,j+1,k+0));
            vtu_data.offset.PushBack(vtu_data.connect.Dim());
            vtu_data.types.PushBack(12);
          }
        }
      }
      vtu_data.WriteVTK("Bext_offsurf_err", sctl::Comm::Self());
    }
  }
}

int main() {
  sctl::Profile::Enable(true);
  for (long digits = 3; digits <= 12; digits+=3) {
    //test<double>(digits, 5, false, 1, 4, biest::SurfType::AxisymNarrow, 2*digits, 7*digits, 20, 20);

    test<double>(digits, 5, false, 20, 20, biest::SurfType::W7X_, 12*digits, 32*digits, 40, 40);
    test<double>(digits, 5, true, 10, 20, biest::SurfType::W7X_, 6*digits, 32*digits, 20, 40);
  }
  return 0;
}
