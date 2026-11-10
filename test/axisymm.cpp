// Copyright (C) 2021-2026 The Simons Foundation, Inc. and the virtual-casing authors (see NOTICE)
// Authors: Dhairya Malhotra, Antoine Cerfon
// SPDX-License-Identifier: Apache-2.0

#include <virtual-casing.hpp>

template <class Real> void test(const Real R0, const Real a, const Real kappa, int digits, int NFP, long Nt, long Np, long src_Nt, long src_Np, long trg_Nt, long trg_Np) {
  // Construct the surface
  std::vector<Real> X(3*Nt*Np);
  //X = VirtualCasingTestData<Real>::SurfaceCoordinates(NFP, Nt, Np, biest::SurfType::W7X_);
  for (long t = 0; t < Nt; t++) { // toroidal direction
    for (long p = 0; p < Np; p++) { // poloidal direction
      Real theta = 2*sctl::const_pi<Real>()*p/Np;
      Real phi   = 2*sctl::const_pi<Real>()*t/(NFP*Nt);

      Real R = sctl::sqrt<Real>(R0*R0 + 2*a*R0*sctl::cos<Real>(theta));
      Real x = R * sctl::cos<Real>(phi);
      Real y = R * sctl::sin<Real>(phi);
      Real z = kappa*a*R0/R*sctl::sin<Real>(theta);

      X[(0*Nt+t)*Np+p] = x;
      X[(1*Nt+t)*Np+p] = y;
      X[(2*Nt+t)*Np+p] = z;
    }
  }

  // Generate B field
  std::vector<Real> B(3*src_Nt*src_Np);
  for (long t = 0; t < src_Nt; t++) { // toroidal direction
    for (long p = 0; p < src_Np; p++) { // poloidal direction
      Real theta = 2*sctl::const_pi<Real>()*p/src_Np;
      Real phi   = 2*sctl::const_pi<Real>()*t/(NFP*src_Nt);

      B[(0*src_Nt+t)*src_Np+p] = -a*sctl::sin<Real>(theta) * sctl::cos<Real>(phi);
      B[(1*src_Nt+t)*src_Np+p] = -a*sctl::sin<Real>(theta) * sctl::sin<Real>(phi);
      B[(2*src_Nt+t)*src_Np+p] = a*kappa*sctl::cos<Real>(theta) + kappa*a*a*sctl::sin<Real>(theta)*sctl::sin<Real>(theta) / (1+2*a*sctl::cos<Real>(theta));
    }
  }

  // Setup
  VirtualCasing<Real> virtual_casing;
  virtual_casing.Setup(digits, NFP, false, Nt, Np, X, src_Nt, src_Np, trg_Nt, trg_Np);

  // Compute Bext field
  std::vector<Real> Bext = virtual_casing.ComputeBext(B);

  // Compute Bext . normal
  std::vector<Real> normal, Bext_dot_n(trg_Nt*trg_Np);
  normal = virtual_casing.GetNormal(NFP, false, trg_Nt, trg_Np);
  for (long i = 0; i < trg_Nt*trg_Np; i++) {
    Real dot_prod = 0;
    for (long k = 0; k < 3; k++) {
      dot_prod += Bext[k*trg_Nt*trg_Np+i] * normal[k*trg_Nt*trg_Np+i];
    }
    Bext_dot_n[i] = dot_prod;
  }

  std::cout<<"\n\nBext = \n";
  for (long k = 0; k < 3; k++) { // Print Bext_
    for (long t = 0; t < trg_Nt; t++) {
      for (long p = 0; p < trg_Np; p++) {
        printf("%10.5f", (double)Bext[(k*trg_Nt+t)*trg_Np+p]);
      }
      std::cout<<'\n';
    }
  }

  std::cout<<"\n\nBext . normal = \n";
  for (long t = 0; t < trg_Nt; t++) { // Print Bext_dot_n
    for (long p = 0; p < trg_Np; p++) {
      printf("%10.5f", (double)Bext_dot_n[t*trg_Np+p]);
    }
    std::cout<<'\n';
  }

  // Generate VTK visualization
  biest::WriteVTK("B", NFP, false, Nt, Np, sctl::Vector<Real>(X), src_Nt, src_Np, sctl::Vector<Real>(B), biest::Parity::Odd);
  biest::WriteVTK("Bext", NFP, false, Nt, Np, sctl::Vector<Real>(X), trg_Nt, trg_Np, sctl::Vector<Real>(Bext), biest::Parity::Odd);
  biest::WriteVTK("Bext_dot_n", NFP, false, Nt, Np, sctl::Vector<Real>(X), trg_Nt, trg_Np, sctl::Vector<Real>(Bext_dot_n), biest::Parity::Odd);
}

int main() {
  long digits = 10;
  double R0 = 1, a = 1./3, kappa = 1;

  const long NFP = 20, Nt = 1, Np = 15;
  test<double>(R0, a, kappa, digits, NFP, Nt, Np, Nt, Np, Nt, Np);

  return 0;
}
