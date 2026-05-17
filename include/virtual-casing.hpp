#include <utility>
#include <biest.hpp>

#ifndef _BIEST_VIRTUAL_CASING_HPP_
#define _BIEST_VIRTUAL_CASING_HPP_

template <class Real> class VirtualCasing {
  static constexpr sctl::Integer COORD_DIM = 3;

  public:

    VirtualCasing();

    /**
     * Setup the VirtualCasing object.
     *
     * @param[in] digits number of decimal digits of accuracy.
     *
     * @param[in] NFP number of toroidal field periods. The surface as well as
     * the magnetic field must have this toroidal periodic symmetry.
     *
     * @param[in] half_period whether the surface and data are defined on half field period.
     *
     * @param[in] Nt surface discretization order in toroidal direction
     * (in one field period or half field period when half_period=true).
     *
     * @param[in] Np surface discretization order in poloidal direction.
     *
     * @param[in] X the surface coordinates in the order {x11, x12, ..., x1Np,
     * x21, x22, ... , xNtNp, y11, ... , z11, ...}.
     *
     * @param[in] src_Nt input B-field discretization order in toroidal direction
     * (in one field period or half field period when half_period=true).
     *
     * @param[in] src_Np input B-field discretization order in poloidal direction.
     *
     * @param[in] trg_Nt output Bext-field discretization order in toroidal direction
     * (in one field period or half field period when half_period=true).
     *
     * @param[in] trg_Np output Bext-field discretization order in poloidal direction.
     *
     * The grids for the input and output data differ in an important
     * way depending on whether or not stellarator symmetry is
     * exploited, i.e. on half_period. For this discussion, consider
     * the toroidal angle phi and poloidal angle theta to have period
     * 1 (not 2*pi).
     *
     * If you do not exploit stellarator symmetry, then half_period
     * is set to false.  In this case the grids in the toroidal angles
     * begin at phi=0.  The grid spacing is 1 / (NFP * Nt), and there
     * is no point at the symmetry plane phi = 1 / NFP.
     *
     * If you do wish to exploit stellarator symmetry, set half_period
     * to true. The grid spacing is 1 / (2 * NFP * Nt).
     * In this case the toroidal grids are each shifted by
     * half a grid point, so there is no grid point at phi = 0. The
     * phi grid for the surface shape has points at 0.5 / (2 * NFP * Nt),
     * 1.5 / (2 * NFP * Nt), ..., (Nt - 0.5) / (2 * NFP * Nt). The phi grid
     * for the output B_external has grid points at 0.5 / (2 * NFP *
     * trg_Nt), 1.5 / (2 * NFP * trg_Nt), ..., (trg_Nt - 0.5) / (2 * NFP *
     * Nt), and similarly for the input B field with trg -> src.
     *
     * The rationale for these conventions is that in both the
     * stellarator-symmetric and non-stellarator-symmetric cases,
     * integration over the surface can be achieved with spectral
     * accuracy on these grids using uniform weights.
     *
     * Regardless of half_period, the poloidal grid always ranges
     * uniformly over [0, 1), with the first grid point at theta = 0,
     * and no grid point at theta = 1.
     *
     * The resolution parameters for the surface shape (Nt, Np), input
     * magnetic field (src_Nt, src_Np), and output external field
     * (trg_Nt, trg_Np) do not need to be related to each other in any
     * particular way. For example, there is no performance penalty if
     * src_Nt and trg_Nt are relatively prime.
     */
    void Setup(const sctl::Integer digits, const sctl::Integer NFP, const bool half_period, const sctl::Long Nt, const sctl::Long Np, const std::vector<Real>& X, const sctl::Long src_Nt, const sctl::Long src_Np, const sctl::Long trg_Nt = -1, const sctl::Long trg_Np = -1);

    /**
     * Recover the Bext component from the total field B = Bext + Bint by
     * applying the virtual-casing principle:
     * Bext = B/2 + gradG[B . n] + BiotSavart[n x B]
     *
     * Here, Bext is the magnetic field due to currents in the exterior of the surface,
     * and Bint is the magnetic field due to currents in the interior of the surface.
     *
     * @param[in] B the total magnetic field on the surface due to all currents.
     * B = {Bx11, Bx12, ..., Bx1Np, Bx21, Bx22, ... , BxNtNp, By11, ... , Bz11, ...},
     * where Nt and Np are the number of discretizations in toroidal and
     * poloidal directions.
     *
     * @return the component of magnetic field on the surface due to
     * currents in the exterior of the surface, computed using the
     * virtual-casing principle.
     */
    std::vector<Real> ComputeBext(const std::vector<Real>& B) const;

    /**
     * Recover the Bext component from the total field B = Bext + Bint by
     * applying the virtual-casing principle (for off-surface target points):
     * Bext = gradG[B . n] + BiotSavart[n x B]
     *
     * Here, Bext is the magnetic field due to currents in the exterior of the surface,
     * and Bint is the magnetic field due to currents in the interior of the surface.
     *
     * @param[in] B the total magnetic field on the surface due to all currents.
     * B = {Bx11, Bx12, ..., Bx1Np, Bx21, Bx22, ... , BxNtNp, By11, ... , Bz11, ...},
     * where Nt and Np are the number of discretizations in toroidal and
     * poloidal directions.
     *
     * @param[in] Xt the vector of coordinates for the off-surface evaluation
     * points in the order {x1, x2, ..., xn, y1, ..., z1, ..., Zn}.
     *
     * @param[in] max_Nt restrict upsampling to max_Nt modes (in a field-period) in toroidal direction.
     *
     * @param[in] max_Np restrict upsampling to max_Np modes in poloidal direction.
     *
     * @return the component of magnetic field at the evaluation points due to
     * currents in the exterior of the surface, computed using the
     * virtual-casing principle.
     */
    std::vector<Real> ComputeBextOffSurf(const std::vector<Real>& B, const std::vector<Real>& Xt, const sctl::Long max_Nt=-1, const sctl::Long max_Np=-1) const;

    /**
     * Recover the GradBext component from the total field B = Bext + Bint by
     * applying the virtual-casing principle:
     * Bext = B/2 + gradG[B . n] + BiotSavart[n x B]
     *
     * @param[in] B the total magnetic field on the surface due to all currents.
     * B = {Bx11, Bx12, ..., Bx1Np, Bx21, Bx22, ... , BxNtNp, By11, ... , Bz11, ...},
     * where Nt and Np are the number of discretizations in toroidal and
     * poloidal directions.
     *
     * @return the gradient of the component of magnetic field on the surface
     * due to currents in the exterior of the surface, computed using the
     * virtual-casing principle.
     */
    std::vector<Real> ComputeGradBext(const std::vector<Real>& B) const;

    /**
     * Recover the Bint component from the total field B = Bext + Bint by
     * applying the virtual-casing principle:
     * Bint = B/2 - gradG[B . n] - BiotSavart[n x B]
     *
     * Here, Bext is the magnetic field due to currents in the exterior of the surface,
     * and Bint is the magnetic field due to currents in the interior of the surface.
     *
     * @param[in] B the total magnetic field on the surface due to all currents.
     * B = {Bx11, Bx12, ..., Bx1Np, Bx21, Bx22, ... , BxNtNp, By11, ... , Bz11, ...},
     * where Nt and Np are the number of discretizations in toroidal and
     * poloidal directions.
     *
     * @return the component of magnetic field on the surface due to
     * currents in the interior of the surface.
     */
    std::vector<Real> ComputeBint(const std::vector<Real>& B) const;

    /**
     * Recover the Bint component from the total field B = Bext + Bint (for off-surface target points).
     *
     * Here, Bext is the magnetic field due to currents in the exterior of the surface,
     * and Bint is the magnetic field due to currents in the interior of the surface.
     *
     * @param[in] B the total magnetic field on the surface due to all currents.
     * B = {Bx11, Bx12, ..., Bx1Np, Bx21, Bx22, ... , BxNtNp, By11, ... , Bz11, ...},
     * where Nt and Np are the number of discretizations in toroidal and
     * poloidal directions.
     *
     * @param[in] Xt the vector of coordinates for the off-surface evaluation
     * points in the order {x1, x2, ..., xn, y1, ..., z1, ..., Zn}.
     *
     * @param[in] max_Nt restrict upsampling to max_Nt modes (in a field-period) in toroidal direction.
     *
     * @param[in] max_Np restrict upsampling to max_Np modes in poloidal direction.
     *
     * @return the component of magnetic field at the evaluation points due to
     * currents in the interior of the surface.
     */
    std::vector<Real> ComputeBintOffSurf(const std::vector<Real>& B, const std::vector<Real>& Xt, const sctl::Long max_Nt=-1, const sctl::Long max_Np=-1) const;

    //std::vector<Real> ComputeGradBint(const std::vector<Real>& B) const; // TODO: requires Hedgehog quadrature normal to be flipped

    /**
     * Returns the surface normal vectors.
     *
     * @param[in] NFP number of toroidal field periods. The result will be on
     * one field period (or half period when half_period=true).
     *
     * @param[in] half_period whether the result should be on half field period.
     *
     * @param[in] Nt surface discretization order in toroidal direction
     * (in one field period or half field period when half_period=true).
     *
     * @param[in] Np surface discretization order in poloidal direction.
     *
     * @return the surface normal vectors in the order {nx11, nx12, ..., nx1Np,
     * nx21, nx22, ... , nxNtNp, ny11, ... , nz11, ...}.
     */
    std::vector<Real> GetNormal(const sctl::Integer NFP, const bool half_period, const sctl::Long Nt, const sctl::Long Np) const;

  private:

    std::vector<Real> ComputeB(const std::vector<Real>& B, bool ext) const;

    std::vector<Real> ComputeBOffSurf(const std::vector<Real>& B, const std::vector<Real>& Xt, const sctl::Long max_Nt, const sctl::Long max_Np, bool ext) const;

    std::vector<Real> ComputeGradB(const std::vector<Real>& B, bool ext) const;


    static void DotProd(sctl::Vector<Real>& AdotB, const sctl::Vector<Real>& A, const sctl::Vector<Real>& B);

    static void CrossProd(sctl::Vector<Real>& AcrossB, const sctl::Vector<Real>& A, const sctl::Vector<Real>& B);

    sctl::Comm comm_;
    mutable biest::FieldPeriodBIOp<Real,COORD_DIM,3,3> BiotSavartFxU;
    mutable biest::FieldPeriodBIOp<Real,COORD_DIM,1,3> LaplaceFxdU;
    mutable biest::FieldPeriodBIOp<Real,COORD_DIM,3,9, 8> BiotSavartFxdU;
    mutable biest::FieldPeriodBIOp<Real,COORD_DIM,1,9, 8> LaplaceFxd2U;
    sctl::Vector<biest::Surface<Real>> Svec;
    bool half_period_;
    sctl::Integer NFP_, digits_;
    sctl::Long src_Nt_, src_Np_;
    sctl::Long trg_Nt_, trg_Np_;
    mutable sctl::Long quad_Nt_, quad_Np_;
    mutable sctl::Long grad_quad_Nt_, grad_quad_Np_;
    mutable sctl::Vector<Real> dX, normal;
    mutable bool dosetup, dosetup_grad;
};

/**
 * Generate data for testing class VirtualCasing.
 */
template <class Real> class VirtualCasingTestData {
  static constexpr int COORD_DIM = 3;

  public:

    /**
     * Generate nodal coordinates for toroidal surfaces.
     *
     * @param[in] NFP number of toroidal field periods.
     *
     * @param[in] half_period whether the returned surface coordinates should be on half field period.
     *
     * @param[in] Nt surface discretization order in toroidal direction
     * (in one field period or half field period when half_period=true).
     *
     * @param[in] Np surface discretization order in poloidal direction.
     *
     * @param[in] surf_type prebuilt surface geometries. Possible values
     * biest::SurfType::{AxisymCircleWide, AxisymCircleNarrow, AxisymWide,
     * AxisymNarrow, RotatingEllipseWide, RotatingEllipseNarrow, Quas3, LHD,
     * W7X, Stell}
     *
     * @return the surface coordinates in the order {x11, x12, ..., x1Np,
     * x21, x22, ... , xNtNp, y11, ... , z11, ...}. The coordinates correspond
     * to the surface in the toroidal angle interval [0, 2*pi/NFP).
     */
    static std::vector<Real> SurfaceCoordinates(const sctl::Integer NFP, const bool half_period, const sctl::Long Nt, const sctl::Long Np, const biest::SurfType surf_type = biest::SurfType::AxisymNarrow);

    /**
     * Generate B field data to be used with class VirtualCasing.
     *
     * @param[in] NFP number of toroidal field periods.
     *
     * @param[in] half_period whether the result should be on half field period.
     *
     * @param[in] Nt surface discretization order in toroidal direction
     * (in one field period or half field period when half_period=true).
     *
     * @param[in] Np surface discretization order in poloidal direction.
     *
     * @param[in] X the surface coordinates in the order {x11, x12, ..., x1Np,
     * x21, x22, ... , xNtNp, y11, ... , z11, ...}.
     *
     * @param[in] trg_Nt output B-field discretization order in toroidal direction
     * (in one field period or half field period when half_period=true).
     *
     * @param[in] trg_Np output B-field discretization order in poloidal direction.
     *
     * @return Bext and Bint, magnetic fields generated by an internal current loop and
     *  an external current loop respectively.
     */
    static std::tuple<std::vector<Real>, std::vector<Real>> BFieldData(const sctl::Integer NFP, const bool half_period, const sctl::Long Nt, const sctl::Long Np, const std::vector<Real>& X, const sctl::Long trg_Nt, const sctl::Long trg_Np);

    /**
     * Generate B field data to be used with class VirtualCasing.
     *
     * @param[in] NFP number of toroidal field periods.
     *
     * @param[in] half_period whether the result should be on half field period.
     *
     * @param[in] Nt surface discretization order in toroidal direction
     * (in one field period or half field period when half_period=true).
     *
     * @param[in] Np surface discretization order in poloidal direction.
     *
     * @param[in] X the surface coordinates in the order {x11, x12, ..., x1Np,
     * x21, x22, ... , xNtNp, y11, ... , z11, ...}.
     *
     * @param[in] Xtrg the target coordinates in the order {x11, x12, ..., x1Np,
     * x21, x22, ... , xNtNp, y11, ... , z11, ...}.
     *
     * @return Bext and Bint, magnetic fields generated by an internal current loop and
     *  an external current loop respectively.
     */
    static std::tuple<std::vector<Real>, std::vector<Real>> BFieldDataOffSurf(const sctl::Integer NFP, const bool half_period, const sctl::Long Nt, const sctl::Long Np, const std::vector<Real>& X, const std::vector<Real>& Xtrg);

    /**
     * Generate gradient of B field data to be used with class VirtualCasing.
     *
     * @param[in] NFP number of toroidal field periods.
     *
     * @param[in] half_period whether the result should be on half field period.
     *
     * @param[in] Nt surface discretization order in toroidal direction
     * (in one field period or half field period when half_period=true).
     *
     * @param[in] Np surface discretization order in poloidal direction.
     *
     * @param[in] X the surface coordinates in the order {x11, x12, ..., x1Np,
     * x21, x22, ... , xNtNp, y11, ... , z11, ...}.
     *
     * @param[in] trg_Nt output B-field discretization order in toroidal direction
     * (in one field period or half field period when half_period=true).
     *
     * @param[in] trg_Np output B-field discretization order in poloidal direction.
     *
     * @return GradBext and GradBint, magnetic fields generated by an internal current loop and
     *  an external current loop respectively.
     */
    static std::tuple<std::vector<Real>, std::vector<Real>> GradBFieldData(const sctl::Integer NFP, const bool half_period, const sctl::Long Nt, const sctl::Long Np, const std::vector<Real>& X, const sctl::Long trg_Nt, const sctl::Long trg_Np);
};

#include <virtual-casing.txx>

#endif
