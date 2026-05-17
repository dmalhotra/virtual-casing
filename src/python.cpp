#include <vector>
#include <virtual-casing.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;


PYBIND11_MODULE(virtual_casing, m) {
    m.doc() = "Virtual casing principle for magnetic field computation.";

    py::class_<sctl::Vector<double>>(m, "SCTLDoubleVector")
         .def(py::init<>())
         .def(py::init<const std::vector<double> & >(), py::arg("data"));

    py::class_<sctl::Vector<int>>(m, "SCTLIntVector")
        .def(py::init<>())
        .def(py::init<const std::vector<int> & >(), py::arg("data"));

    py::enum_<biest::SurfType>(m, "SurfType",
        "Prebuilt surface geometry types for testing.")
        .value("AxisymCircleWide", biest::SurfType::AxisymCircleWide)           // 125 x 50
        .value("AxisymCircleNarrow", biest::SurfType::AxisymCircleNarrow)       // 250 x 50
        .value("AxisymWide", biest::SurfType::AxisymWide)                       // 125 x 50
        .value("AxisymNarrow", biest::SurfType::AxisymNarrow)                   // 250 x 50
        .value("RotatingEllipseWide", biest::SurfType::RotatingEllipseWide)     // 125 x 50
        .value("RotatingEllipseNarrow", biest::SurfType::RotatingEllipseNarrow) // 250 x 50
        .value("Quas3", biest::SurfType::Quas3)                                 // 250 x 50
        .value("LHD", biest::SurfType::LHD)                                     // 250 x 50
        .value("W7X", biest::SurfType::W7X)                                     // 250 x 45
        .value("Stell", biest::SurfType::Stell)                                 // 250 x 50
        .value("W7X_", biest::SurfType::W7X_);

    py::class_<VirtualCasing<double>>(m, "VirtualCasing",
        "Compute the external or internal magnetic field using the virtual casing principle.")
        .def(py::init<>())
        .def("setup", &VirtualCasing<double>::Setup,
            R"pbdoc(
            Setup the VirtualCasing object.

            Parameters
            ----------
            digits : int
                Number of decimal digits of accuracy.
            NFP : int
                Number of toroidal field periods. The surface as well as the
                magnetic field must have this toroidal periodic symmetry.
            half_period : bool
                Whether the surface and data are defined on half field period.
            Nt : int
                Surface discretization order in toroidal direction (in one field period).
            Np : int
                Surface discretization order in poloidal direction.
            X : list of float
                The surface coordinates in the order {x11, x12, ..., x1Np,
                x21, x22, ..., xNtNp, y11, ..., z11, ...}.
            src_Nt : int
                Input B-field discretization order in toroidal direction (in one field period).
            src_Np : int
                Input B-field discretization order in poloidal direction.
            trg_Nt : int, optional
                Output Bext-field discretization order in toroidal direction
                (in one field period). Default is -1 (same as src_Nt).
            trg_Np : int, optional
                Output Bext-field discretization order in poloidal direction.
                Default is -1 (same as src_Np).

            Notes
            -----
            The grids for input and output data differ depending on whether
            stellarator symmetry is exploited (half_period). For this discussion,
            consider the toroidal angle phi and poloidal angle theta to have
            period 1 (not 2*pi).

            If half_period is False, the grids in toroidal angles begin at phi=0.
            The grid spacing is 1/(NFP*Nt), with no point at phi = 1/NFP.

            If half_period is True, the toroidal grids are shifted by half a grid
            point, so there is no grid point at phi=0. The phi grid for the surface
            has points at 0.5/(NFP*Nt), 1.5/(NFP*Nt), ..., (Nt-0.5)/(NFP*Nt).

            The poloidal grid always ranges uniformly over [0, 1), with the first
            grid point at theta=0 and no grid point at theta=1.

            The resolution parameters for the surface shape (Nt, Np), input magnetic
            field (src_Nt, src_Np), and output external field (trg_Nt, trg_Np) do
            not need to be related to each other in any particular way.
            )pbdoc",
            py::arg("digits"), py::arg("NFP"), py::arg("half_period"),
            py::arg("Nt"), py::arg("Np"), py::arg("X"),
            py::arg("src_Nt"), py::arg("src_Np"),
            py::arg("trg_Nt") = -1, py::arg("trg_Np") = -1)
        .def("compute_external_B", &VirtualCasing<double>::ComputeBext,
            py::arg("B_total"),
            R"pbdoc(
            Compute the magnetic field due to currents external to the surface.

            Recover the Bext component from the total field B = Bext + Bint:
            Bext = B/2 + gradG[B . n] + BiotSavart[n x B]

            Here, Bext is the magnetic field due to currents in the exterior of the
            surface, and Bint is the magnetic field due to currents in the interior.

            Parameters
            ----------
            B : list of float
                The total magnetic field on the surface due to all currents.
                B = {Bx11, Bx12, ..., Bx1Np, Bx21, Bx22, ..., BxNtNp, By11, ..., Bz11, ...},
                where Nt and Np are the number of discretizations in toroidal and
                poloidal directions.

            Returns
            -------
            list of float
                The component of magnetic field on the surface due to currents in
                the exterior of the surface.
            )pbdoc")
        .def("compute_external_B_offsurf", &VirtualCasing<double>::ComputeBextOffSurf,
            py::arg("B_total"), py::arg("Xt"), py::arg("max_Nt") = -1, py::arg("max_Np") = -1,
            R"pbdoc(
            Compute the magnetic field due to currents external to the surface at off-surface points.

            Recover the Bext component from the total field B = Bext + Bint:
            Bext = gradG[B . n] + BiotSavart[n x B]

            Parameters
            ----------
            B : list of float
                The total magnetic field on the surface due to all currents.
                B = {Bx11, Bx12, ..., Bx1Np, Bx21, Bx22, ..., BxNtNp, By11, ..., Bz11, ...}.
            Xt : list of float
                The coordinates for off-surface evaluation points in the order
                {x1, x2, ..., xn, y1, ..., z1, ..., zn}.
            max_Nt : int, optional
                Restrict upsampling to max_Nt modes (in a field-period) in toroidal
                direction. Default is -1 (no restriction).
            max_Np : int, optional
                Restrict upsampling to max_Np modes in poloidal direction.
                Default is -1 (no restriction).

            Returns
            -------
            list of float
                The component of magnetic field at the evaluation points due to
                currents in the exterior of the surface.
            )pbdoc")
        .def("compute_external_gradB", &VirtualCasing<double>::ComputeGradBext,
            py::arg("B_total"),
            R"pbdoc(
            Compute the gradient of the magnetic field due to currents external to the surface.

            Recover GradBext from the total field B = Bext + Bint using the
            virtual casing principle.

            Parameters
            ----------
            B : list of float
                The total magnetic field on the surface due to all currents.
                B = {Bx11, Bx12, ..., Bx1Np, Bx21, Bx22, ..., BxNtNp, By11, ..., Bz11, ...}.

            Returns
            -------
            list of float
                The gradient of the magnetic field component on the surface due to
                currents in the exterior of the surface.
            )pbdoc")
        .def("compute_internal_B", &VirtualCasing<double>::ComputeBint,
            py::arg("B_total"),
            R"pbdoc(
            Compute the magnetic field due to currents internal to the surface.

            Recover the Bint component from the total field B = Bext + Bint:
            Bint = B/2 - gradG[B . n] - BiotSavart[n x B]

            Here, Bext is the magnetic field due to currents in the exterior of the
            surface, and Bint is the magnetic field due to currents in the interior.

            Parameters
            ----------
            B : list of float
                The total magnetic field on the surface due to all currents.
                B = {Bx11, Bx12, ..., Bx1Np, Bx21, Bx22, ..., BxNtNp, By11, ..., Bz11, ...}.

            Returns
            -------
            list of float
                The component of magnetic field on the surface due to currents in
                the interior of the surface.
            )pbdoc")
        .def("compute_internal_B_offsurf", &VirtualCasing<double>::ComputeBintOffSurf,
            py::arg("B_total"), py::arg("Xt"), py::arg("max_Nt") = -1, py::arg("max_Np") = -1,
            R"pbdoc(
            Compute the magnetic field due to currents internal to the surface at off-surface points.

            Recover the Bint component from the total field B = Bext + Bint.

            Parameters
            ----------
            B : list of float
                The total magnetic field on the surface due to all currents.
                B = {Bx11, Bx12, ..., Bx1Np, Bx21, Bx22, ..., BxNtNp, By11, ..., Bz11, ...}.
            Xt : list of float
                The coordinates for off-surface evaluation points in the order
                {x1, x2, ..., xn, y1, ..., z1, ..., zn}.
            max_Nt : int, optional
                Restrict upsampling to max_Nt modes (in a field-period) in toroidal
                direction. Default is -1 (no restriction).
            max_Np : int, optional
                Restrict upsampling to max_Np modes in poloidal direction.
                Default is -1 (no restriction).

            Returns
            -------
            list of float
                The component of magnetic field at the evaluation points due to
                currents in the interior of the surface.
            )pbdoc");

    py::class_<VirtualCasingTestData<double>>(m, "VirtualCasingTestData",
        "Generate test data for class VirtualCasing.")
        .def(py::init<>())
        .def_static("surface_coordinates", &VirtualCasingTestData<double>::SurfaceCoordinates,
            py::arg("NFP"), py::arg("half_period"), py::arg("Nt"), py::arg("Np"),
            py::arg("surf_type") = biest::SurfType::AxisymNarrow,
            R"pbdoc(
            Generate nodal coordinates for toroidal surfaces.

            Parameters
            ----------
            NFP : int
                Number of toroidal field periods.
            half_period : bool
                Whether the returned surface coordinates should be on half field period.
            Nt : int
                Surface discretization order in toroidal direction (in one field period).
            Np : int
                Surface discretization order in poloidal direction.
            surf_type : SurfType, optional
                Prebuilt surface geometry type. Default is SurfType.AxisymNarrow.
                Possible values: AxisymCircleWide, AxisymCircleNarrow, AxisymWide,
                AxisymNarrow, RotatingEllipseWide, RotatingEllipseNarrow, Quas3,
                LHD, W7X, Stell.

            Returns
            -------
            list of float
                The surface coordinates in the order {x11, x12, ..., x1Np,
                x21, x22, ..., xNtNp, y11, ..., z11, ...}. The coordinates
                correspond to the surface in the toroidal angle interval [0, 2*pi/NFP).
            )pbdoc")
        .def_static("magnetic_field_data", &VirtualCasingTestData<double>::BFieldData,
            py::arg("NFP"), py::arg("half_period"), py::arg("Nt"), py::arg("Np"),
            py::arg("X"), py::arg("trg_Nt"), py::arg("trg_Np"),
            R"pbdoc(
            Generate B field data for testing with class VirtualCasing.

            Parameters
            ----------
            NFP : int
                Number of toroidal field periods.
            half_period : bool
                Whether the result should be on half field period.
            Nt : int
                Surface discretization order in toroidal direction (in one field period).
            Np : int
                Surface discretization order in poloidal direction.
            X : list of float
                The surface coordinates in the order {x11, x12, ..., x1Np,
                x21, x22, ..., xNtNp, y11, ..., z11, ...}.
            trg_Nt : int
                Output B-field discretization order in toroidal direction (in one field period).
            trg_Np : int
                Output B-field discretization order in poloidal direction.

            Returns
            -------
            tuple of (list of float, list of float)
                Bext and Bint, magnetic fields generated by an internal current loop
                and an external current loop respectively.
            )pbdoc")
        .def_static("magnetic_field_grad_data", &VirtualCasingTestData<double>::GradBFieldData,
            py::arg("NFP"), py::arg("half_period"), py::arg("Nt"), py::arg("Np"),
            py::arg("X"), py::arg("trg_Nt"), py::arg("trg_Np"),
            R"pbdoc(
            Generate gradient of B field data for testing with class VirtualCasing.

            Parameters
            ----------
            NFP : int
                Number of toroidal field periods.
            half_period : bool
                Whether the result should be on half field period.
            Nt : int
                Surface discretization order in toroidal direction (in one field period).
            Np : int
                Surface discretization order in poloidal direction.
            X : list of float
                The surface coordinates in the order {x11, x12, ..., x1Np,
                x21, x22, ..., xNtNp, y11, ..., z11, ...}.
            trg_Nt : int
                Output B-field discretization order in toroidal direction (in one field period).
            trg_Np : int
                Output B-field discretization order in poloidal direction.

            Returns
            -------
            tuple of (list of float, list of float)
                GradBext and GradBint, gradients of magnetic fields generated by an
                internal current loop and an external current loop respectively.
            )pbdoc");

}

