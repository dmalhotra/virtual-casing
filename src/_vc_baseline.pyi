"""
Virtual casing principle for magnetic field computation.
"""
from __future__ import annotations
import collections.abc
import typing
__all__: list[str] = ['SCTLDoubleVector', 'SCTLIntVector', 'SurfType', 'VirtualCasing', 'VirtualCasingTestData']
class SCTLDoubleVector:
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, data: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
class SCTLIntVector:
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, data: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
class SurfType:
    """
    Prebuilt surface geometry types for testing.

    Members:

      AxisymCircleWide

      AxisymCircleNarrow

      AxisymWide

      AxisymNarrow

      RotatingEllipseWide

      RotatingEllipseNarrow

      Quas3

      LHD

      W7X

      Stell

      W7X_
    """
    AxisymCircleNarrow: typing.ClassVar[SurfType]  # value = <SurfType.AxisymCircleNarrow: 1>
    AxisymCircleWide: typing.ClassVar[SurfType]  # value = <SurfType.AxisymCircleWide: 0>
    AxisymNarrow: typing.ClassVar[SurfType]  # value = <SurfType.AxisymNarrow: 3>
    AxisymWide: typing.ClassVar[SurfType]  # value = <SurfType.AxisymWide: 2>
    LHD: typing.ClassVar[SurfType]  # value = <SurfType.LHD: 7>
    Quas3: typing.ClassVar[SurfType]  # value = <SurfType.Quas3: 6>
    RotatingEllipseNarrow: typing.ClassVar[SurfType]  # value = <SurfType.RotatingEllipseNarrow: 5>
    RotatingEllipseWide: typing.ClassVar[SurfType]  # value = <SurfType.RotatingEllipseWide: 4>
    Stell: typing.ClassVar[SurfType]  # value = <SurfType.Stell: 9>
    W7X: typing.ClassVar[SurfType]  # value = <SurfType.W7X: 8>
    W7X_: typing.ClassVar[SurfType]  # value = <SurfType.W7X_: 10>
    __members__: typing.ClassVar[dict[str, SurfType]]  # value = {'AxisymCircleWide': <SurfType.AxisymCircleWide: 0>, 'AxisymCircleNarrow': <SurfType.AxisymCircleNarrow: 1>, 'AxisymWide': <SurfType.AxisymWide: 2>, 'AxisymNarrow': <SurfType.AxisymNarrow: 3>, 'RotatingEllipseWide': <SurfType.RotatingEllipseWide: 4>, 'RotatingEllipseNarrow': <SurfType.RotatingEllipseNarrow: 5>, 'Quas3': <SurfType.Quas3: 6>, 'LHD': <SurfType.LHD: 7>, 'W7X': <SurfType.W7X: 8>, 'Stell': <SurfType.Stell: 9>, 'W7X_': <SurfType.W7X_: 10>}
    @typing.overload
    def __eq__(self, other: SurfType) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    @typing.overload
    def __ne__(self, other: SurfType) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class VirtualCasing:
    """
    Compute the external or internal magnetic field using the virtual casing principle.
    """
    def __init__(self) -> None:
        ...
    def compute_external_B(self, B_total: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> list[float]:
        """
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
        """
    def compute_external_B_offsurf(self, B_total: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], Xt: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], max_Nt: typing.SupportsInt | typing.SupportsIndex = -1, max_Np: typing.SupportsInt | typing.SupportsIndex = -1) -> list[float]:
        """
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
        """
    def compute_external_gradB(self, B_total: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> list[float]:
        """
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
        """
    def compute_internal_B(self, B_total: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> list[float]:
        """
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
        """
    def compute_internal_B_offsurf(self, B_total: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], Xt: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], max_Nt: typing.SupportsInt | typing.SupportsIndex = -1, max_Np: typing.SupportsInt | typing.SupportsIndex = -1) -> list[float]:
        """
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
        """
    def setup(self, digits: typing.SupportsInt | typing.SupportsIndex, NFP: typing.SupportsInt | typing.SupportsIndex, half_period: bool, Nt: typing.SupportsInt | typing.SupportsIndex, Np: typing.SupportsInt | typing.SupportsIndex, X: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], src_Nt: typing.SupportsInt | typing.SupportsIndex, src_Np: typing.SupportsInt | typing.SupportsIndex, trg_Nt: typing.SupportsInt | typing.SupportsIndex = -1, trg_Np: typing.SupportsInt | typing.SupportsIndex = -1) -> None:
        """
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
        """
class VirtualCasingTestData:
    """
    Generate test data for class VirtualCasing.
    """
    @staticmethod
    def magnetic_field_data(NFP: typing.SupportsInt | typing.SupportsIndex, half_period: bool, Nt: typing.SupportsInt | typing.SupportsIndex, Np: typing.SupportsInt | typing.SupportsIndex, X: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], trg_Nt: typing.SupportsInt | typing.SupportsIndex, trg_Np: typing.SupportsInt | typing.SupportsIndex) -> tuple[list[float], list[float]]:
        """
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
        """
    @staticmethod
    def magnetic_field_grad_data(NFP: typing.SupportsInt | typing.SupportsIndex, half_period: bool, Nt: typing.SupportsInt | typing.SupportsIndex, Np: typing.SupportsInt | typing.SupportsIndex, X: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], trg_Nt: typing.SupportsInt | typing.SupportsIndex, trg_Np: typing.SupportsInt | typing.SupportsIndex) -> tuple[list[float], list[float]]:
        """
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
        """
    @staticmethod
    def surface_coordinates(NFP: typing.SupportsInt | typing.SupportsIndex, half_period: bool, Nt: typing.SupportsInt | typing.SupportsIndex, Np: typing.SupportsInt | typing.SupportsIndex, surf_type: SurfType = ...) -> list[float]:
        """
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
        """
    def __init__(self) -> None:
        ...
