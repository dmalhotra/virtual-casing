"""Public API of the package.

__init__.py re-exports one instruction-set variant chosen at import time, which
static analysis cannot follow, so the surface is declared here instead. The
variants are identical apart from the instructions they are compiled to, and
_vc_baseline is built on every platform, so it stands in for all of them.
"""

from virtual_casing._vc_baseline import SCTLDoubleVector as SCTLDoubleVector
from virtual_casing._vc_baseline import SCTLIntVector as SCTLIntVector
from virtual_casing._vc_baseline import SurfType as SurfType
from virtual_casing._vc_baseline import VirtualCasing as VirtualCasing
from virtual_casing._vc_baseline import VirtualCasingTestData as VirtualCasingTestData

__isa__: str

__all__: list[str] = [
    "SCTLDoubleVector",
    "SCTLIntVector",
    "SurfType",
    "VirtualCasing",
    "VirtualCasingTestData",
]
