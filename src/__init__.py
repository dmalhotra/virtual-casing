"""Virtual casing principle for magnetic field computation.

The extension is compiled once per x86-64 instruction-set level; the variant
matching this CPU is selected here. Set VIRTUAL_CASING_ISA to baseline, v3 or
v4 to override (useful for benchmarking and for testing a variant that this
CPU could not otherwise reach).
"""

import importlib
import os

_ORDER = ("v4", "v3", "baseline")


def _available():
    """Variants this CPU can execute, most capable first."""
    try:
        from . import _cpu
    except ImportError:  # single-variant build (e.g. arm64)
        return ["baseline"]
    ok = ["baseline"]
    if _cpu.has_avx2():
        ok.insert(0, "v3")
    if _cpu.has_avx512():
        ok.insert(0, "v4")
    return ok


def _select():
    forced = os.environ.get("VIRTUAL_CASING_ISA")
    candidates = [forced] if forced else _available()
    if forced and forced not in _ORDER:
        raise ValueError(
            "VIRTUAL_CASING_ISA must be one of %s, got %r" % (", ".join(_ORDER), forced)
        )
    for name in candidates:
        try:
            return name, importlib.import_module("._vc_" + name, __name__)
        except ImportError:
            continue  # variant not built for this platform
    raise ImportError("no usable virtual_casing extension variant was found")


__isa__, _impl = _select()

globals().update({k: v for k, v in vars(_impl).items() if not k.startswith("_")})
