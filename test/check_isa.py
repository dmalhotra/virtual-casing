"""Check that every ISA variant which was built AND is runnable here can be selected.

The set of variants depends on the target: arm64 builds only baseline, and a CPU
without AVX-512 must not be asked to load v4 (that would be an illegal
instruction, not an exception).
"""

import importlib.util
import os
import subprocess
import sys

import virtual_casing as vc

ALL = ("baseline", "v3", "v4")

built = {n for n in ALL if importlib.util.find_spec("virtual_casing._vc_" + n) is not None}

runnable = {"baseline"}
try:
    from virtual_casing import _cpu

    if _cpu.has_avx2():
        runnable.add("v3")
    if _cpu.has_avx512():
        runnable.add("v4")
except ImportError:
    pass

print("auto-selected: %s | built: %s | runnable: %s"
      % (vc.__isa__, sorted(built), sorted(runnable)))

assert vc.__isa__ in built, "selected a variant that was not built"

for name in sorted(built & runnable):
    env = dict(os.environ, VIRTUAL_CASING_ISA=name)
    r = subprocess.run([sys.executable, "-c", "import virtual_casing as vc; print(vc.__isa__)"],
                       env=env, capture_output=True, text=True)
    assert r.returncode == 0 and r.stdout.strip() == name, \
        "forcing %s failed: rc=%d out=%r err=%r" % (name, r.returncode, r.stdout, r.stderr)
    print("  forced %s: OK" % name)

print("ISA dispatch OK")
