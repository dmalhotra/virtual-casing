// Baseline-compiled probe used by virtual_casing/__init__.py to choose an
// instruction-set variant. __builtin_cpu_supports checks OS enablement via
// XGETBV, not just CPUID, so a kernel that has not enabled AVX-512 reports false.
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace {

bool has_avx2() {
#if defined(__x86_64__) || defined(__i386__)
  __builtin_cpu_init();
  return __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
#else
  return false;
#endif
}

bool has_avx512() { // the x86-64-v4 subset
#if defined(__x86_64__) || defined(__i386__)
  __builtin_cpu_init();
  return __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512bw")
      && __builtin_cpu_supports("avx512cd") && __builtin_cpu_supports("avx512dq")
      && __builtin_cpu_supports("avx512vl");
#else
  return false;
#endif
}

} // namespace

PYBIND11_MODULE(_cpu, m) {
  m.doc() = "CPU feature probe for instruction-set dispatch.";
  m.def("has_avx2", &has_avx2);
  m.def("has_avx512", &has_avx512);
}
