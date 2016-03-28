from numpy.distutils.core import Extension, setup

swe_sources = ["fortran/ftcs.f95", "fortran/ftcs_tlm.f95", \
    "fortran/ftcs_adj.f95"]

setup(name="swe_fortran", ext_modules=[Extension(name="swe_fortran",
    sources=swe_sources)])
