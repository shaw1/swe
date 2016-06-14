from numpy.distutils.core import Extension, setup

swe_topog_sources = ["fortran/ftcs_topog.f95", \
    "fortran/ftcs_topog_tlm.f95", "fortran/ftcs_topog_adj.f95"]

swe_flat_sources = ["fortran/ftcs_flat.f95", \
    "fortran/ftcs_flat_tlm.f95", "fortran/ftcs_flat_adj.f95"]

setup(name="swe_topog_fortran", ext_modules=[Extension(name="swe_topog_fortran",
    sources=swe_topog_sources)])

setup(name="swe_flat_fortran", ext_modules=[Extension(name="swe_flat_fortran",
    sources=swe_flat_sources)])
