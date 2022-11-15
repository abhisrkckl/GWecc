from setuptools import setup, Extension
import os

GWecc_sources = (
    "antenna_pattern.cpp",
    "EccentricResiduals_Anl.cpp",
    "EccentricResiduals_Num.cpp",
    "EccentricResiduals_Adb.cpp",
    "EccentricResiduals_PM.cpp",
    "EccentricResiduals.cpp",
    "GWecc.cpp",
    "Evolve.cpp",
    "FourierWaveform.cpp",
    "mikkola.cpp",
    "OrbitalEvolution.cpp",
    "PN.cpp",
    "Precompute_Orbit.cpp",
    "EccentricWaveform.cpp",
    "FeStat.cpp",
    "FeStat_Adb.cpp",
)
src_dir = "GWecc/"
current_dir = os.path.dirname(os.path.realpath(__file__))
include_dir = current_dir + "/" + "GWecc/"
GWecc_sources = [src_dir + src_file for src_file in GWecc_sources] + [
    "enterprise_GWecc/GWecc.i"
]

GWecc_cpp_module = Extension(
    "enterprise_GWecc._GWecc",
    sources=GWecc_sources,
    include_dirs=[include_dir, current_dir],
    libraries=["gsl", "gslcblas"],
    swig_opts=["-c++", "-threads"],
    extra_compile_args=["-std=c++17", "-Wno-unused-result"],
)

setup(
    name="enterprise_GWecc",
    version="0.1.5",
    description="Computes pulsar TOA delays due to gravitational waves from eccentric supermassive binary sources.",
    author="Abhimanyu Susobhanan",
    author_email="abhisrkckl@gmail.com",
    ext_modules=[GWecc_cpp_module],
    py_modules=[
        "enterprise_GWecc.GWecc",
        "enterprise_GWecc.enterprise_GWecc",
        "enterprise_GWecc.enterprise_GWecc_cosmoz",
        "enterprise_GWecc.block",
        "enterprise_GWecc.spline",
        "enterprise_GWecc.Fe_stat",
    ],
    # scripts = ['examples/Example1.py'],
    install_requires=["numpy", "astropy", "matplotlib", "enterprise-pulsar"],
)
