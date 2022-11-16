from setuptools import setup, Extension
import os

GWecc_sources = (
    "antenna_pattern.cpp",
    "postcircular_residuals_1.cpp",
    "postcircular_residuals_2.cpp",
    "numerical_residuals.cpp",
    "adiabatic_residuals.cpp",
    "peters_mathews_residuals.cpp",
    "eccentric_residuals.cpp",
    "eccentric_waveform.cpp",
    "GWecc.cpp",
    "mikkola.cpp",
    "orbital_evolution_1.cpp",
    "orbital_evolution_2.cpp",
    "post_newtonian.cpp",
    "precompute_orbit.cpp",
    "FeStat.cpp",
    "FeStat_Adb.cpp",
    "waveform_vars.cpp"
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
    extra_compile_args=["-std=c++17", "-Wno-unused-result", "-Werror"],
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
