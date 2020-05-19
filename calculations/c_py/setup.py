from distutils.core import setup, Extension

setup(name='covid_genomics', ext_modules=[Extension('covid_genomics', extra_compile_args=['-O3', '-std=c++17', '-Wall', '-Werror', '-Wextra'], sources=['covid_genomics.cpp'])])

