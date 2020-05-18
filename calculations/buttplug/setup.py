from distutils.core import setup, Extension

setup(name='buttplug', ext_modules=[Extension('buttplug', extra_compile_args=['-O3', '-std=c++17', '-Wall', '-Werror', '-Wextra'], sources=['buttplug.cpp'])])

