from distutils.core import setup, Extension

module1 = Extension('copt',
                    include_dirs = ['../../include','../../../Eigen','../../../libs/SuiteSparse/include','..'],
                    libraries = ['blas','umfpack','gfortran'],
                    sources = ['py_least_squares_wrapper.cpp'],
                    extra_compile_args=['-std=c++11'])

setup (name = 'COPT',
       version = '1.0',
       description = 'Least Squares Method',
       long_description = '''
The module contains three basic approaches for solving least square method.
''',
       ext_modules = [module1])