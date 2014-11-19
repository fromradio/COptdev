from distutils.core import setup, Extension

module1 = Extension('copt',
                    include_dirs = ['../../include'],
                    libraries = ['cblas','blas'],
                    library_dirs = [],
                    sources = ['py_least_squares_wrapper.cpp'])

setup (name = 'COPT',
       version = '1.0',
       description = 'Least Squares Method',
       long_description = '''
The module contains three basic approaches for solving least square method.
''',
       ext_modules = [module1])