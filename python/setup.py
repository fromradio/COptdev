from distutils.core import setup, Extension

module1 = Extension('geVec',
                    include_dirs = ['E:/COptdev/include/Base','E:/COptdev/include/ThirParty','E:/COptdev/Algorithms','E:/COptdev/Differential','E:/COptdev/FunctionRepository','E:/COptdev/LeastSquares','E:/COptdev/include','E:/libs/Eigen'],
                    libraries = ['cblas','blas'],
                    library_dirs = ['E:/libs'],
                    sources = ['geVec.cpp'])

setup (name = 'PackageName',
       version = '1.0',
       description = 'This is a demo package',
       long_description = '''
This is really just a demo package.
''',
       ext_modules = [module1])