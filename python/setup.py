from distutils.core import setup, Extension

module1 = Extension('call',
                    include_dirs = ['E:/COptdev/include/Base','E:/COptdev/include/ThirParty','E:/COptdev/include/Algorithms','E:/COptdev/include/FunctionRepository','E:/COptdev/include','E:/libs/eigen'],
                    libraries = ['cblas','blas'],
                    library_dirs = ['E:/libs'],
                    sources = ['call.cpp'])

setup (name = 'PackageName',
       version = '1.0',
       description = 'This is a demo package',
       long_description = '''
This is really just a demo package.
''',
       ext_modules = [module1])