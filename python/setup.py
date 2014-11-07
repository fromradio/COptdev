#! /usr/bin/python 



from distutils.core import setup, Extension 



module1 = Extension('geVec', sources = ['geVec.cpp']) 



setup (name = 'PackageName', version = '2.6', description = 'This is a demo package', ext_modules = [module1]) 
