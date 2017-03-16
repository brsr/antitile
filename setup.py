#! /usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

setup(name='antitile',
      version='20170316',
      description='Tiling subdivision, for use with Antiprism',
      url='http://github.com/brsr/antitile',
      author='B R S Recht',
      author_email='brstone@ufl.edu',
      license='MIT',
      packages=['antitile'],
      keywords="tiling tesselation geodesic grid goldberg polyhedra " +
               "sphere polyhedron geometry buckyball",
      install_requires=['numpy', 'scipy', 'matplotlib'],
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Topic :: Artistic Software',
        'Topic :: Multimedia :: Graphics :: 3D Modeling',
        'Topic :: Scientific/Engineering :: Mathematics'
      ],
      scripts = ['bin/view_off.py', 'bin/balloon.py', 'bin/subtile.py'],
      include_package_data=True,
      zip_safe=True,
      test_suite='tests')
