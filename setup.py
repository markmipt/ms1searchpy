#!/usr/bin/env python

'''
setup.py file for ms1searchpy
'''

from setuptools import setup, find_packages

version = open('VERSION').readline().strip()

setup(
    name                 = 'ms1searchpy',
    version              = version,
    description          = '''A proteomics search engine for LC-MS1 spectra.''',
    long_description     = (''.join(open('README.md').readlines())),
    author               = 'Mark Ivanov',
    author_email         = 'pyteomics@googlegroups.com',
    install_requires     = ['pyteomics[XML]', 'scipy', 'numpy', 'sklearn'],
    classifiers          = ['Intended Audience :: Science/Research',
                            'Programming Language :: Python :: 2.7',
                            'Topic :: Education',
                            'Topic :: Scientific/Engineering :: Bio-Informatics',
                            'Topic :: Scientific/Engineering :: Chemistry',
                            'Topic :: Scientific/Engineering :: Physics'],
    license              = 'License :: OSI Approved :: Apache Software License',
    packages             = find_packages(),
    scripts              = ['search.py']
    )