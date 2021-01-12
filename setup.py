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
    long_description_content_type = 'text/markdown',
    author               = 'Mark Ivanov',
    author_email         = 'pyteomics@googlegroups.com',
    install_requires     = ['pyteomics[XML]', 'scipy', 'numpy', 'sklearn', 'matplotlib', 'pandas', 'seaborn', 'lightgbm'],
    classifiers          = ['Intended Audience :: Science/Research',
                            'Programming Language :: Python :: 3.8',
                            'Topic :: Education',
                            'Topic :: Scientific/Engineering :: Bio-Informatics',
                            'Topic :: Scientific/Engineering :: Chemistry',
                            'Topic :: Scientific/Engineering :: Physics'],
    license              = 'License :: OSI Approved :: Apache Software License',
    packages             = find_packages(),
    package_data         = {'ms1searchpy': ['Dinosaur/*']},
    entry_points         = {'console_scripts': ['ms1searchpy = ms1searchpy.search:run',
                                                'ms1combine = ms1searchpy.combine:run',
                                                'ms1tonsaf = ms1searchpy.ms1tonsaf:run',
                                                'ms1todiffacto = ms1searchpy.ms1todiffacto:run',]}
    )
