#!/usr/bin/env python

'''
setup.py file for ms1searchpy
'''
import os
from setuptools import setup, find_packages, Extension

version = open('VERSION').readline().strip()

def make_extensions():
    is_ci = bool(os.getenv("CI", ""))
    include_diagnostics = False
    try:
        import numpy
    except ImportError:
        print("C Extensions require `numpy`")
        raise
    from Cython.Build import cythonize
    cython_directives = {
        'embedsignature': True,
        "profile": include_diagnostics
    }
    macros = []
    if include_diagnostics:
        macros.append(("CYTHON_TRACE_NOGIL", "1"))
    if is_ci and include_diagnostics:
        cython_directives['linetrace'] = True

    extensions = cythonize([
        Extension(name='ms1searchpy.cutils', sources=['ms1searchpy/cutils.pyx'],
                    include_dirs=[numpy.get_include(), ])
    ], compiler_directives=cython_directives)
    return extensions


def do_setup(cext=True):
    setup(
        name                 = 'ms1searchpy',
        version              = version,
        description          = '''A proteomics search engine for LC-MS1 spectra.''',
        long_description     = (''.join(open('README.md').readlines())),
        author               = 'Mark Ivanov',
        author_email         = 'pyteomics@googlegroups.com',
        install_requires     = ['pyteomics[XML]', 'scipy', 'numpy', 'sklearn', 'matplotlib', 'pandas', 'seaborn'],
        ext_modules=make_extensions() if cext else None,
        classifiers          = ['Intended Audience :: Science/Research',
                                'Programming Language :: Python :: 2.7',
                                'Topic :: Education',
                                'Topic :: Scientific/Engineering :: Bio-Informatics',
                                'Topic :: Scientific/Engineering :: Chemistry',
                                'Topic :: Scientific/Engineering :: Physics'],
        license              = 'License :: OSI Approved :: Apache Software License',
        packages             = find_packages(),
        package_data         = {'ms1searchpy': ['Dinosaur/*']},
        entry_points         = {'console_scripts': ['ms1searchpy = ms1searchpy.search:run', ]}
        )

try:
    do_setup(True)
except Exception as err:
    print("*" * 60)
    print("Could not compile C Extensions due to %r, attempting pure Python installation." % (err,))
    print("*" * 60)
    do_setup(False)
    print("Could not compile C Extensions due to %r, speedups are not enabled." % (err,))