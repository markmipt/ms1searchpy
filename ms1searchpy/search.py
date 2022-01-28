from sys import argv
from . import main, utils
import argparse

def run():
    parser = argparse.ArgumentParser(
        description='Search proteins using LC-MS spectra',
        epilog='''

    Example usage
    -------------
    $ search.py input.mzML input2.mzML -d human.fasta -ad 1 -fdr 5.0
    -------------
    ''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('files', help='input mzML or .tsv files with peptide features', nargs='+')
    parser.add_argument('-d', help='path to protein fasta file', required=True)
    parser.add_argument('-ptol', help='precursor mass tolerance in ppm', default=10.0, type=float)
    parser.add_argument('-fdr', help='protein fdr filter in %%', default=1.0, type=float)
    parser.add_argument('-i', help='minimum number of isotopes', default=2, type=int)
    parser.add_argument('-ci', help='minimum number of isotopes for mass and RT calibration', default=4, type=int)
    parser.add_argument('-ts', help='Two-stage RT training: 0 - turn off, 1 - turn one, 2 - turn on and use additive model in the first stage (Default)', default=2, type=int)
    parser.add_argument('-sc', help='minimum number of scans for peptide feature', default=2, type=int)
    parser.add_argument('-lmin', help='min length of peptides', default=7, type=int)
    parser.add_argument('-lmax', help='max length of peptides', default=30, type=int)
    parser.add_argument('-e', help='cleavage rule in quotes!. X!Tandem style for cleavage rules: "[RK]|{P}" for trypsin,\
     "[X]|[D]" for asp-n or "[RK]|{P},[K]|[X]" for mix of trypsin and lys-c', default='[RK]|{P}')
    parser.add_argument('-mc', help='number of missed cleavages', default=0, type=int)
    parser.add_argument('-cmin', help='min precursor charge', default=1, type=int)
    parser.add_argument('-cmax', help='max precursor charge', default=4, type=int)
    parser.add_argument('-fmods', help='fixed modifications. in mass1@aminoacid1,mass2@aminoacid2 format', default='57.021464@C')
    parser.add_argument('-ad', help='add decoy', default=0, type=int)
    parser.add_argument('-ml', help='use machine learning for PFMs', default=1, type=int)
    parser.add_argument('-prefix', help='decoy prefix', default='DECOY_')
    parser.add_argument('-nproc',   help='number of processes', default=1, type=int)
    parser.add_argument('-elude', help='path to elude binary file. If empty, the built-in additive model will be used for RT prediction', default='')
    parser.add_argument('-deeplc', help='path to deeplc', default='')
    parser.add_argument('-deeplc_library', help='path to deeplc library', default='')
    parser.add_argument('-pl', help='path to list of peptides for RT calibration', default='')
    parser.add_argument('-mcalib',   help='mass calibration: 2 - group by ion mobility and RT, 1 - by RT, 0 - no calibration', default=2, type=int)
    args = vars(parser.parse_args())

    main.process_file(args)

if __name__ == '__main__':
    run()
