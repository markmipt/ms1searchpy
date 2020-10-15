from sys import argv
from . import main, utils
import argparse

def run():
    parser = argparse.ArgumentParser(
        description='Search proteins using LC-MS spectra',
        epilog='''

    Example usage
    -------------
    $ search.py input.mzML -d human.fasta
    -------------
    ''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', help='input mzML or .tsv file with peptide features')
    parser.add_argument('-d', help='path to protein fasta file', required=True)
    parser.add_argument('-ptol', help='precursor mass tolerance in ppm', default=10.0, type=float)
    parser.add_argument('-rtt', help='retention time threshold in sigma', default=2.0, type=float)
    parser.add_argument('-fdr', help='protein fdr filter in %%', default=1.0, type=float)
    parser.add_argument('-i', help='minimum number of isotopes', default=2, type=int)
    parser.add_argument('-ci', help='minimum number of isotopes for mass and RT calibration', default=4, type=int)
    parser.add_argument('-ts', help='Two-stage RT training: 0 - turn off, 1 - turn one, 2 - turn on and use additive model in the first stage (Default)', default=2, type=int)
    parser.add_argument('-sc', help='minimum number of scans for peptide feature', default=3, type=int)
    parser.add_argument('-lmin', help='min length of peptides', default=7, type=int)
    parser.add_argument('-lmax', help='max length of peptides', default=30, type=int)
    parser.add_argument('-e', help='cleavage rule in quotes!. X!Tandem style for cleavage rules', default='[RK]|{P}')
    parser.add_argument('-mc', help='number of missed cleavages', default=0, type=int)
    parser.add_argument('-cmin', help='min precursor charge', default=1, type=int)
    parser.add_argument('-cmax', help='max precursor charge', default=5, type=int)
    parser.add_argument('-fmods', help='fixed modifications. in mass1@aminoacid1,mass2@aminoacid2 format', default='57.021464@C')
    parser.add_argument('-ad', help='add decoy', default=0, type=int)
    parser.add_argument('-ml', help='use machine learning for PFMs', default=1, type=int)
    parser.add_argument('-prefix', help='decoy prefix', default='DECOY_')
    parser.add_argument('-nproc',   help='number of processes', default=1, type=int)
    parser.add_argument('-elude', help='path to elude binary file. If empty, the built-in additive model will be used for RT prediction', default='')
    parser.add_argument('-deeplc', help='path to deeplc', default='')
    parser.add_argument('-pl', help='path to list of peptides for RT calibration', default='')
    args = vars(parser.parse_args())

    main.process_file(args)
    print('The search is finished.')

if __name__ == '__main__':
    run()
