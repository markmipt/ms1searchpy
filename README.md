ms1searchpy - a proteomics search engine for MS1 spectra
---------------------------------------------------------------

The mzML, FASTA and .cfg files are required for basic operation of the script.
Cfg file contains settings for the algorithm. For an efficient usage of retention time, retention coefficients for additive model training are required. They can be obtained using either MPscore software (https://bitbucket.org/markmipt/mp-score) on MS/MS data or script getRC.py included here.

Algorithm can be run with following command:

    python search.py path_to_mzML path_to_cfg

The script output contains files:
    all identified proteins (filename_proteins_full.csv)
    filtered proteins (filename_proteins.csv)
    all matched peptide match fingerprints (filename_PFMs.csv)

To combine results of MS1 searches use the following command:

    python search.py path_to_output_proteins_full.csv path_to_cfg

Dependencies
------------

- pyteomics
- numpy
- scipy

Links
-----

- BitBucket repo & issue tracker: https://bitbucket.org/markmipt/ms1searchpy
- Mailing list: pyteomics@googlegroups.com