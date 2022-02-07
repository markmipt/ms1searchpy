ms1searchpy - a DirectMS1 proteomics search engine for LC-MS1 spectra
=====================================================================

`ms1searchpy` consumes LC-MS data (**mzML**) or peptide features (**tsv**) and performs protein identification and quantitation.

Basic usage
-----------
Basic command for protein identification:

    ms1searchpy *.mzML -d path_to.FASTA

OR

    ms1searchpy *_peptideFeatures.tsv -d path_to.FASTA

Read further for detailed info, including quantitative analysis.

Citing ms1searchpy
------------------
Ivanov et al. Boosting MS1-only Proteomics with Machine Learning Allows 2000 Protein Identifications in Single-Shot Human Proteome Analysis Using 5 min HPLC Gradient. https://doi.org/10.1021/acs.jproteome.0c00863

Ivanov et al. DirectMS1: MS/MS-free identification of 1000 proteins of cellular proteomes in 5 minutes. https://doi.org/10.1021/acs.analchem.9b05095

Installation
------------

Using pip:

    pip install ms1searchpy

It is recommended to additionally install [DeepLC](https://github.com/compomics/DeepLC); you may also want to install
[diffacto](https://github.com/statisticalbiotechnology/diffacto):

    pip install deeplc diffacto

This should work on recent versions of Python (3.8-3.10).

Usage tutorial: protein identification
--------------------------------------

The script used for protein identification is called `ms1searchpy`. It needs input files (mzML or tsv) and a FASTA database.

Input files
...........

If mzML are provided, ms1searchpy will invoke [biosaur2](https://github.com/markmipt/biosaur2) to generate the features table.
You can also use other software like [Dinosaur](https://github.com/fickludd/dinosaur) or [Biosaur](https://github.com/abdrakhimov1/Biosaur),
but [biosaur2](https://github.com/markmipt/biosaur2) is recommended. You can also make it yourself,
the table must contain columns 'massCalib', 'rtApex', 'charge' and 'nIsotopes' columns.

**How to get mzML files**:

To get mzML from RAW files, you can use [Proteowizard MSConvert](https://proteowizard.sourceforge.io/download.html)...

    msconvert path_to_file.raw -o path_to_output_folder --mzML --filter "peakPicking true 1-" --filter "MS2Deisotope" --filter "zeroSamples removeExtra" --filter "threshold absolute 1 most-intense"

...or [compomics ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser), which produces suitable files
with default parameters.

RT predictor
............

For protein identification, `ms1searchpy` needs a retention time prediction model. The recommended one is [DeepLC](https://github.com/compomics/DeepLC),
but you can also use the Elude predictor from [Percolator](https://github.com/percolator/percolator) or the built-in additive model (default).

Examples
........

    ms1searchpy test.mzML -d sprot_human.fasta -deeplc deeplc -ad 1

This command will run `ms1searchpy` with DeepLC RT predictor available as `deeplc` (should work if you install DeepLC
alongside `ms1searchpy`. `-ad 1` creates a shuffled decoy database for FDR estimation.
You should use it only once and just use the created database for other searches.

    ms1searchpy test.features.tsv -d sprot_human_shuffled.fasta -deeplc env_deeplc/bin/deeplc

Here, instead of mzML file, a file with peptide features is used. Also, DeepLC is installed in a separate environment, so
a path is specified.

Output files
............

`ms1searchpy` produces several tables:
 - findetified proteins, FDR-filtered (`sample_proteins.tsv`);
 - all identified proteins (`sample_proteins_full.tsv`) - this is the main result;
 - all identified proteins based on all PFMs (`sample_proteins_full_noexclusion.tsv`);
 - all matched peptide match fingerprints, or peptide-feature matches (`sample_PFMs.tsv`);
 - all PFMs with features prepared for Machnine Learning (`sample_PFMs_ML.tsv`);
 - log file with estimated mass and RT accuracies (`sample_log.txt`).

Combine results from replicates
...............................

You can combine the results from several replicate runs with `ms1combine` by feeding it `_PFMs_ML.tsv` tables:

    ms1combine sample_rep_*_PFMs_ML.tsv

Usage tutorial: Quantitation
----------------------------

After obtaining the protein identification results, you can proceed to compare your samples using LFQ.

Using diffacto
..............

Here's an example where we use Bourne Shell syntax for brevity. Each sample contains three replicates:

    ms1todiffacto -dif diffacto -S1 sample1_r{1,2,3}.proteins.tsv -S2 sample2_r{1,2,3}.proteins.tsv -norm median -out diffacto_output.tsv -min_samples 3

`ms1todiffacto` prepares input file for [diffacto](https://github.com/statisticalbiotechnology/diffacto) from ms1searchpy output and to automatically runs diffacto.


Links
-----

- GitHub repo & issue tracker: https://github.com/markmipt/ms1searchpy
- Mailing list: markmipt@gmail.com

- Diffacto repo: https://github.com/statisticalbiotechnology/diffacto
- DeepLC repo: https://github.com/compomics/DeepLC
