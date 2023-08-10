# ms1searchpy - a DirectMS1 proteomics search engine for LC-MS1 spectra

`ms1searchpy` consumes LC-MS data (**mzML**) or peptide features (**tsv**) and performs protein identification and quantitation.

## Basic usage

Basic command for protein identification:

    ms1searchpy *.mzML -d path_to.FASTA

or

    ms1searchpy *_peptideFeatures.tsv -d path_to.FASTA

Read further for detailed info, including quantitative analysis.

## Citing ms1searchpy

Ivanov et al. DirectMS1Quant: Ultrafast Quantitative Proteomics with MS/MS-Free Mass Spectrometry. https://pubs.acs.org/doi/10.1021/acs.analchem.2c02255

Ivanov et al. Boosting MS1-only Proteomics with Machine Learning Allows 2000 Protein Identifications in Single-Shot Human Proteome Analysis Using 5 min HPLC Gradient. https://doi.org/10.1021/acs.jproteome.0c00863

Ivanov et al. DirectMS1: MS/MS-free identification of 1000 proteins of cellular proteomes in 5 minutes. https://doi.org/10.1021/acs.analchem.9b05095

## Installation

Using pip:

    pip install ms1searchpy


## Usage tutorial: protein identification

The script used for protein identification is called `ms1searchpy`. It needs input files (mzML or tsv) and a FASTA database.

### Input files

If mzML are provided, ms1searchpy will invoke [biosaur2](https://github.com/markmipt/biosaur2) to generate the features table.
You can also use other software like [Dinosaur](https://github.com/fickludd/dinosaur) or [Biosaur](https://github.com/abdrakhimov1/Biosaur),
but [biosaur2](https://github.com/markmipt/biosaur2) is recommended. You can also make it yourself,
the table must contain columns 'massCalib', 'rtApex', 'charge' and 'nIsotopes' columns.

#### How to get mzML files

To get mzML from RAW files, you can use [Proteowizard MSConvert](https://proteowizard.sourceforge.io/download.html)...

    msconvert path_to_file.raw -o path_to_output_folder --mzML --filter "peakPicking true 1-" --filter "MS2Deisotope" --filter "zeroSamples removeExtra" --filter "threshold absolute 1 most-intense"

...or [compomics ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser), which produces suitable files
with default parameters.

### RT predictor

For protein identification, `ms1searchpy` needs a retention time prediction model. The recommended one is [DeepLC](https://github.com/compomics/DeepLC)  (default),
but you can also use built-in additive model.

### Examples

    ms1searchpy test.mzML -d sprot_human.fasta -deeplc 1 -ad 1

This command will run `ms1searchpy` with DeepLC RT predictor available as `deeplc` . `-ad 1` creates a shuffled decoy database for FDR estimation.
You should use it only once and just use the created database for other searches.

    ms1searchpy test.features.tsv -d sprot_human_shuffled.fasta -deeplc 1

Here, instead of mzML file, a file with peptide features is used. Also, DeepLC is installed in a separate environment, so
a path is specified.

### Output files

`ms1searchpy` produces several tables:
 - findetified proteins, FDR-filtered (`sample.features_proteins.tsv`);
 - all identified proteins (`sample.features_proteins_full.tsv`) - this is the main result;
 - all identified proteins based on all PFMs (`sample.features_proteins_full_noexclusion.tsv`);
 - all matched peptide match fingerprints, or peptide-feature matches (`sample.features_PFMs.tsv`);
 - all PFMs with features prepared for Machnine Learning (`sample.features_PFMs_ML.tsv`);
 - number of theoretical peptides per protein (`sample.features_protsN.tsv`);
 - log file with estimated mass and RT accuracies (`sample.features_log.txt`).

### Combine results from replicates

You can combine the results from several replicate runs with `ms1combine` by feeding it `_PFMs_ML.tsv` tables:

    ms1combine sample_rep_*.features_PFMs_ML.tsv

## Usage tutorial: Quantitation

After obtaining the protein identification results, you can proceed to compare your samples using LFQ.

### Using directms1quant

New LFQ method designed specifically for DirectMS1 is invoked like this:

    directms1quant -S1 sample1_r{1,2,3}.features_proteins_full.tsv -S2 sample2_r{1,2,3}.features_proteins_full.tsv

It produces a filtered table of significantly changed proteins with p-values and fold changes,
as well as the full protein table and a separate file simply listing all
IDs of significantly modified proteins (e.g. for easy copy-paste into a StringDB search window).

## Links

- GitHub repo & issue tracker: https://github.com/markmipt/ms1searchpy
- Mailing list: markmipt@gmail.com

- DeepLC repo: https://github.com/compomics/DeepLC
