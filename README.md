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

It is recommended to additionally install [DeepLC](https://github.com/compomics/DeepLC) version either 1.1.2 (official) or 1.1.2.2 (unofficial fork with small changes) . Newer version has some issues right now.

    pip install deeplc==1.1.2

Or

    pip install https://github.com/markmipt/DeepLC/archive/refs/heads/alternative_best_model.zip

This should work on recent versions of Python (3.8-3.10).

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

For protein identification, `ms1searchpy` needs a retention time prediction model. The recommended one is [DeepLC](https://github.com/compomics/DeepLC),
but you can also use built-in additive model (default).

### Examples

    ms1searchpy test.mzML -d sprot_human.fasta -deeplc 1 -ad 1

This command will run `ms1searchpy` with DeepLC RT predictor available as `deeplc` (should work if you install DeepLC
alongside `ms1searchpy`. `-ad 1` creates a shuffled decoy database for FDR estimation.
You should use it only once and just use the created database for other searches.

    ms1searchpy test.features.tsv -d sprot_human_shuffled.fasta -deeplc 1

Here, instead of mzML file, a file with peptide features is used.

### Output files

`ms1searchpy` produces several tables:
 - identified proteins, FDR-filtered (`sample.features_proteins.tsv`) - this is the main result;
 - all identified proteins (`sample.features_proteins_full.tsv`);
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


### Multi-condition protein profiling using directms1quantmulti

You can make a quantitation for complex projects using script directms1quantmulti. The example below is shown for our project of time-series profiling of glioblastoma cell line under interferon treatment.

Script takes a table with details for all project files. An example of a sample file table is available here in the examples folder. It should contain the following columns:

File Name - filename of raw file. For example, “QEHFX_JB_000379”.

group - sample group of file. In our example, there are K (Control group), IFN30 (treatment with 30 units/ml of interferon) and IFN1000 groups. The first group mentioned in the table will be used as control for pairwise directms1quant runs.  

condition - sample subgroup of file. In our example, there are multiple time points after treatment, such as 0h, 30min, 1h, 2h, etc. By default, only the same conditions will be used for pairwise comparisons. For example, IFN30 0h vs K 0h; IFN1000 0h vs K 0h, etc.

vs - column for specific condition comparison. For example, in our case, we did not have control samples at the 30 min time point. Thus, we would like to proceed directms1quant runs for IFN30 30 min vs K 0h; and IFN1000 30 min vs K 0h comparisons. Thus, for the 30 min IFN30 and IFN1000 files we put “0h” in the “vs” column. See example table for details.

replicate - column for replicate number of specific condition and sample group.

BatchMS - column for mass-spectrometry Batch. This parameter is used for extra normalization within a batch.


The script consists of four different stages and you can rerun the script without rerunning previous stages (“-start_stage” option).

Stage 1 is a set of pairwise DirectMS1Quant runs for different interferon treatment conditions versus control samples.

Stage 2 is preparation of peptide LFQ table for all files using the results obtained in the previous step.

Stage 3 is preparation of the protein LFQ table. Only the peptides labeled by DirectMS1Quant as significantly different between samples in at least X pairwise comparisons are used for protein quantitation. The X parameter is controlled by “min_signif_for_pept” option.

Stage 4 is preparation of LFQ profiling figures for proteins specified in the file under “proteins_for_figure” option. The file should be a tsv table with column “dbname” containing protein database names in the swiss-prot format. Any default directms1quant output table with differentially expressed proteins can be used here.


Example of script usage::

    directms1quantmulti -db ~/fasta_folder/sprot_human_shuffled.fasta -pdir ~/folder_with_ms1searchpy_results/ -samples ~/samples.csv -min_signif_for_pept 2 -out DQmulti_2024 -pep_min_non_missing_samples 0.75 -start_stage 1 -proteins_for_figure ~/custom_list_of_proteins.csv -figdir ~/output_figure_folder/

## Links

- GitHub repo & issue tracker: https://github.com/markmipt/ms1searchpy
- Mailing list: markmipt@gmail.com

- DeepLC repo: https://github.com/compomics/DeepLC
