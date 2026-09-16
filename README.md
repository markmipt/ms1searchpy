# ms1searchpy - a DirectMS1 proteomics search engine for LC-MS1 spectra

`ms1searchpy` consumes peptide isotopic envelopes from LC-MS1 data (**tsv**) and performs protein identification and quantitation. It is recommended to use biosaur2 for the peptide ion envelopes detection.

## Basic usage

Basic command for protein identification:

    biosaur2 filename1.mzML
    ms1searchpy filename1.features.tsv -d path_to_fasta_with_decoys.FASTA -deeplc 1
    
To use msallsearchpy functionality add path to mzML file

    ms1searchpy filename1.features.tsv -d path_to_fasta_with_decoys.FASTA -deeplc 1 -ms2mzml filename1.mzML

To create peptide-level shuffled decoy database add "-ad 1" option:

    ms1searchpy filename1.features.tsv -d path_to_fasta_without_decoys.FASTA -ad 1 -deeplc 1

To speed up RT predictions for subsequent searches, save the RT prediction results to a file:

    ms1searchpy filename1.features.tsv -d path_to_fasta_with_decoys.FASTA -deeplc 1 -deeplc_library /home/user1/deeplc_lib_path.lib


Read further for detailed info, including quantitative analysis.

## Citing ms1searchpy

Ivanov et al. DirectMS1Quant: Ultrafast Quantitative Proteomics with MS/MS-Free Mass Spectrometry. https://pubs.acs.org/doi/10.1021/acs.analchem.2c02255

Ivanov et al. Boosting MS1-only Proteomics with Machine Learning Allows 2000 Protein Identifications in Single-Shot Human Proteome Analysis Using 5 min HPLC Gradient. https://doi.org/10.1021/acs.jproteome.0c00863

Ivanov et al. DirectMS1: MS/MS-free identification of 1000 proteins of cellular proteomes in 5 minutes. https://doi.org/10.1021/acs.analchem.9b05095

## Installation

We suggest to use a standalone virtual enviroment with Python version 3.10.11. An example on how to create and to
activate such enviroment is shown below:

    pyenv install 3.10.11
    pyenv virtualenv 3.10.11 virt_ms1searchpy
    pyenv activate virt_ms1searchpy

After that, install ms1searchpy:

    pip install ms1searchpy

It will automatically install unofficial fork of [DeepLC](https://github.com/compomics/DeepLC), as well as [Identipy](https://github.com/levitsky/identipy) search engine. If you need an option to use [MS2PIP](https://github.com/compomics/ms2pip) for MS/MS spectra processing, install ms2pip:

    pip install ms2pip==4.2.0

## Usage tutorial: protein identification

The script used for protein identification is called `ms1searchpy`. It needs input files (tsv) and a FASTA database.

### Input files

You need a file with peptide ion envelopes. The default way is to use [biosaur2](https://github.com/markmipt/biosaur2) to generate the features table. You can also use other software like [Dinosaur](https://github.com/fickludd/dinosaur) or [Biosaur](https://github.com/abdrakhimov1/Biosaur), but [biosaur2](https://github.com/markmipt/biosaur2) is recommended. You can also make it yourself, the table must contain columns 'massCalib', 'rtApex', 'charge' and 'nIsotopes' columns.
All of the mentioned feature detection algorithms work with mzML files as input.

#### How to get mzML files

To get mzML from RAW files, you can use [Proteowizard MSConvert](https://proteowizard.sourceforge.io/download.html)...

    msconvert path_to_file.raw -o path_to_output_folder --mzML --filter "peakPicking true 1-" --filter "MS2Deisotope" --filter "zeroSamples removeExtra" --filter "threshold absolute 1 most-intense"

...or [compomics ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser), which produces suitable files
with default parameters.

### RT predictor

For protein identification, `ms1searchpy` needs a retention time prediction model. The recommended one is [DeepLC](https://github.com/compomics/DeepLC),
but you can also use built-in additive model (default or with an option "-deeplc 0").

### Examples

    biosaur2 test.mzML -minlh 3
    ms1searchpy test.features.tsv -d sprot_human.fasta -deeplc 1 -ad 1

The first command will run `biosaur2` to detect all peptide isotopic clusters which are visible in at least 3 consecutive MS1 scans. The second command will run `ms1searchpy` with DeepLC RT predictor. `-ad 1` creates a shuffled decoy database for FDR estimation. You should use it only once and use the created database for other searches within the project to
proceed quantitative analysis.

### Output files

`ms1searchpy` produces several tables:
 - identified proteins, FDR-filtered (`sample.features_proteins.tsv`) - this is the main result;
 - all peptide-feature matches (PFMs) (`sample.features_PFMs.tsv`);
 - all PFMs with extended columns used for Machnine Learning (`sample.features_PFMs_ML.tsv`);
 - all proteins (`sample.features_proteins_full.tsv`);
 - all proteins with scores based on all PFMs (`sample.features_proteins_full_noexclusion.tsv`);
 - number of theoretical peptides per protein (`sample.features_protsN.tsv`);
 - log file with estimated mass and RT accuracies (`sample.features_log.txt`).

### msallsearchpy (for DDA and DIA data)

The workflow we called msallsearchpy is the basic algorithm of ms1searchpy, where some of the matched PFMs are extended with useful MS/MS-based information. It means that some of the PFMs became more reliable identifications, but the PFMs with no matched fragments are still fully used in the protein identification and quantification process. To use msallsearchpy functionality add path to mzML file:

    ms1searchpy test.features.tsv -d sprot_human.fasta -deeplc 1 -ad 1 -ms2mzml test.mzML

It will automatically detect all MS/MS spectra within isolation window for PFM and calculate different MS/MS-based scores.
See details in Ivanov et al (doi: UNPUBLISHED).


### Using directms1quant

New LFQ method designed specifically for DirectMS1 is invoked like this:

    directms1quant -S1 sample1_r{1,2,3}.features_proteins_full.tsv -S2 sample2_r{1,2,3}.features_proteins_full.tsv

It produces a filtered table of significantly changed proteins with p-values and fold changes,
as well as the full protein table and a separate file simply listing all
IDs of significantly modified proteins (e.g. for easy copy-paste into a StringDB search window).

It was designed to automatically set a fold change threshold and produces the results with well-controlled quantitative
FDR according to our tests against multiple benchmark datasets (LFQ Bench, UPS-E.coli, TPP experiments, etc).

### Combine results from replicates

If you want, you can combine the results from several replicate runs.

The simplest method is to average the protein scores obtained from multiple runs and filter the results again using decoys:

    ms1combine_proteins sample_rep_1.features_proteins_full.tsv sample_rep_2.features_proteins_full.tsv sample_rep_3.features_proteins_full.tsv

The second method involves combining all PFM results following the machine learning stage and recalculating the protein scores:

    ms1combine sample_rep_1.features_PFMs_ML.tsv sample_rep_2.features_PFMs_ML.tsv sample_rep_3.features_PFMs_ML.tsv

### Using Group-specific FDR for metaproteomics

Group-specific FDR for metaproteomics should be used for accurate estimation of protein identified among the different groups. The command ms1groups should be used for that:

     ms1groups F04.features_PFMs_ML.tsv -d F04_top15_shuffled.fasta -out group_statistics_by -fdr 5.0 -groups genus

It produces a table with the number of identified proteins for each group using group-specific FDR. This is basically multiple DirectMS1 searches with small protein databases containing only a single group and combining results all together. However, using the mentioned “ms1groups” script and preliminary DirectMS1 search, two problems are solved: small statistics for all mass/RT/Machine Learning calibration procedures within DirectMS1 workflow for low-populated groups and computational time. Currently supported groups are 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'domain'. The groups are automatically extracted using ete3 Python module and NCBI Taxonomic database. Also, the script supports groups dbname and OX: the the former is a taxonomy in swiss-prot protein name (_HUMAN, _YEAST, etc.) and the latter is the taxonomy provided by 'OX=' from protein description in the fasta file.

## Usage tutorial: Quantitation

After obtaining the protein identification results, you can proceed to compare your samples using LFQ.

### Multi-condition protein profiling using directms1quantmulti

You can make a quantitation for complex projects using script directms1quantmulti. The example below is shown for our project of time-series profiling of glioblastoma cell line under interferon treatment.

Script takes a tab-separated table (.tsv) with details for all project files. An example of a sample file table is available here in the examples folder. It should contain the following columns:

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

    directms1quantmulti -db ~/fasta_folder/sprot_human_shuffled.fasta -pdir ~/folder_with_ms1searchpy_results/ -samples ~/samples.csv -min_signif_for_pept 2 -out DQmulti_2024 -pep_min_non_missing_samples 0.75 -start_stage 1 -proteins_for_figure ~/custom_list_of_proteins.tsv -figdir ~/output_figure_folder/

## Links

- GitHub repo & issue tracker: https://github.com/markmipt/ms1searchpy
- Mailing list: markmipt@gmail.com

- DeepLC repo: https://github.com/compomics/DeepLC
