# gCNV Analysis Tools

Some analysis tools for the exome GATK-gCNV pipeline. This is NOT the GATK-gCNV
pipeline itself. Most of these workflows rely on post-processing steps internal
to the Talkowski lab and are not applicable to the output of GATK
GermlineCNVCaller. They are also designed to be run on Terra.

There is a WDL workflow for each analysis and an R script that contains the
core analysis code. It is possible to run the analysis just using the R script
if the inputs are formatted correctly.

## Repository Structure

| Path         | Description                       |
| ------------ | --------------------------------- |
| src/scripts/ | core analysis scripts             |
| src/gelpers/ | R package with helper code        |
| wdl/         | WDL workflows                     |
| Dockerfile   | Dockerfile for all the workflows  |
| Makefile     | Makefile with some useful recipes |

## Workflows

There are workflows for the main analysis and some helper workflows to help
construct some of the inputs for the main workflows.

### `MakeSampleBatchMap`

- workflow: `wdl/MakeSampleBatchMap.wdl`

The purpose of this workflow is to create a file that links batch IDs to sample
IDs. The batch ID refers to the ID that is in the final GATK-gCNV callset,
identifying the samples that were processed together using a single model.
Typically it is derived from the sample set IDs on Terra. The workflow should be
run once per sample set. The outputs are a TSV with batch IDs in the first
column (the same ID repeated) and sample IDs in the second column, and the batch
ID as a string.

### `ConcatSampleBatchMaps`

- workflow: `wdl/ConcatSampleBatchMaps.wdl`

This workflow should be run on a sample set set and simply concatenates all the
files produced by `MakeSampleBatchMap` for each of the constituent sample sets.

### `FixDCR`

- workflow: `wdl/FixDCR.wdl`

Modifies an old version of the batch-level denoised copy ratio (dCR) file to be
compatible with these workflows. Probably will not be needed.

### `MergeDcr`

- workflow: `wdl/MergeDcr.wdl`

Merges per-sample dCR files into per-batch dCR files. Probably will be not be
needed.

### Common Inputs

There are some inputs that are commonly used and are described here.

#### `callset`

The callset refers to the final table of CNV calls produced by the GATK-gCNV
pipeline.

#### `pedigree`

Pedigree file based off the GATK PED format. The GATK PED format is described
[here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format).
Briefly, it is a six column TSV describing the family ID, sample ID,
paternal ID, maternal ID, sex, and phenotype, in that order. The sex is 1 for
male and 2 for female. All other values are treated as unknown. The phenotype is
1 for unaffected and 2 for affected. All other values are treated as unknown.

#### `intervals`

This is the file of genomic intervals used to collect read counts. It is a
three column TSV with the chromosome, start, and end of each intervals. All
intervals are 1-based start and inclusive end. The intervals must be disjoint.

#### `dcr_files` and `dcr_indices`

The dCR and dCR index files. These are usually the per-batch ones.

#### `batch_ids`

The IDs used to identify batches of samples that used the same CNV calling
model. These must match the batch IDs in `callset`.

#### `sample_set_ids`

The IDs used for the sample sets on Terra.

#### `sample_batch_map`

A two-column TSV of the format produced by `ConcatSampleBatchMaps`.

#### `runtime_docker`

The URI for the Docker image created from `Dockerfile`. All of the workflows use
the same Docker image.

### *de novo* Calling

Annotate de novo CNVs in a callset.

- workflow: `wdl/AnnotateDeNovo.wdl`
- script: `src/scripts/annotate_denovo_cnv.R`

#### Inputs

##### `recal_freq`

Should the variant site frequencies be recalibrated? The default is to set the
site frequency of every variant to the frequency just among the parents. If
`false`, there must be a column named "sf" in the callset with the site
frequencies.

##### `hq_cols`

A comma-separated list of column names in the input callset that must all be
true (as parsed by `data.table::fread`) for a CNV to be considered high quality.
This is relevant as only high quality variants will be checked for *de novo*
status. The default is `PASS_SAMPLE,PASS_QS`.

##### `max_freq`

The maximum site frequency (as a decimal fraction) a CNV can have. All variants
above this threshold will be ignored. The default is 0.01.

##### `skip_allosomes`

Should CNVs on allosomes be skipped. The default is to not skip. It is not
necessary to set this option just because there aren't any CNVs on allosomes in
the callset.

#### Description

The source code is the best place to look for the de novo calling steps.
Briefly, there is a check for overlap between CNVs in the parents and the child
and then further filtering of those without overlap using the dCR values.

#### Outputs

The input callset file will have a column added with the predicted inheritance
status of each variant: `denovo`, `inherited`, `mosaic_father`, etc.

### Autosomal Ploidy Calling

Predict the ploidy of autosomes.

- workflow: `wdl/CheckPloidy.wdl`
- script (ploidy calling): `src/scripts/check_ploidy.R`
- script (merging ploidy tables): `src/scripts/merge_ploidy_matrices.R`

#### Inputs

##### `contigs`

The list of chromosomes to check for aneuploidies.

##### `is_hg19`

Are the chromosomes in the dCR files named according to the hg19 naming scheme
(i.e. 1, 2, etc).

#### Description

For each chromosome in `contigs`:
1. dCR values are extracted from the dCR files
2. joined across all the batches by interval coordinates
3. grouped into 36 bins
4. the mean sample dCR is computed for each bin and plotted across the
   chromosome
5. the mean sample dCR across all the bins is computed and used as the predicted
   ploidy
6. any sample with predicted ploidy less than 1.5 or greater than 2.5 is
   treated as having an aneuploidy
7. any sample with a predicted aneuploidy has its dCR intervals (the original
   from the dCR files) grouped into 101 bins and those are plotted
8. the table of the predicted ploidy at the chromosome, the plots of ploidy and
   any predicted aneuploidies are collected in the output

The workflow combines the outputs across `contigs`.

#### Outputs

##### `ploidy_plots`

A tar file with the ploidy plots, including any aneuploidies.

##### `ploidy_matrix`

A TSV with the predicted ploidy of every sample/chromosome pair.

##### `aneuploidies`

A TSV with the predicted ploidy of any sample/chromosome pair predicted to be an
aneuploidy.

### Allosomal Ploidy Calling

Predict the ploidy of allosomes. Requires that the intervals used to run
GATK-gCNV contained intervals on the allosomes.

- workflow: `wdl/CheckSex.wdl`
- script: `src/scripts/check_sex.R`

#### Inputs

Nothing unique to this workflow.

#### Description

Predicts the ploidy of the allosomes by computing their mean dCR.

##### Outputs

A TSV of the predicted ploidy of every sample at chrX and chrY.

### CNV Evidence Plotting

Make plots of CNV evidence.

- workflow: `wdl/PlotCNVEvidence.wdl`
- script: `src/scripts/plot_cnv_evidence.R`

#### Inputs

##### `calls_per_split`

Number of CNVs to process per shard.

#### Description

Plots a sample's dCR values for every interval overlapping a CNV.

#### Outputs

A tar file of plots.

### *de novo* CNV Evidence Plotting

Make visualizations of potential *de novo* CNVs.

- workflow: `wdl/PlotDeNovoEvidence.wdl`
- script: `src/scripts/plot_denovo_evidence.R`

#### Inputs

##### `denovos`

GATK-gCNV callset with inheritance status annotated.

##### Description

Plots the dCR values for every interval overlapping a CNV for all members of a
trio.

#### Outputs

A tar file of plots.
