# MGscan Metagenomics: Usage

## Table of Contents

* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Running the pipeline](#running-the-pipeline)
* [Updating the pipeline](#updating-the-pipeline)
* [Reproducibility](#reproducibility)
* [Main arguments](#main-arguments)
   * [`-profile`](#-profile)
   * [`--design`](#--design)
   * [`--database`](#--database)
* [Trimming and filtering](#trimming-and-filtering)
   * [`--preprocessing_qc_tool`](#preprocessing-qc-tool)
   * [`--save_preprocessed_reads`](#save-preprocessed-reads)
   * [`--shortread_qc_tool`](#shortread-qc-tool)
   * [`--shortread_qc_skipadaptertrim`](#--shortread-qc-skipadaptertrim)
   * [`--shortread_qc_adapter1`](#--shortread-qc-adapter1)
   * [`--shortread_qc_adapter2`](#--shortread-qc-adapter2)
   * [`--shortread_qc_adapterlist`](#--shortread-qc-adapterlist)
   * [`--shortread_qc_mergepairs`](#--shortread-qc-mergepairs)
   * [`--shortread_qc_includeunmerged`](#--shortread-qc-includeunmerged)
   * [`--shortread_qc_minlength`](#--shortread_qc_minlength)
   * [`--shortread_complexityfilter_tool`](#--shortread-complexityfilter-tool)
   * [`--shortread_complexityfilter_entropy`](#--shortread-complexityfilter-entropy)
   * [`--shortread_complexityfilter_bbduk_windowsize`](#--shortread-complexityfilter-bbduk-windowsize)
   * [`--shortread_complexityfilter_bbduk_mask`](#--shortread-complexityfilter-bbduk-mask)
   * [`--shortread_complexityfilter_fastp_threshold`](#--shortread-complexityfilter-fastp-threshold)
   * [`--shortread_complexityfilter_prinseqplusplus_mode`](#--shortread-complexityfilter-prinseqplusplus-mode)
   * [`--shortread_complexityfilter_prinseqplusplus_dustscore`](#--shortread-complexityfilter-prinseqplusplus-dustscore)
   * [`--save_complexityfiltered_reads`](#--save-complexityfiltered-reads)
* [Preprocessing Long Read QC Options](#preprocessing-long-read-qc-options)
   * [`--perform_longread_qc`](#--perform-longread-qc)
   * [`--longread_qc_skipadaptertrim`](#--longread-qc-skipadaptertrim)
   * [`--longread_qc_skipqualityfilter`](#--longread-qc-skipqualityfilter)
   * [`--longread_qc_qualityfilter_minlength`](#--longread-qc-quality-filter-minlength)
   * [`--longread_qc_qualityfilter_keeppercent`](#--longread-qc-qualityfilter-keeppercent)
   * [`--longread_qc_qualityfilter_targetbases`](#--longread-qc-qualityfiler-targetbases)
* [Preprocessing Host Removal Options](#preprocessing-host-removal-options)
   * [`--perform_shortread_hostremoval`](#--perform-shortread-hostremoval)
   * [`--perform_longread_hostremoval`](#--perform-longread-hostremoval)
   * [`--hostremoval_reference`](#--hostremoval-reference)
   * [`--shortread_hostremoval_index`](#--shortread-hostremoval-index)
   * [`--longread_hostremoval_index`](#--longread-hostremoval-index)
   * [`--save_hostremoval_index`](#--save-hostremoval-index)
   * [`--save_hostremoval_bam`](#--save-hostremoval-bam)
   * [`--save_hostremoval_unmapped`](#--save-hostremoval-unmapped)
* [Preprocessing Run Merging Options](#preprocessing-run-merging-options)
   * [`--perform_runmerging`](#--perform-runmerging)
   * [`--save_runmerged_reads`](#--save-runmerged-reads)
* [Profiling Options](#profiling-options)
   * [`--centrifuge_save_reads`](#--centrifuge-save-reads)
   * [`--diamond_output_format`](#--diamond-output-format)
   * [`--diamond_save_reads`](#--diamond-save-reads)
   * [`--kaiju_taxon_rank`](#--kaiju-taxon-rank)
   * [`--kraken2_save_reads`](#--kraken2-save-reads)
   * [`--kraken2_save_readclassification`](#--kraken2-save-readclassification)
   * [`--kraken2_save_minimizers`](#--kraken2-save-minimizers)
   * [`--krakenuniq_save_reads`](#--krakenuniq-save-reads)
   * [`--krakenuniq_ram_chunk_size`](#--krakenuniq-ram-chunk-size)
   * [`--krakenuniq_save_readclassification`](#--krakenuniq-save-readclassification)
   * [`--malt_mode`](#--malt-mode)
   * [`--malt_save_reads`](#--malt-save-reads)
   * [`--malt_generate_megansummary`](#--malt-generate-megansummary)
   * [`--sourmash_kmersize`](#--sourmash-kmersize)
   * [`--sourmash_threshold_bp`](#--sourmash-threshold-bp)
   * [`--sourmash_trim_low_abund`](#--sourmash-trim-low-abund)
   * [`--motus_use_relative_abundance`](#--motus-use-relative-abundance)
   * [`--motus_save_mgc_read_counts`](#--motus-save-mgc-read-counts)
   * [`--motus_remove_ncbi_ids`](#--motus-remove-ncbi-ids)
* [Postprocessing and visualization options](#postprocessing-and-visualization-options)
   * [`--lowread_filter`](#--lowread-filter)
   * [`--min_frequency`](#--min-frequency)
   * [`--min_samples`](#--min-samples)
   * [`--qiime_tax_agglom_min`](#--qiime-tax-agglom-min)
   * [`--qiime_tax_agglom_max`](#--qiime-tax-agglom-max)
   * [`--ancombc_fdr_cutoff`](#--ancombc-fdr-cutoff)
   * [`--skip_heatmap`](#--skip-heatmap)
   * [`--skip_alpha_rarefaction`](#--skip-alpha-rarefaction)
   * [`--top_taxa`](#--top-taxa)
   * [`--group_of_interest`](#--group-of-interest)
   * [`--skip_alphadiversity`](#--skip-alphadiversity)
   * [`--skip_individalpha`](#--skip-individalpha)
   * [`--skip_betadiversity`](#--skip-betadiversity)
   * [`--run_profile_standardisation`](#--run-profile-standardisation)
   * [`--standaridisation_motus_generation`](#--standaridisation-motus-generation)
   * [`--run_krona`](#--run-krona)
   * [`--krona_taxonomy_directory`](#--krona-taxonomy-directory)
   * [`--standardisation_taxpasta_format`](#--standardisation-taxpasta-format)
   * [`--taxpasta_taxonomy_dir`](#--taxpasta-taxonomy-dir)
   * [`--taxpasta_add_name`](#--taxpasta-add-name)
   * [`--taxpasta_add_rank`](#--taxpasta-add-rank)
   * [`--taxpasta_add_lineage`](#--taxpasta-add-lineage)
   * [`--taxpasta_add_idlineage`](#--taxpasta-add-idlineage)
* [Max job request options](#max-jobs-request-options)
   * [`--max_cpus`](#--max-cpus)
   * [`--max_memory`](#--max-memory)
   * [`--max_time`](#--max-time)
* [Generic options](#generic-options)
   * [`--help`](#--help)
   * [`--version`](#--version)
   * [`--publish_dir_mode`](#--publish-dir-mode)
   * [`--email_on_fail`](#--email-on-fail)
   * [`--plaintext_email`](#--plaintext-email)
   * [`--max_multiqc_email_size`](#--max-multiqc-email-size)
   * [`--monochrome_logs`](#--monochrome-logs)
   * [`--hook_url`](#--hook-url)
   * [`--multiqc_config`](#--multiqc-config)
   * [`--multiqc_logo`](#--multiqc-logo)
   * [`--multiqc_methods_description`](#--multiqc-methods-description)
   * [`--tracedir`](#--tracedir)
   * [`--validate_params`](#--validate-params)
   * [`--show_hidden_params`](#--show-hidden-params)
   * [`--ignore_failed_samples`](#--ignore-failed-samples)
   * [`--report_name`](#--report-name)
   * [`--genome`](#--genome)
   * [`--igenomes_base`](#--igenomes-base)
   * [`--igenomes_ignore`](#--igenomes-ignore)
* [AMR Options](#amr-options)
   * [`--run_amr`](#--run-amr)
   * [`--resistome_threshold`](#--resistome-threshold)
   * [`--amr_index_files`](#--amr-index-files)


## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the full MGscan Metagenomics pipeline on awsbatch is as follows:

### Using AWS Batch

```bash
nextflow run Zymo-Research/aladdin-shotgun \
    -profile awsbatch \
    --design "<path to design CSV file>" \
    --database sourmash-zymo-2024 \
    --run_amr true \
    -work-dir "<work dir on S3>" \
    --awsregion "<AWS Batch region> \
    --awsqueue "<SQS ARN>" \
    --outdir "<output dir on S3>" \
    -r "0.0.14" \
    -name "<report title>"
```

### Using Docker locally

```bash
nextflow run Zymo-Research/aladdin-shotgun \
    -profile docker,dev \
    --design "<path to design CSV file>" \
    --database sourmash-zymo-2024 \
    --run_amr true
    -name "<report title>"
```

### Using SLURM on ZymoCloud

```bash
nextflow run Zymo-Research/aladdin-shotgun \
    -profile slurm \
    --design "<path to design CSV file>" \
    --database sourmash-zymo-2024 \
    --run_amr true \
    -work-dir "<work dir on ZymoCloud>" \
    --partition "<partition name on ZymoCloud>" \
    --outdir "<output dir on ZymoCloud>" \
    -r "0.0.14" \
    -name "<report title>"
```

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull Zymo-Research/aladdin-shotgun
```

### Reproducibility
It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [Zymo-Research/aladdin-shotgun releases page](https://github.com/Zymo-Research/aladdin-shotgun/releases) and find the latest version number - numeric only (eg. `0.0.14`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 0.0.14`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Main arguments

### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. While there are multiple profiles listed below, this pipeline has only been tested with `docker`, `slurm`, and `awsbatch`.

* `awsbatch`
  * A generic configuration profile to be used with AWS Batch.
* `conda`
  * A generic configuration profile to be used with [conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub
* `slurm`
  * A generic configuration profile to be used with slurm clusters. When the slurm profile is used, `apptainer`, which pulls software from Dockerhub, is enabled.
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

When using `awsbatch` profile, one must supply [other options related to AWS batch](#aws-batch-specific-parameters), and supply the locations of [work directory](#-work-dir) and [output directory](#--outdir) on AWS S3.

### `--design`
You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 4 columns, and a header row as shown in the examples in the subsections below.

```bash
--design 'path/to/data/sample_sheet.csv'
```

#### Full analysis
The parameter `--design` is required. It must be a CSV file with the following format.
```
sample,read_1,read_2,group,run_accession
sample1,s1_run1_R1.fastq.gz,s1_run1_R2.fastq.gz,groupA,run1
sample1,s1_run2_R1.fastq.gz,s1_run2_R2.fastq.gz,groupA,run2
sample2,s2_run1_R1.fastq.gz,,groupB,,
sample3,s3_run1_R1.fastq.gz,s3_run1_R2.fastq.gz,groupB,,
```
   - The header line must be present. 
   - The columns "sample", "read_1", "read_2", "group" must be present. Column "run_accession" is optional.
   - The column "sample" contains the name/label for each sample. It can be duplicate. When duplicated, it means the same sample has multiple sequencing runs. In those cases, a different value for "run_accession" is expected. See "sample1" in above example. Sample names must contain only alphanumerical characters or underscores, and must start with a letter.
   - The columns "read_1", "read_2" refers to the paths, including S3 paths, of Read 1 and 2 of Illumina paired-end data. They must be ".fastq.gz" or ".fq.gz" files. When your data are single-end Illumina or PacBio data, simply use "read_1" column, and leave "read_2" column empty. FASTA files from Nanopore data are currently not supported.
   - The column "group" contains the group name/label for comparison purposes in the diversity analysis. If you don't have/need this information, simply leave the column empty, but this column must be present regardless. Same rules for legal characters of sample names apply here too. 
   - The column "run_accesssion" is optional. It is only required when there are duplicates in the "sample" column. This is to mark different run names for the sample. 

## Reference databases
This pipeline includes software tools and pre-prepared reference databases used to assign taxonomy to your metagenomic reads. Supported databases are specified with the argument `--database`, included below.

### `--database`
Select the tool you would like to use for taxonomic profiling and its corresponding reference database. Your choice will dictate the ranking and nomenclature of the taxonomy profile. We provide several popular reference databases to choose from, but recommend the sourmash-compatible database from Zymo Research.

Sourmash prepared GenBank and GTDB databases are available on the sourmash [site](https://sourmash.readthedocs.io/en/latest/databases.html).

Instructions for downloading the MetaPhlAn databases are available on the MetaPhlAn [site](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4#pre-requisites).

* Zymo Jun 2024 Custom Database (`--database sourmash-zymo-2024`)
   * GenBank March 2022 Fungi
   * GenBank March 2022 Protozoa
   * GenBank March 2022 Viral
   * GTDB rs214 DNA database
   * Zymo Custom Host Genomes
   * Zymo Custom Genomes

* Zymo Jul 2023 Custom Database (`--database sourmash-zymo-2023`)
   * GenBank March 2022 Fungi
   * GenBank March 2022 Protozoa
   * GenBank March 2022 Viral
   * GTDB rs207 DNA database
   * Zymo Custom Host Genomes
   * Zymo Custom Genomes

* Metaphlan4.0 Database (`--metaphlan4-db`)
   * Metaphlan mpa_vOct22_CHOCOPhlAnSGB_202212 database
  
This pipeline has currently been tested with [sourmash](https://sourmash.readthedocs.io/en/latest/) and [MetaPhlAn4](https://huttenhower.sph.harvard.edu/metaphlan/).

## Trimming and filtering

### `--preprocessing_qc_tool` 
Specify the tool used for quality control of raw sequencing reads to be FastQC or falco. Falco is designed as a drop-in replacement for FastQC but written in C++ for faster computation. We particularly recommend using falco when using long reads (due to reduced memory constraints), however it is also applicable for short reads.

### `--save_preprocessed_reads`
This saves the FASTQ output from the following tools:

- fastp
- AdapterRemoval
- Porechop
- Filtlong

Depending on the parameters you set, these reads will be a mixture of: adapter clipped, quality trimmed, pair-merged, and length filtered.

### `--shortread_qc_tool`
Specify which tool to use for short-read quality control. This will remove adapters, trim low quality bases, remove reads that are too short, etc. Choose 'DO_NOT_RUN' if you don't want this step performed. Default: fastp

### `--shortread_qc_skipadaptertrim`
Skip the removal of sequencing adapters. \n\nThis often can be useful to speed up run-time of the pipeline when analysing data downloaded from public databases such as the ENA or SRA, as adapters should already be removed (however we recommend to check FastQC results to ensure this is the case).

### `--shortread_qc_adapter1`
Specify a custom forward or R1 adapter sequence to be removed from reads.

If not set, the selected short-read QC tool's defaults will be used.

> Modifies tool parameter(s):
> - fastp: `--adapter_sequence`. fastp default: `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`
> - AdapterRemoval: `--adapter1`. AdapteRemoval2 default: `AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG`

### `--shortread_qc_adapter2`
Specify a custom reverse or R2 adapter sequence to be removed from reads.

If not set, the selected short-read QC tool's defaults will be used.

> Modifies tool parameter(s):
> - fastp: `--adapter_sequence`. fastp default: `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`
> - AdapterRemoval: `--adapter1`. AdapteRemoval2 default: `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT`

### `--shortread_qc_adapterlist`
Allows to supply a file with a list of adapter (combinations) to remove from all files. 

Overrides the --shortread_qc_adapter1/--shortread_qc_adapter2 parameters .

For AdapterRemoval this consists of a two column table with a `.txt` extension: first column represents forward strand, second column for reverse strand. You must supply all possible combinations, one per line, and this list is applied to all files. See AdapterRemoval documentation for more information.

For fastp this consists of a standard FASTA format with a `.fasta`/`.fa`/`.fna`/`.fas` extension. The adapter sequence in this file should be at least 6bp long, otherwise it will be skipped. fastp trims the adapters present in the FASTA file one by one.

> Modifies AdapterRemoval parameter: --adapter-list
> Modifies fastp parameter: --adapter_fasta

### `--shortread_qc_mergepairs`
Turn on the merging of read-pairs of paired-end short read sequencing data. 
> Modifies tool parameter(s):
> - AdapterRemoval: `--collapse`
> - fastp: `-m --merged_out`

### `--shortread_qc_includeunmerged`
Turns on the inclusion of unmerged reads in resulting FASTQ file from merging paired-end sequencing data when using `fastp` and/or `AdapterRemoval`. For `fastp` this means the unmerged read pairs are directly included in the output FASTQ file. For `AdapterRemoval`, additional output files containing unmerged reads are all concatenated into one file by the workflow.

Excluding unmerged reads can be useful in cases where you prefer to have very short reads (e.g. aDNA), thus excluding longer-reads or possibly faulty reads where one of the pair was discarded.

> Adds `fastp` option: `--include_unmerged`

### `--shortread_qc_minlength`
Specifying a mimum read length filtering can speed up profiling by reducing the number of short unspecific reads that need to be match/aligned to the database.

> Modifies tool parameter(s):
> - removed from reads `--length_required`
> - AdapterRemoval: `--minlength`

### `--shortread_complexityfilter_tool`
Specify which tool to use for complexity filtering. This will remove low complexity or highly repetitive sequences that are often not informative. Choose 'DO_NOT_RUN' if you don't want this step performed. Default: bbduk

### `--shortread_complexityfilter_entropy`
Specify the minimum 'entropy' value for complexity filtering for BBDuk or PRINSEQ++.

Note that this value will only be used for PRINSEQ++ if `--shortread_complexityfilter_prinseqplusplus_mode` is set to `entropy`.

Entropy here corresponds to the amount of sequence variation exists within the read. Higher values correspond to more variety, and thus will likely reslut in more specific matching to a taxon's reference genome. The trade off here is fewer reads (or abundance information) available for having a confident identification.


> Modifies tool parameter(s):
> - BBDuk: `entropy=`
> - PRINSEQ++:  `-lc_entropy`
                    
### `--shortead_complexityfilter_bbduk_windowsize`
Specify the window size to calculate the level entropy within for BBDuk.

> Modifies tool parameter(s):
> - BBDuk: `entropywindow=`

### `--shortread_complexityfilter_bbduk_mask`
Turn on masking of low-complexity reads (i.e., replacement with `N`) rather than removal.

> Modifies tool parameter(s)
> - BBDuk: `entropymask=`

### `--shortread_complexityfilter_fastp_threshold`
Specify the minimum sequence complexity value for fastp. This value corresponds to the percentage of bases that is different from it's adjacent bases.

> Modifies tool parameter(s):
> - removed from reads `--complexity_threshold`

### `--shortread_complexityfilter_prinseqplusplus_mode`
Specify the complexity filter mode for PRINSEQ++

### `--shortread_complexityfilter_prinseqplusplus_dustscore`
Specify the minimum dust score below which low-complexity reads will be removed. A DUST score is based on how often different tri-nucleotides occur along a read.

> Modifies tool parameter(s):
> - PRINSEQ++: `--lc_dust`

### `--save_complexityfiltered_reads`
Specify whether to save the final complexity filtered reads in your results directory (`--outdir`).

## Preprocessing Long Read QC Options

### `--perform_longread_qc`
Turns on long read quality control steps (adapter clipping, length and/or quality filtering.)

Removing adapters (if present) is recommend to reduce false-postive hits that may occur from 'dirty' or 'contaminated' reference genomes in a profiling database that contain accidentially incorporated adapter sequences.

Length filtering, and quality filtering can speed up alignment by reducing the number of unspecific reads that need to be aligned.


### `--longread_qc_skipadaptertrim`
Skip removal of adapters by Porechop. This can be useful in some cases to speed up run time - particularly when you are running data downloading from public databases such as the ENA/SRA that should already have adapters removed. We recommend that you check your FastQC results this is indeed the case.

### `--longread_qc_skipqualityfilter`
Skip removal of quality filtering with Filtlong. This will skip length, percent reads, and target bases filtering (see other `--longread_qc_qualityfilter_*` parameters).

### `--longread_qc_qualityfilter_minlength`
Specify the minimum of length of reads to be kept for downstream analysis.

> Modifies tool parameter(s):
> - Filtlong: `--min_length`

### `--longread_qc_qualityfilter_keeppercent`
Throw out the remaining percentage of reads outside the value. This is measured by bp, not by read count. So this option throws out the worst e.g. 10% of read bases if the parameter is set to `90`.  _Modified from [Filtlong documentation](https://github.com/rrwick/Filtlong)_

> Modifies tool parameter(s):
> - Filtlong: `--keep_percent`

### `--longread_qc_qualityfilter_targetbases`
Removes the worst reads until only the specified value of bases remain, useful for very large read sets. If the input read set is less than the specified value, this setting will have no effect. _Modified from [Filtlong documentation](https://github.com/rrwick/Filtlong)_
> Modifies tool parameter(s):
> - Filtlong: `--keep_percent`

## Preprocessing Host Removal Options

### `--perform_shortread_hostremoval`
Turns on the ability to remove short-reads from the that derived from a known organism, using Bowtie2 and samtools.

This subworkflow is useful to remove reads that may come from a host, or a known contamination like the human reference genome. Human DNA contamination of (microbial) reference genomes is well known, so removal of these prior profiling both reduces the risks of false positives, and in _some cases_ a faster runtime (as less reads need to be profiled).

Alternatively, you can include the reference genome within your profiling databases and can turn off this subworkflow, with the trade off of a larger taxonomic profiling database.

### `--perform_longread_hostremoval`
Turns on the ability to remove long-reads from the that derived from a known organism, using minimap2 and samtools

This subworkflow is useful to remove reads that may come from a host, or a known contamination like the human reference genome. Human DNA contamination of (microbial) reference genomes is well known, so removal of these prior profiling both reduces the risks of false positives, and in _some cases_ a faster runtime (as less reads need to be profiled).

Alternatively, you can include the reference genome within your profiling databases and can turn off this subworkflow, with the trade off of a larger taxonomic profiling database.

### `--hostremoval_reference`
Specify a path to the FASTA file (optionally gzipped) of the reference genome of the organism to be removed.\n\nIf you have two or more host organisms or contaminants you wish to remove, you can concatenate the FASTAs of the different taxa into a single one to provide to the pipeline.


### `--shortread_hostremoval_index`
Specify the path to a _directory_ containing pre-made Bowtie2 reference index files (i.e. the directory containing `.bt1`, `.bt2` files etc.). These should sit in the same directory alongside the the reference file specified in `--hostremoval_reference`.

Specifying premade indices can speed up runtime of the host-removal step, however if not supplied the pipeline will generate the indices for you.

### `--longread_hostremoval_index`
Specify path to a pre-made Minimap2 index file (.mmi) of the host removal reference file given to `--hostremoval_reference`.

Specifying a premade index file can speed up runtime of the host-removal step, however if not supplied the pipeline will generate the indices for you.


### `--save_hostremoval_index`
Save the output files of the in-built indexing of the host genome.

This is recommend to be turned on if you plan to use the same reference genome multiple times, as supplying the directory or file to `--shortread_hostremoval_index` or `--longread_hostremoval_index` respectively can speed up runtime of future runs. Once generated, we recommend you place this file _outside_ of your run results directory in a central 'cache' directory you and others using your machine can access and supply to the pipeline.

### `--save_hostremoval_bam`
Save the reads mapped to the reference genome and off-target reads in BAM format as output by the respective hostremoval alignment tool.

This can be useful if you wish to perform other analyses on the host organism (such as host-microbe interaction), however, you should consider whether the default mapping parameters of Bowtie2 (short-read) or minimap2 (long-read) are optimised to your context.

### `--save_hostremoval_unmapped`
Save only the reads NOT mapped to the reference genome in FASTQ format (as exported from `samtools view` and `bam2fq`).

This can be useful if you wish to perform other analyses on the off-target reads from the host mapping, such as manual profiling or _de novo_ assembly.

## Preprocessing Run Merging Options

### `--perform_runmerging`
Turns on the concatenation of sequencing runs or libraries with the same sample name.

This can be useful to ensure you get a single profile per sample, rather than one profile per run or library. Note that in some cases comparing profiles of independent _libraries_ may be useful, so this parameter may not always be suitable.

### `--save_runmerged_reads`
Save the run- and library-concatenated reads of a given sample in FASTQ format.

> Only samples that went through the run-merging step of the pipeline will be stored in the resulting directory.

> If you wish to save the files that go to the classification/profiling steps for samples that _did not_ go through run merging, you must supply the appropriate upstream `--save_<preprocessing_step>` flag.

## Profiling Options

### `--centrifuge_save_reads`
Save mapped (SAM, FASTQ) and unmapped (FASTQ) reads from alignment step of centrifuge in your output results directory.

> Modifies tool parameter(s):
> - centrifuge: `--un-gz`, `--al-gz`, `--un-conc-gz`, `--al-conc-gz`, `--out-fmt`

### `--diamond_output_format`
DIAMOND can produce output in a number of different formats, you can specify here which to produce.

Note that DIAMOND can only produce one format at a time, and depending on which you pick, some downstream steps may not be executed. For example, selecting `daa` or `sam` will mean you will not get a tabular taxonomic profile as with the other tools.

Will be overriden by `--diamond_save_reads.`

> Modifies tool parameter(s):
> - diamond blastx: `--outfmt`
                    
### `--diamond_save_reads`
Save aligned reads in SAM format from alignment step of DIAMOND in your output results directory.

Note this explicitly overrides `--diamond_output_format` to produce the SAM file, and no taxon table will be generated.

> Modifies tool parameter(s):
> - DIAMOND: `--outfmt`
                    
### `--kaiju_taxon_rank`
Specify the taxonomic level(s) to be displayed in the resulting Kaiju taxon table, as generated by the kaiju2table helper tool.

This can be either a single level (e.g. `species`), or a comma separated list to display the full taxonomic path (e.g. `superkingdom,phylum,class,order,family,genus,species.`).
> Modifies tool parameter(s):
> - kaiju2table: `-l`
                    
### `--kraken2_save_reads`
Save reads that do and do not have a taxonomic classification in your output results directory in FASTQ format.
> Modifies tool parameter(s):
> - kraken2: `--classified-out` and `--unclassified-out`

### `--kraken2_save_readclassification`
Save a text file that contains a list of each read that had a taxonomic assignment, with information on specific taxonomic taxonomic assignment that that read recieved.

> Modifies tool parameter(s):
> - kraken2: `--output`

### `--kraken2_save_minimizers`
Turn on saving minimizer information in the kraken2 report thus increasing to an eight column layout.

Adds `--report-minimizer-data` to the kraken2 command.

### `--krakenuniq_save_reads`
Save reads that do and do not have a taxonomic classification in your output results directory in FASTQ format.

> Modifies tool parameter(s):
> - krakenuniq: `--classified-out` and `--unclassified-out`

### `--krakenuniq_ram_chunk_size`
nf-core/taxprofiler utilises a 'low memory' option for KrakenUniq that can reduce the amount of RAM the process requires using the `--preloaded` option.

A further extension to this option is that you can specify how large each chunk of the database should be that gets loaded into memory at any one time. You can specify the amount of RAM to chunk the database to with this parameter, and is particularly useful for people with limited computational resources.

More information about this parameter can be seen [here](https://github.com/fbreitwieser/krakenuniq/blob/master/README.md#new-release-v07).

> Modifies KrakenUniq parameter: --preload-size

### `--krakenuniq_save_readclassification`
Save a text file that contains a list of each read that had a taxonomic assignment, with information on specific taxonomic taxonomic assignment that that read recieved.

> Modifies tool parameter(s):
> - krakenuniq: `--output`

### `--malt_mode`
Specify which version of MALT alignment to use.

BlastN is generally recommended (nucleotide-nucleotide alignment), but particularly for very short reads (such as aDNA), whereas BlastX mode is similar to DIAMOND and will translate the nucleotide to amino acid sequences. Note each type of alignment mode requires different parameters during database construction. Refer to the MALT manual for more information.

> Modifies tool parameter(s):
> - malt-run: `-mode`

### `--malt_save_reads`
Turns on saving of MALT aligned reads in SAM format.

Note that the SAM format produce by MALT is not completely valid, and may not work with downstream tools.

> Modifies tool parameter(s):
> - malt-run: `--alignments`, `-za`

### `--malt_generate_megansummary`
Turns on saving of MALT output in an additional MEGAN summary file (`.megan`) that can be loaded into the MEGAN metagenomic exploration tool.

Note: this file is generated not directly from MALT but rather then MEGAN utility script `rma2info`.
> Modifies tool parameter(s):
> - rma2info: `-es`

### `--sourmash_kmersize`
Sourmash breaks reads down into k-mers which are compared against database for taxonomic assignment. Sourmash recommends k=31 for most cases. k=51 increase specificity at the expense of sensitivity, therefore should be used when most species in your samples have very closely related genomes in the database. k=21 is useful when species in your sample have no good reference genomes and you just want classification on the family level or above.

### `--sourmash_threshold_bp`
Sourmash will only report a match to a genome in the database when the estimated overlap between sample and the reference is at least this many base pairs. Choose a smaller number to increase sensitivity to low abundant species or small genomes like virus. Default is 5,000bp.

### `--sourmash_trim_low_abund`
Runs trim-low-abund.py from the khmer software package, which removes low abundance, potentially erroneous kmers from data before taxonomy assignment. Only applies when sourmash databases were selected.

### `--motus_use_relative_abundance`
Turn on printing relative abundance instead of counts.

### `--motus_save_mgc_read_counts`
Turn on saving the mgc reads count.

### `--motus_remove_ncbi_ids`
Turn on removing NCBI taxonomic IDs.

## Postprocessing and visualization options

### `--lowread_filter`
Applies a read count filter for samples before downstream analysis. Samples must have at least this many reads with assigned taxonomy to be included in analyses such as alpha and beta diversity. Default is 50,000.

### `--min_frequency`
Removes taxa with fewer than this many total reads assigned among all samples. Default setting of 1 means this filter is disabled, and all taxa will be kept.


### `--min_samples`
Remove taxa that were detected in fewer than this many samples. Default setting of 1 means this filter is disabled, and all taxa will be kept.

### `--qiime_tax_agglom_min`
There are seven levels of taxonomy rankings: kingdom[1], phylum[2], class[3], order[4], family[5], genus[6], and species[7]. ANCOMBC and heatmap analysis will be conducted at the levels between this minimum parameter and the qiime maximum taxa parameter. Default is 7, species level.

### `--qiime_tax_agglom_max`
There are seven levels of taxonomy rankings: kingdom[1], phylum[2], class[3], order[4], family[5], genus[6], and species[7]. ANCOMBC and heatmap analysis will be conducted at the levels between this maximum parameter and the qiime minimum taxa parameter. Default is 7, species level.

### `--ancombc_fdr_cutoff`
FDR threshold for a taxa to be considered significant and displayed on the ANCOM-BC plot

### `--skip_heatmap`
Remove off heatmap plot from MultiQC report

### `--skip_alpha_rarefaction`
Turn off alpha rarefaction analysis and remove Alpha Rarefaction plot from MultiQC report

### `--top_taxa`
Specify the maximum number of top taxa from each sample you want to compare in QIIME_HEATMAP

### `--group_of_interest`
Group(s) of taxa that are of specific interest. Abundances of these taxa will be summarized in a separate file and plotted in a separate section in the report to highlight these taxa.

### `--skip_alphadiversity`
Turn off alpha diversity analysis and remove alpha diversity plot from MultiQC report

### `--skip_individalpha`
Turn off the chart showing individual alpha diversity per sample.

### `--skip_betadiversity`
Turn off beta diversity analysis and remove beta diversity plot from MultiQC reports

### `--run_profile_standardisation`
Turns on standardisation of output OTU tables across all tools; each into a TSV format following the following scheme:

|TAXON    | SAMPLE_A | SAMPLE_B |
|---------|----------|----------|
| taxon_a | 32       | 123      |
| taxon_b | 1        | 5        |

This currently only is generated for mOTUs.

### `--standaridisation_motus_generation`
Turn on the saving of the taxonomic output in BIOM format (`.biom`) in the results directory of your pipeline run, instead of the default TSV format.

Note this file is from the output of the `motus merge` command.

> Modifies tool parameter(s):
> - `-B -o`

### `--run_krona`
Turn on the generation of Krona interactive pie-chart HTMLs for a selection of profilers.

The tools currently supported are:
- centrifuge
- kraken2
- kaiju
- MALT
  
### `--krona_taxonomy_directory`
Specify a path to a Krona taxonomy database directory (i.e. a directory containing a krona generated `.tab` file).

This is only required for generating Krona plots of MALT output.

Note this taxonomy database must be downloaded and generated with the `updateTaxonomy.sh` script from the krona-tools package.

### `--standardisation_taxpasta_format`
The desired output format.

### `--taxpasta_taxonomy_dir` 
This arguments provides the path to the directory containing taxdump files. At least nodes.dmp and names.dmp are required. A merged.dmp file is optional.

Modifies tool parameter(s):
-taxpasta: `--taxpasta_taxonomy_dir`

### `--taxpasta_add_name`
The standard output format of taxpasta is a two-column table including the read counts and the integer taxonomic ID. The taxon name can be added as additional information to the output table.

Modifies tool parameter(s):
- taxpasta: `--taxpasta_add_name`

### `--taxpasta_add_rank`
The standard output format of taxpasta is a two-column table including the read counts and the integer taxonomic ID. The taxon rank can be added as additional information to the output table.

Modifies tool parameter(s):
- taxpasta: `--taxpasta_add_rank`

### `--taxpasta_add_lineage`
The standard output format of taxpasta is a two-column table including the read counts and the integer taxonomic ID. The taxon's entire lineage with the taxon names separated by semi-colons can be added as additional information to the output table.

Modifies tool parameter(s):
- taxpasta: `--taxpasta_add_lineage`

### `--taxpasta_add_idlineage`
The standard output format of taxpasta is a two-column table including the read counts and the integer taxonomic ID. The taxon's entire lineage with the taxon identifiers separated by semi-colons can be added as additional information to the output table.

Modifies tool parameter(s):
- taxpasta: `--taxpasta_add_idlineage`

## Max job request options

### `--max_cpus`
Maximum number of CPUs that can be requested for any single job. Use to set an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`

### `--max_memory`
Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`

### `--max_time`
Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`

## Generic options

### `--help`
Display help text.

### `--version`
Display version and exit.

### `--publish_dir_mode`
The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.

### `--email_on_fail`
An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.

### `--plaintext_email`
Send plain-text email instead of HTML.

### `--max_multiqc_email_size`
File size limit when attaching MultiQC reports to summary emails.

### `--monochrome_logs`
Do not use coloured log outputs.

### `--hook_url`
Incoming hook URL for messaging service. Currently, MS Teams and Slack are supported.

### `--multiqc_config`
Custom config file to supply to MultiQC.

### `--multiqc_logo`
Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file

### `--multiqc_methods_description`
Custom MultiQC yaml file containing HTML including a methods description.

### `--tracedir`
Directory to keep pipeline Nextflow logs and reports.

### `--validate_params`
Boolean whether to validate parameters against the schema at runtime

### `--show_hidden_params`
By default, parameters set as _hidden_ in the schema are not shown on the command line when a user runs with `--help`. Specifying this option will tell the pipeline to show all parameters

### `--ignore_failed_samples`
Whether to ignore samples that fail QC and taxonomy profiling steps and carry on with good samples. Default True. Set to False to force the pipeline to abort when any samples fails.

### `--report_name`
Name for report title and report file name.

### `--genome`
If using a reference genome configured in the pipeline using iGenomes, use this parameter to give the ID for the reference. This is then used to build the full paths for all required reference genome files e.g. `--genome GRCh38`. \n\nSee the [nf-core website docs](https://nf-co.re/usage/reference_genomes) for more details.

### `--igenomes_base`
Directory / URL base for iGenomes references.

### `--igenomes_ignore`
Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.

## AMR Options

### `--run_amr`
Whether or not to conduct analysis for antimicrobial resistance genes.

### `--resistome_threshold`
The minimum proportion of nucleotides that have aligned reads to for any AMR gene to be reported.

### `--amr_index_files`
Path to files related to AMR analysis.

















####################################################################### TAXPROFILER USAGE.MD #############################################################################################################
# nf-core/taxprofiler: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/taxprofiler/usage](https://nf-co.re/taxprofiler/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

nf-core/taxprofiler is a pipeline for highly-parallelised taxonomic classification and profiling of shotgun metagenomic data across multiple tools simultaneously. In addition to multiple classification and profiling tools, at the same time it allows you to performing taxonomic classification and profiling across multiple databases and settings per tool, as well as produces standardised output tables to allow immediate cross comparison of results between tools.

To run nf-core/taxprofiler, at a minimum two you require two inputs:

- a sequencing read samplesheet
- a database samplesheet

Both contain metadata and paths to the data of your input samples and databases.

When running nf-core/taxprofiler, every step and tool is 'opt in'. To run a given classifier or profiler you must make sure to supply both a database in your `<database>.csv` and supply `--run_<profiler>` flag to your command. Omitting either will result in the profiling tool not executing.

nf-core/taxprofiler also includes optional pre-processing (adapter clipping, merge running etc.) or post-processing (visualisation) steps. These are also opt in with a `--perform_<step>` flag. In some cases, the pre- and post-processing steps may also require additional files. Please check the parameters tab of this documentation for more information.

Please see the rest of this page for information about how to prepare input samplesheets and databases and how to run Nextflow pipelines. See the [parameters](https://nf-co.re/taxprofiler/parameters) documentation for more information about specific options the pipeline also offers.

## Samplesheet inputs

nf-core/taxprofiler can accept as input raw or preprocessed single- or paired-end short-read (e.g. Illumina) FASTQ files, long-read FASTQ files (e.g. Oxford Nanopore), or FASTA sequences (available for a subset of profilers).

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 6 columns, and a header row as shown in the examples below. Furthermother, nf-core/taxprofiler also requires a second comma-separated file of 3 columns with a header row as in the examples below.

This samplesheet is then specified on the command line as follows:

```console
--input '[path to samplesheet file]' --databases '[path to database sheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate different runs FASTQ files of the same sample before performing profiling, when `--perform_runmerging` is supplied. Below is an example for the same sample sequenced across 3 lanes:

```console
sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta
2612,run1,ILLUMINA,2612_run1_R1.fq.gz,,
2612,run2,ILLUMINA,2612_run2_R1.fq.gz,,
2612,run3,ILLUMINA,2612_run3_R1.fq.gz,2612_run3_R2.fq.gz,

```

> ⚠️ Runs of the same sample sequenced on Illumina platforms with a combination of single and paired-end data will **not** be run-wise concatenated, unless pair-merging is specified. In the example above, `run3` will be profiled independently of `run1` and `run2` if pairs are not merged.

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 6 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data, as well as long-read FASTA files may look something like the one below. This is for 6 samples, where `2612` has been sequenced twice.

```console
sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta
2611,ERR5766174,ILLUMINA,,,/<path>/<to>/fasta/ERX5474930_ERR5766174_1.fa.gz
2612,ERR5766176,ILLUMINA,/<path>/<to>/fastq/ERX5474932_ERR5766176_1.fastq.gz,/<path>/<to>/fastq/ERX5474932_ERR5766176_2.fastq.gz,
2612,ERR5766180,ILLUMINA,/<path>/<to>/fastq/ERX5474936_ERR5766180_1.fastq.gz,,
2613,ERR5766181,ILLUMINA,/<path>/<to>/fastq/ERX5474937_ERR5766181_1.fastq.gz,/<path>/<to>/fastq/ERX5474937_ERR5766181_2.fastq.gz,
ERR3201952,ERR3201952,OXFORD_NANOPORE,/<path>/<to>/fastq/ERR3201952.fastq.gz,,
```

> ⚠️ Input FASTQ and FASTA files _must_ be gzipped

> ⚠️ While one can include both short-read and long-read data in one run, we recommend that you split these across _two_ pipeline runs and database sheets (see below). This will allow classification optimisation for each data type, and make MultiQC run-reports more readable (due to run statistics having vary large number differences).

| Column                | Description                                                                                                                                                                                              |
| --------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`              | Unique sample name [required].                                                                                                                                                                           |
| `run_accession`       | Run ID or name unique for each (pairs of) file(s) .Can also supply sample name again here, if only a single run was generated [required].                                                                |
| `instrument_platform` | Sequencing platform reads generated on, selected from the EBI ENA [controlled vocabulary](https://www.ebi.ac.uk/ena/portal/api/controlledVocab?field=instrument_platform) [required].                    |
| `fastq_1`             | Path or URL to sequencing reads or for Illumina R1 sequencing reads in FASTQ format. GZipped compressed files accepted. Can be left empty if data in FASTA is specifed. Cannot be combined with `fasta`. |
| `fastq_2`             | Path or URL to Illumina R2 sequencing reads in FASTQ format. GZipped compressed files accepted. Can be left empty if single end data. Cannot be combined with `fasta`.                                   |
| `fasta`               | Path or URL to long-reads or contigs in FASTA format. GZipped compressed files accepted. Can be left empty if data in FASTA is specifed. Cannot be combined with `fastq_1` or `fastq_2`.                 |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Full database sheet

nf-core/taxprofiler supports multiple databases being classified/profiled against in parallel for each tool.

Databases can be supplied either in the form of a compressed `.tar.gz` archive of a directory containing all relevant database files or the path to a directory on the filesystem.

> ⚠️ nf-core/taxprofiler does not provide any databases by default, nor does it currently generate them for you. This must be performed manually by the user. See below for more information of the expected database files.

The pipeline takes the paths and specific classification/profiling parameters of the tool of these databases as input via a four column comma-separated sheet.

> ⚠️ To allow user freedom, nf-core/taxprofiler does not check for mandatory or the validity of non-file database parameters for correct execution of the tool - excluding options offered via pipeline level parameters! Please validate your database parameters (cross-referencing [parameters](https://nf-co.re/taxprofiler/parameters, and the given tool documentation) before submitting the database sheet! For example, if you don't use the default read length - Bracken will require `-r <read_length>` in the `db_params` column.

An example database sheet can look as follows, where 7 tools are being used, and `malt` and `kraken2` will be used against two databases each.

`kraken2` will be run twice even though only having a single 'dedicated' database because specifying `bracken` implies first running `kraken2` on the `bracken` database, as required by `bracken`.

```console
tool,db_name,db_params,db_path
malt,malt85,-id 85,/<path>/<to>/malt/testdb-malt/
malt,malt95,-id 90,/<path>/<to>/malt/testdb-malt.tar.gz
bracken,db1,;-r 150,/<path>/<to>/bracken/testdb-bracken.tar.gz
kraken2,db2,--quick,/<path>/<to>/kraken2/testdb-kraken2.tar.gz
krakenuniq,db3,,/<path>/<to>/krakenuniq/testdb-krakenuniq.tar.gz
centrifuge,db1,,/<path>/<to>/centrifuge/minigut_cf.tar.gz
metaphlan4,db1,,/<path>/<to>/metaphlan4/metaphlan_database/
motus,db_mOTU,,/<path>/<to>/motus/motus_database/
```

For Bracken, if you wish to supply any parameters to either the Kraken or Bracken step you **must** have a _semi-colon_ `;` list as in `db_params`. This is to allow to specify the Kraken2 parameters before, and Bracken parameters after the `;` as Bracken is a two step process. This is particularly important if you supply a Bracken database with a non-default read length parameter. If you do not have any parameters to specify, you can leave this as empty.

Column specifications are as follows:

| Column      | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `tool`      | Taxonomic profiling tool (supported by nf-core/taxprofiler) that the database has been indexed for [required]. Please note that `bracken` also implies running `kraken2` on the same database.                                                                                                                                                                                                                                                                                                                                                                                                              |
| `db_name`   | A unique name per tool for the particular database [required]. Please note that names need to be unique across both `kraken2` and `bracken` as well, even if re-using the same database.                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `db_params` | Any parameters of the given taxonomic classifier/profiler that you wish to specify that the taxonomic classifier/profiling tool should use when profiling against this specific database. Can be empty to use taxonomic classifier/profiler defaults. Must not be surrounded by quotes [required]. We generally do not recommend specifying parameters here that turn on/off saving of output files or specifying particular file extensions - this should be already addressed via pipeline parameters. For Bracken databases, must at a minimum contain a `;` separating Kraken2 from Bracken parameters. |
| `db_path`   | Path to the database. Can either be a path to a directory containing the database index files or a `.tar.gz` file which contains the compressed database directory with the same name as the tar archive, minus `.tar.gz` [required].                                                                                                                                                                                                                                                                                                                                                                       |

> 💡 You can also specify the same database directory/file twice (ensuring unique `db_name`s) and specify different parameters for each database to compare the effect of different parameters during classification/profiling.

nf-core/taxprofiler will automatically decompress and extract any compressed archives for you.

The (uncompressed) database paths (`db_path`) for each tool are expected to contain the contents of:

- [**Bracken**:](#bracken-custom-database) output of the combined `kraken2-build` and `bracken-build` process.
- [**Centrifuge**:](#centrifuge-custom-database) output of `centrifuge-build`.
- [**DIAMOND**:](#diamond-custom-database) output of `diamond makedb`.
- [**Kaiju**:](#kaiju-custom-database) output of `kaiju-makedb`.
- [**Kraken2**:](#kraken2-custom-database) output of `kraken2-build` command(s).
- [**KrakenUniq**:](#krakenuniq-custom-database) output of `krakenuniq-build` command(s).
- [**MALT**](#malt-custom-database) output of `malt-build`.
- [**MetaPhlAn4**:](#metaphlan4-custom-database) output of with `metaphlan --install` or downloaded from links on the [MetaPhlAn4 wiki](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4.0#customizing-the-database).
- [**mOTUs**:](#motus-custom-database) is composed of code and database together.

Click the links in the list above for short quick-reference tutorials how to generate custom databases for each tool.

## Running the pipeline

The typical command for running the pipeline is as follows:

```console
nextflow run nf-core/taxprofiler --input samplesheet.csv --databases databases.csv --outdir <OUTDIR> -profile docker --run_<TOOL1> --run_<TOOL2>
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

When running nf-core/taxprofiler, every step and tool is 'opt in'. To run a given classifier/profiler you must make sure to supply both a database in your `<database>.csv` and supply `--run_<profiler>` flag to your command. Omitting either will result in the classification/profiling tool not executing. If you wish to perform pre-processing (adapter clipping, merge running etc.) or post-processing (visualisation) steps, these are also opt in with a `--perform_<step>` flag. In some cases, the pre- and post-processing steps may also require additional files. Please check the parameters tab of this documentation for more information.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Sequencing quality control

[`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. nf-core taxprofiler offers [`falco`](https://github.com/smithlabcode/falco) as an drop-in replacement, with supposedly better improvement particularly for long reads.

### Preprocessing Steps

nf-core/taxprofiler offers four main preprocessing steps for preprocessing raw sequencing reads:

- [**Read processing**](#read-processing): adapter clipping and pair-merging.
- [**Complexity filtering**](#complexity-filtering): removal of low-sequence complexity reads.
- [**Host read-removal**](#host-read-removal): removal of reads aligning to reference genome(s) of a host.
- [**Run merging**](#run-merging): concatenation of multiple FASTQ chunks/sequencing runs/libraries of a sample.

#### Read Processing

Raw sequencing read processing in the form of adapter clipping and paired-end read merging can be activated via the `--perform_shortread_qc` or `--perform_longread_qc` flags.

It is highly recommended to run this on raw reads to remove artifacts from sequencing that can cause false positive identification of taxa (e.g. contaminated reference genomes) and/or skews in taxonomic abundance profiles. If you have public data, normally these should have been corrected for, however you should still check that these steps have indeed been already performed.

There are currently two options for short-read preprocessing: [`fastp`](https://github.com/OpenGene/fastp) or [`adapterremoval`](https://github.com/MikkelSchubert/adapterremoval).

For adapter clipping, you can either rely on the tool's default adapter sequences, or supply your own adapters (`--shortread_qc_adapter1` and `--shortread_qc_adapter2`)
By default, paired-end merging is not activated. In this case paired-end 'alignment' against the reference databases is performed where supported, and if not, supported pairs will be independently classified/profiled. If paired-end merging is activated you can also specify whether to include unmerged reads in the reads sent for classification/profiling (`--shortread_qc_mergepairs` and `--shortread_qc_includeunmerged`).
You can also turn off clipping and only perform paired-end merging, if requested. This can be useful when processing data downloaded from the ENA, SRA, or DDBJ (`--shortread_qc_skipadaptertrim`).
Both tools support length filtering of reads and can be tuned with `--shortread_qc_minlength`. Performing length filtering can be useful to remove short (often low sequencing complexity) sequences that result in unspecific classification and therefore slow down runtime during classification/profiling, with minimal gain.

There is currently one option for long-read Oxford Nanopore processing: [`porechop`](https://github.com/rrwick/Porechop).

For both short-read and long-read preprocessing, you can optionally save the resulting processed reads with `--save_preprocessed_reads`.

#### Complexity Filtering

Complexity filtering can be activated via the `--perform_shortread_complexityfilter` flag.

Complexity filtering is primarily a run-time optimisation step. It is not necessary for accurate taxonomic classification/profiling, however it can speed up run-time of each tool by removing reads with low-diversity of nucleotides (e.g. with mono-nucleotide - `AAAAAAAA`, or di-nucleotide repeats `GAGAGAGAGAGAGAG`) that have a low-chance of giving an informative taxonomic ID as they can be associated with many different taxa. Removing these reads therefore saves computational time and resources.

There are currently three options for short-read complexity filtering: [`bbduk`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/), [`prinseq++`](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus), and [`fastp`](https://github.com/OpenGene/fastp#low-complexity-filter).

There is one option for long-read quality filtering: [`Filtlong`](https://github.com/rrwick/Filtlong)

The tools offer different algorithms and parameters for removing low complexity reads and quality filtering. We therefore recommend reviewing the pipeline's [parameter documentation](https://nf-co.re/taxprofiler/parameters) and the documentation of the tools (see links above) to decide on optimal methods and parameters for your dataset.

You can optionally save the FASTQ output of the run merging with the `--save_complexityfiltered_reads`. If running with `fastp`, complexity filtering happens inclusively within the earlier shortread preprocessing step. Therefore there will not be an independent pipeline step for complexity filtering, and no independent FASTQ file (i.e. `--save_complexityfiltered_reads` will be ignored) - your complexity filtered reads will also be in the `fastp/` folder in the same file(s) as the preprocessed read.

> ⚠️ For nanopore data: we do not recommend performing any read preprocessing or complexity filtering if you are using ONTs Guppy toolkit for basecalling and post-processing.

#### Host-Read Removal

Removal of possible-host reads from FASTQ files prior classification/profiling can be activated with `--perform_shortread_hostremoval` or `--perform_longread_hostremoval`.

Similarly to complexity filtering, host-removal can be useful for runtime optimisation and reduction in misclassified reads. It is not always necessary to report classification of reads from a host when you already know the host of the sample, therefore you can gain a run-time and computational advantage by removing these prior typically resource-heavy classification/profiling with more efficient methods. Furthermore, particularly with human samples, you can reduce the number of false positives during classification/profiling that occur due to host-sequence contamination in reference genomes on public databases.

nf-core/taxprofiler currently offers host-removal via alignment against a reference genome with Bowtie2 for short reads and minimap2 for long reads, and the use of the unaligned reads for downstream classification/profiling.

You can supply your reference genome in FASTA format with `--hostremoval_reference`. You can also optionally supply a directory containing pre-indexed Bowtie2 index files with `--shortread_hostremoval_index` or a minimap2 `.mmi` file for `--longread_hostremoval_index`, however nf-core/taxprofiler will generate these for you if necessary. Pre-supplying the index directory or files can greatly speed up the process, and these can be re-used.

> 💡 If you have multiple taxa or sequences you wish to remove (e.g., the host genome and then also PhiX - common quality-control reagent during sequencing) you can simply concatenate the FASTAs of each taxa or sequences into a single reference file.

#### Run Merging

For samples that may have been sequenced over multiple runs, or for FASTQ files split into multiple chunks, you can activate the ability to merge across all runs or chunks with `--perform_runmerging`.

For more information how to set up your input samplesheet, see [Multiple runs of the same sample](#multiple-runs-of-the-same-sample).

Activating this functionality will concatenate the FASTQ files with the same sample name _after_ the optional preprocessing steps and _before_ classification/profiling. Note that libraries with runs of different pairing types will **not** be merged and this will be indicated on output files with a `_se` or `_pe` suffix to the sample name accordingly.

You can optionally save the FASTQ output of the run merging with the `--save_runmerged_reads`.

#### Classification and Profiling

The following sections provide tips and suggestions for running the different taxonomic classification and profiling tools _within the pipeline_. For advice and/or guidance whether you should run a particular tool on your specific data, please see the documentation of each tool!

An important distinction between the different tools in included in the pipeline is classification versus profiling. Taxonomic _classification_ is concerned with simply detecting the presence of species in a given sample. Taxonomic _profiling_ involves additionally estimating the _abundance_ of each species.

Note that not all taxonomic classification tools (e.g. Kraken, MALT, Kaiju) performs _profiling_, but all taxonomic profilers (e.g. MetaPhlAn, mOTUs, Bracken) must perform some form of _classification_ prior to profiling.

For advice as to which tool to run in your context, please see the documentation of each tool.

> 🖊️ If you would like to change this behaviour, please contact us on the [nf-core slack](https://nf-co.re/join) and we can discuss this.

Not all tools currently have dedicated tips, suggestions and/or recommendations, however we welcome further contributions for existing and additional tools via pull requests to the [nf-core/taxprofiler repository](https://github.com/nf-core/taxprofiler)!

##### Bracken

You must make sure to also activate Kraken2 to run Bracken in the pipeline.

It is unclear whether Bracken is suitable for running long reads, as it makes certain assumptions about read lengths. Furthemore, during testing we found issues where Bracken would fail on the long-read test data.

Therefore currently nf-core/taxprofiler does not run Bracken on data specified as being sequenced with `OXFORD_NANOPORE` in the input samplesheet.

##### Centrifuge

Centrifuge currently does not accept FASTA files as input, therefore no output will be produced for these input files.

##### DIAMOND

DIAMOND only allows output of a single file format at a time, therefore parameters such `--diamond_save_reads` supplied will result in only aligned reads in SAM format will be produced, no taxonomic profiles will be available. Be aware of this when setting up your pipeline runs, depending on your particular use case.

##### Kaiju

Currently, no specific tips or suggestions.

##### Kraken2

Currently, no specific tips or suggestions.

##### KrakenUniq

Currently, no specific tips or suggestions.

##### MALT

MALT does not support paired-end reads alignment (unlike other tools), therefore nf-core/taxprofiler aligns these as indepenent files if read-merging is skipped. If you skip merging, you can sum or average the results of the counts of the pairs.

Krona can only be run on MALT output if path to Krona taxonomy database supplied to `--krona_taxonomy_directory`. Therefore if you do not supply the a Krona directory, Krona plots will not be produced for MALT.

##### MetaPhlAn4

MetaPhlAn4 currently does not accept FASTA files as input, therefore no output will be produced for these input files.

##### mOTUs

mOTUs currently does not accept FASTA files as input, therefore no output will be produced for these input files.

#### Post Processing

##### Visualisation

nf-core/taxprofiler supports generation of Krona interactive pie chart plots for the following compatible tools.

- Kraken2
- Centrifuge
- Kaiju
- MALT

> ⚠️ MALT KRONA plots cannot be generated automatically, you must also specify a Krona taxonomy directory with `--krona_taxonomy_directory` if you wish to generate these.

##### Multi-Table Generation

In addition to per-sample profiles, the pipeline also supports generation of 'native' multi-sample taxonomic profiles (i.e., those generated by the taxonomic profiling tools themselves or additional utility scripts provided by the tool authors).

These are executed on a per-database level. I.e., you will get a multi-sample taxon table for each database you provide for each tool and will be placed in the same directory as the directories containing the per-sample profiles.

The following tools will produce multi-sample taxon tables:

- **Bracken** (via bracken's `combine_bracken_outputs.py` script)
- **Centrifuge** (via KrakenTools' `combine_kreports.py` script)
- **Kaiju** (via Kaiju's `kaiju2table` tool)
- **Kraken2** (via KrakenTools' `combine_kreports.py` script)
- **MetaPhlAn4** (via MetaPhlAn's `merge_metaphlan_tables.py` script)
- **mOTUs** (via the `motus merge` command)

Note that the multi-sample tables from these folders are not inter-operable with each other as they can have different formats.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/taxprofiler
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/taxprofiler releases page](https://github.com/nf-core/taxprofiler/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

#### For beginners

A first step to bypass this error, you could try to increase the amount of CPUs, memory, and time for the whole pipeline. Therefor you can try to increase the resource for the parameters `--max_cpus`, `--max_memory`, and `--max_time`. Based on the error above, you have to increase the amount of memory. Therefore you can go to the [parameter documentation of rnaseq](https://nf-co.re/rnaseq/3.9/parameters) and scroll down to the `show hidden parameter` button to get the default value for `--max_memory`. In this case 128GB, you than can try to run your pipeline again with `--max_memory 200GB -resume` to skip all process, that were already calculated. If you can not increase the resource of the complete pipeline, you can try to adapt the resource for a single process as mentioned below.

#### Advanced option on process level

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so, based on the search results, the file we want is `modules/nf-core/star/align/main.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB.
Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB.
The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

> **NB:** We specify the full process name i.e. `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

### Updating containers (advanced users)

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

     ```nextflow
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.0.5'
         }
     }
     ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Tutorials

### Retrieving databases or building custom databases

Not all taxonomic profilers provide ready-made or default databases. Here we will give brief guidance on how to build custom databases for each supported taxonomic profiler.

You should always consult the documentation of each tool for more information, as here we only provide short minimal-tutorials as quick reference guides (with no guarantee they are up to date).

The following tutorials assumes you already have the tool available (e.g. installed locally, or via conda, docker etc.), and you have already downloaded the FASTA files you wish to build into a database.

#### Bracken custom database

Bracken does not require an independent database nor not provide any default databases for classification/profiling, but rather builds upon Kraken2 databases. See [Kraken2](#kraken2-custom-database) for more information on how to build these.

In addition to a Kraken2 database, you also need to have the (average) read lengths (in bp) of your sequencing experiment, the K-mer size used to build the Kraken2 database, and Kraken2 available on your machine.

```bash
bracken-build -d <KRAKEN_DB_DIR> -k <KRAKEN_DB_KMER_LENGTH> -l <READLENGTH>
```

> 🛈 You can speed up database construction by supplying the threads parameter (`-t`).

> 🛈 If you do not have Kraken2 in your `$PATH` you can point to the binary with `-x /<path>/<to>/kraken2`.

<details markdown="1">
<summary>Expected files in database directory</summary>

- `bracken`
  - `hash.k2d`
  - `opts.k2d`
  - `taxo.k2d`
  - `database.kraken`
  - `database100mers.kmer_distrib`
  - `database100mers.kraken`
  - `database150mers.kmer_distrib`
  - `database150mers.kraken`

</details>

You can follow Bracken [tutorial](https://ccb.jhu.edu/software/bracken/index.shtml?t=manual) for more information.

#### Centrifuge custom database

To build a custom Centrifuge database, a user needs to download taxonomy files, make a custom `seqid2taxid.map` and combine the fasta files together.

In total, you need four components: a tab-separated file mapping sequence IDs to taxonomy IDs (`--conversion-table`), a tab-separated file mapping taxonomy IDs to their parents and rank, up to the root of the tree (`--taxonomy-tree`), a pipe-separated file mapping taxonomy IDs to a name (`--name-table`), and the reference sequences.

An example of custom `seqid2taxid.map`:

```
 NC_001133.9 4392
 NC_012920.1 9606
 NC_001134.8 4392
 NC_001135.5 4392
```

```bash
centrifuge-download -o taxonomy taxonomy
cat *.{fa,fna} > input-sequences.fna
centrifuge-build -p 4 --conversion-table seqid2taxid.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp input-sequences.fna taxprofiler_cf
```

<details markdown="1">
<summary>Expected files in database directory</summary>

- `centrifuge`
  - `<database_name>.<number>.cf`
  - `<database_name>.<number>.cf`
  - `<database_name>.<number>.cf`
  - `<database_name>.<number>.cf`

</details>

For the Centrifuge custom database documentation, see [here](https://ccb.jhu.edu/software/centrifuge/manual.shtml#custom-database).

#### DIAMOND custom database

To create a custom database for DIAMOND, the user should download and unzip the NCBI's taxonomy files and the input FASTA files.

The download and build steps are as follows:

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
unzip taxdmp.zip

## warning: large file!
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz

## warning: takes a long time!
cat ../raw/*.faa | diamond makedb -d testdb-diamond --taxonmap prot.accession2taxid.FULL.gz --taxonnodes nodes.dmp --taxonnames names.dmp

## clean up
rm *dmp *txt *gz *prt *zip
```

<details markdown="1">
<summary>Expected files in database directory</summary>

- `diamond`
  - `<database_name>.dmnd`

</details>

A detailed description can be found [here](https://github.com/bbuchfink/diamond/wiki/1.-Tutorial)

#### Kaiju custom database

To build a kaiju database, you need three components: a FASTA file with the protein sequences ,the NCBI taxonomy dump files, and you need to define the uppercase characters of the standard 20 amino acids you wish to include.

> ⚠️ The headers of the protein fasta file must be numeric NCBI taxon identifiers of the protein sequences.

To download the NCBI taxonomy files, please run the following commands:

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
unzip new_taxdump.zip
```

To build the database, run the following command (the contents of taxdump must be in the same location where you run the command):

```bash
kaiju-mkbwt -a ACDEFGHIKLMNPQRSTVWY -o proteins proteins.faa
kaiju-mkfmi proteins
```

> 🛈 You can speed up database construction by supplying the threads parameter (`-t`).

<details markdown="1">
<summary>Expected files in database directory</summary>

- `kaiju`
  - `kaiju_db_*.fmi`
  - `nodes.dmp`
  - `names.dmp`

</details>

For the Kaiju database construction documentation, see [here](https://github.com/bioinformatics-centre/kaiju#custom-database).

#### Kraken2 custom database

To build a Kraken2 database you need two components: a taxonomy (consisting of `names.dmp`, `nodes.dmp`, and `*accession2taxid`) files, and the FASTA files you wish to include.

To pull the NCBI taxonomy, you can run the following:

```bash
kraken2-build --download-taxonomy --db <YOUR_DB_NAME>
```

You can then add your FASTA files with the following build command.

```bash
kraken2-build --add-to-library *.fna --db <YOUR_DB_NAME>
```

You can repeat this step multiple times to iteratively add more genomes prior building.

Once all genomes are added to the library, you can build the database (and optionally clean it up):

```bash
kraken2-build --build --db <YOUR_DB_NAME>
kraken2-build --clean --db <YOUR_DB_NAME>
```

You can then add the `<YOUR_DB_NAME>/` path to your nf-core/taxprofiler database input sheet.

<details markdown="1">
<summary>Expected files in database directory</summary>

- `kraken2`
  - `opts.k2d`
  - `hash.k2d`
  - `taxo.k2d`

</details>

You can follow the Kraken2 [tutorial](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#custom-databases) for a more detailed description.

#### KrakenUniq custom database

For any KrakenUniq database, you require: taxonomy files, the FASTA files you wish to include, a `seqid2mapid` file, and a k-mer length.

First you must make a `seqid2taxid.map` file which is a two column text file containing the FASTA sequence header and the NCBI taxonomy ID for each sequence:

```
MT192765.1  2697049
```

Then make a directory (`<DB_DIR_NAME>/`), containing the `seqid2taxid.map` file, and your FASTA files in a subdirectory called `library/` (these FASTA files can be symlinked). You must then run the `taxonomy` command on the `<DB_DIR_NAME>/` directory, and then build it.

```bash
mkdir -p <DB_DIR_NAME>/library
mv `seqid2taxid.map` <DB_DIR_NAME>/
mv *.fna  <DB_DIR_NAME>/library
krakenuniq-download --db <DB_DIR_NAME>  taxonomy
krakenuniq-build --db <DB_DIR_NAME> --kmer-len 31
```

> 🛈 You can speed up database construction by supplying the threads parameter (`--threads`) to `krakenuniq-build`.

<details markdown="1">
<summary>Expected files in database directory</summary>

- `krakenuniq`
  - `opts.k2d`
  - `hash.k2d`
  - `taxo.k2d`
  - `database.idx`
  - `taxDB`

</details>

Please see the [KrakenUniq documentation](https://github.com/fbreitwieser/krakenuniq#database-building) for more information.

#### MALT custom database

To build a MALT database, you need the FASTA files to include, and an (unzipped) [MEGAN mapping 'db' file](https://software-ab.informatik.uni-tuebingen.de/download/megan6/) for your FASTA type. In addition to the input directory, output directory, and the mapping file database, you also need to specify the sequence type (DNA or Protein) with the `-s` flag.

```bash
malt-build -i <path>/<to>/<fasta>/*.{fna,fa,fasta} -a2t <path>/<to>/<map>.db -d <YOUR_DB_NAME>/  -s DNA
```

You can then add the `<YOUR_DB_NAME>/` path to your nf-core/taxprofiler database input sheet.

⚠️ MALT generates very large database files and requires large amounts of RAM. You can reduce both by increasing the step size `-st` (with a reduction in sensitivity).

> 🛈 MALT-build can be multi-threaded with `-t` to speed up building.

<details markdown="1">
<summary>Expected files in database directory</summary>

- `malt`
  - `ref.idx`
  - `taxonomy.idx`
  - `taxonomy.map`
  - `index0.idx`
  - `table0.idx`
  - `table0.db`
  - `ref.inf`
  - `ref.db`
  - `taxonomy.tre`

</details>

See the [MALT manual](https://software-ab.informatik.uni-tuebingen.de/download/malt/manual.pdf) for more information.

#### MetaPhlAn4 custom database

MetaPhlAn4 does not allow (easy) construction of custom databases. Therefore we recommend to use the prebuilt database of marker genes that is provided by the developers.

To do this you need to have `MetaPhlAn4` installed on your machine.

```bash
metaphlan --install --bowtie2db <YOUR_DB_NAME>/
```

You can then add the `<YOUR_DB_NAME>/` path to your nf-core/taxprofiler database input sheet.

> 🛈 It is generally not recommended to modify this database yourself, thus this is currently not supported in the pipeline. However, it is possible to customise the existing database by adding your own marker genomes following the instructions [here](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.1#customizing-the-database).

> 🖊️ If using your own database is relevant for you, please contact the nf-core/taxprofiler developers on the [nf-core slack](https://nf-co.re/join) and we will investigate supporting this.

<details markdown="1">
<summary>Expected files in database directory</summary>

- `metaphlan4`
  - `mpa_v30_CHOCOPhlAn_201901.pkl`
  - `mpa_v30_CHOCOPhlAn_201901.pkl`
  - `mpa_v30_CHOCOPhlAn_201901.fasta`
  - `mpa_v30_CHOCOPhlAn_201901.3.bt2`
  - `mpa_v30_CHOCOPhlAn_201901.4.bt2`
  - `mpa_v30_CHOCOPhlAn_201901.1.bt2`
  - `mpa_v30_CHOCOPhlAn_201901.2.bt2`
  - `mpa_v30_CHOCOPhlAn_201901.rev.1.bt2`
  - `mpa_v30_CHOCOPhlAn_201901.rev.2.bt2`
  - `mpa_latest`

</details>

More information on the MetaPhlAn4 database can be found [here](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.1#installation).

#### mOTUs custom database

mOTUs does not provide the ability to construct custom databases. Therefore we recommend to use the the prebuilt database of marker genes provided by the developers.

To do this you need to have `mOTUs` installed on your machine.

```bash
motus downloadDB
```

Then supply the `db_mOTU/` path to your nf-core/taxprofiler database input sheet.

> ⚠️ The `db_mOTU/` directory may be downloaded to somewhere in your Python's `site-package` directory. You will have to find this yourself as the exact location varies depends on installation method.

More information on the mOTUs database can be found [here](https://motu-tool.org/installation.html).

## Troubleshooting and FAQs

### I get a warning during centrifuge_kreport process with exit status 255

When a sample has insufficient hits for abundance estimation, the resulting `report.txt` file will be empty.

When trying to convert this to a kraken-style report, the conversion tool will exit with a status code `255`, and provide a `WARN`.

This is **not** an error nor a failure of the pipeline, just your sample has no hits to the provided database when using centrifuge.
