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
* [Reference databases](#reference-databases)
   * [`--database`](#--database)
* [Trimming and filtering](#trimming-and-filtering)
   * [`--preprocessing_qc_tool`](#--preprocessing_qc_tool)
   * [`--save_preprocessed_reads`](#--save_preprocessed_reads)
   * [`--shortread_qc_tool`](#--shortread_qc_tool)
   * [`--shortread_qc_skipadaptertrim`](#--shortread_qc_skipadaptertrim)
   * [`--shortread_qc_adapter1`](#--shortread_qc_adapter1)
   * [`--shortread_qc_adapter2`](#--shortread_qc_adapter2)
   * [`--shortread_qc_adapterlist`](#--shortread_qc_adapterlist)
   * [`--shortread_qc_mergepairs`](#--shortread_qc_mergepairs)
   * [`--shortread_qc_includeunmerged`](#--shortread_qc_includeunmerged)
   * [`--shortread_qc_minlength`](#--shortread_qc_minlength)
   * [`--shortread_complexityfilter_tool`](#--shortread_complexityfilter_tool)
   * [`--shortread_complexityfilter_entropy`](#--shortread_complexityfilter_entropy)
   * [`--shortread_complexityfilter_bbduk_windowsize`](#--shortread_complexityfilter_bbduk_windowsize)
   * [`--shortread_complexityfilter_bbduk_mask`](#--shortread_complexityfilter_bbduk_mask)
   * [`--shortread_complexityfilter_fastp_threshold`](#--shortread_complexityfilter_fastp_threshold)
   * [`--shortread_complexityfilter_prinseqplusplus_mode`](#--shortread_complexityfilter_prinseqplusplus_mode)
   * [`--shortread_complexityfilter_prinseqplusplus_dustscore`](#--shortread_complexityfilter_prinseqplusplus_dustscore)
   * [`--save_complexityfiltered_reads`](#--save_complexityfiltered_reads)
* [Preprocessing Long Read QC Options](#preprocessing-long-read-qc-options)
   * [`--perform_longread_qc`](#--perform_longread_qc)
   * [`--longread_qc_skipadaptertrim`](#--longread_qc_skipadaptertrim)
   * [`--longread_qc_skipqualityfilter`](#--longread_qc_skipqualityfilter)
   * [`--longread_qc_qualityfilter_minlength`](#--longread_qc_qualityfilter_minlength)
   * [`--longread_qc_qualityfilter_keeppercent`](#--longread_qc_qualityfilter_keeppercent)
   * [`--longread_qc_qualityfilter_targetbases`](#--longread_qc_qualityfilter_targetbases)
* [Preprocessing Host Removal Options](#preprocessing-host-removal-options)
   * [`--perform_shortread_hostremoval`](#--perform_shortread_hostremoval)
   * [`--perform_longread_hostremoval`](#--perform_longread_hostremoval)
   * [`--hostremoval_reference`](#--hostremoval_reference)
   * [`--shortread_hostremoval_index`](#--shortread_hostremoval_index)
   * [`--longread_hostremoval_index`](#--longread_hostremoval_index)
   * [`--save_hostremoval_index`](#--save_hostremoval_index)
   * [`--save_hostremoval_bam`](#--save_hostremoval_bam)
   * [`--save_hostremoval_unmapped`](#--save_hostremoval_unmapped)
* [Preprocessing Run Merging Options](#preprocessing-run-merging-options)
   * [`--perform_runmerging`](#--perform_runmerging)
   * [`--save_runmerged_reads`](#--save_runmerged_reads)
* [Profiling Options](#profiling-options)
   * [`--centrifuge_save_reads`](#--centrifuge_save_reads)
   * [`--diamond_output_format`](#--diamond_output_format)
   * [`--diamond_save_reads`](#--diamond_save_reads)
   * [`--kaiju_taxon_rank`](#--kaiju_taxon_rank)
   * [`--kraken2_save_reads`](#--kraken2_save_reads)
   * [`--kraken2_save_readclassification`](#--kraken2_save_readclassification)
   * [`--kraken2_save_minimizers`](#--kraken2_save_minimizers)
   * [`--krakenuniq_save_reads`](#--krakenuniq_save_reads)
   * [`--krakenuniq_ram_chunk_size`](#--krakenuniq_ram_chunk_size)
   * [`--krakenuniq_save_readclassification`](#--krakenuniq_save_readclassification)
   * [`--malt_mode`](#--malt_mode)
   * [`--malt_save_reads`](#--malt_save_reads)
   * [`--malt_generate_megansummary`](#--malt_generate_megansummary)
   * [`--sourmash_kmersize`](#--sourmash_kmersize)
   * [`--sourmash_threshold_bp`](#--sourmash_threshold_bp)
   * [`--sourmash_trim_low_abund`](#--sourmash_trim_low_abund)
   * [`--motus_use_relative_abundance`](#--motus_use_relative_abundance)
   * [`--motus_save_mgc_read_counts`](#--motus_save_mgc_read_counts)
   * [`--motus_remove_ncbi_ids`](#--motus_remove_ncbi_ids)
* [Visualization and diversity options](#visualization-and-diversity-options)
   * [`--lowread_filter`](#--lowread_filter)
   * [`--min_frequency`](#--min_frequency)
   * [`--min_samples`](#--min_samples)
   * [`--qiime_tax_agglom_min`](#--qiime_tax_agglom_min)
   * [`--qiime_tax_agglom_max`](#--qiime_tax_agglom_max)
   * [`--ancombc_fdr_cutoff`](#--ancombc_fdr_cutoff)
   * [`--top_taxa`](#--top_taxa)
   * [`--group_of_interest`](#--group_of_interest)
   * [`--run_profile_standardisation`](#--run_profile_standardisation)
   * [`--standaridisation_motus_generation`](#--standaridisation_motus_generation)
   * [`--run_krona`](#--run_krona)
   * [`--krona_taxonomy_directory`](#--krona_taxonomy_directory)
   * [`--standardisation_taxpasta_format`](#--standardisation_taxpasta_format)
   * [`--taxpasta_taxonomy_dir`](#--taxpasta_taxonomy_dir)
   * [`--taxpasta_add_name`](#--taxpasta_add_name)
   * [`--taxpasta_add_rank`](#--taxpasta_add_rank)
   * [`--taxpasta_add_lineage`](#--taxpasta_add_lineage)
   * [`--taxpasta_add_idlineage`](#--taxpasta_add_idlineage)
* [AMR Options](#amr-options)
   * [`--run_amr`](#--run_amr)
   * [`--resistome_threshold`](#--resistome_threshold)
   * [`--amr_index_files`](#--amr_index_files)
* [Skipping options](#skipping-options)
   *  [`--skip_heatmap`](#--skip_heatmap)
   *  [`--skip_alpha_rarefaction`](#--skip_alpha_rarefaction)
   *  [`--skip_alphadiversity`](#--skip_alphadiversity)
   *  [`--skip_individalpha`](#--skip_individalpha)
   *  [`--skip_betadiversity`](#--skip_betadiversity)
* [Max job request options](#max-job-request-options)
   * [`--max_cpus`](#--max_cpus)
   * [`--max_memory`](#--max_memory)
   * [`--max_time`](#--max_time)
* [Internal options](#internal-options)
   * [`--multiqc_config`](#--multiqc_config)
   * [`--max_multiqc_email_size`](#--max_multiqc_email_size)
   * [`--validate_params`](#--validate_params)
   * [`--show_hidden_params`](#--show_hidden_params)
   * [`--tracedir`](#--tracedir)
   * [`--publish_dir_mode`](#--publish_dir_mode)
   * [`--multiqc_logo`](#--multiqc_logo)
   * [`--igenomes_base`](#--igenomes_base)
* [Generic options](#generic-options)
   * [`--help`](#--help)
   * [`--version`](#--version)
   * [`--email_on_fail`](#--email_on_fail)
   * [`--plaintext_email`](#--plaintext_email)
   * [`--monochrome_logs`](#--monochrome_logs)
   * [`--hook_url`](#--hook_url)
   * [`--multiqc_methods_description`](#--multiqc_methods_description)
   * [`--ignore_failed_samples`](#--ignore_failed_samples)
   * [`--report_name`](#--report_name)
   * [`--genome`](#--genome)
   * [`--igenomes_ignore`](#--igenomes_ignore)



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
Skip the removal of sequencing adapters.

This often can be useful to speed up run-time of the pipeline when analysing data downloaded from public databases such as the ENA or SRA, as adapters should already be removed (however we recommend to check FastQC results to ensure this is the case).

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

Entropy here corresponds to the amount of sequence variation exists within the read. Higher values correspond to more variety, and thus will likely result in more specific matching to a taxon's reference genome. The trade off here is fewer reads (or abundance information) available for having a confident identification.


> Modifies tool parameter(s):
> - BBDuk: `entropy=`
> - PRINSEQ++:  `-lc_entropy`
                    
### `--shortread_complexityfilter_bbduk_windowsize`
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
Turns on the ability to remove short-reads derived from a known organism, using Bowtie2 and samtools.

This subworkflow is useful to remove reads that may come from a host, or a known contamination like the human reference genome. Human DNA contamination of (microbial) reference genomes is well known, so removal of these prior profiling both reduces the risks of false positives, and in _some cases_ a faster runtime (as less reads need to be profiled).

Alternatively, you can include the reference genome within your profiling databases and can turn off this subworkflow, with the trade off of a larger taxonomic profiling database.

### `--perform_longread_hostremoval`
Turns on the ability to remove long-reads from the that derived from a known organism, using minimap2 and samtools

This subworkflow is useful to remove reads that may come from a host, or a known contamination like the human reference genome. Human DNA contamination of (microbial) reference genomes is well known, so removal of these prior profiling both reduces the risks of false positives, and in _some cases_ a faster runtime (as less reads need to be profiled).

Alternatively, you can include the reference genome within your profiling databases and can turn off this subworkflow, with the trade off of a larger taxonomic profiling database.

### `--hostremoval_reference`
Specify a path to the FASTA file (optionally gzipped) of the reference genome of the organism to be removed.

If you have two or more host organisms or contaminants you wish to remove, you can concatenate the FASTAs of the different taxa into a single one to provide to the pipeline.

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
Save the reads mapped to the host reference genome and off-target reads in BAM format as output by the respective hostremoval alignment tool.

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

## Visualization and diversity options

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

### `--top_taxa`
Specify the maximum number of top taxa from each sample you want to compare in QIIME_HEATMAP

### `--group_of_interest`
Group(s) of taxa that are of specific interest. Abundances of these taxa will be summarized in a separate file and plotted in a separate section in the report to highlight these taxa.

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

## AMR Options

### `--run_amr`
Whether or not to conduct analysis for antimicrobial resistance genes.

### `--resistome_threshold`
The minimum proportion of nucleotides that have aligned reads to for any AMR gene to be reported.

### `--amr_index_files`
Path to files related to AMR analysis.

## Skipping options

### `--skip_heatmap`
Remove off heatmap plot from MultiQC report

### `--skip_alpha_rarefaction`
Turn off alpha rarefaction analysis and remove alpha rarefaction plot from MultiQC report

### `--skip_alphadiversity`
Turn off alpha diversity analysis and remove alpha diversity plot from MultiQC report

### `--skip_individalpha`
Turn off the chart showing individual alpha diversity per sample.

### `--skip_betadiversity`
Turn off beta diversity analysis and remove beta diversity plot from MultiQC reports

## Max job request options

### `--max_cpus`
Maximum number of CPUs that can be requested for any single job. Use to set an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`

### `--max_memory`
Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`

### `--max_time`
Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`

## Internal options
### `--multiqc_config`
Custom config file to supply to MultiQC.

### `--max_multiqc_email_size`
File size limit when attaching MultiQC reports to summary emails.

### `--validate_params`
Boolean whether to validate parameters against the schema at runtime

### `--show_hidden_params`
By default, parameters set as _hidden_ in the schema are not shown on the command line when a user runs with `--help`. Specifying this option will tell the pipeline to show all parameters

### `--tracedir`
Directory to keep pipeline Nextflow logs and reports.

### `--publish_dir_mode`
The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.

### `--multiqc_logo`
Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file

### `--igenomes_base`
Directory / URL base for iGenomes references.

## Generic options
### `--help`
Display help text.

### `--version`
Display version and exit.

### `--email_on_fail`
An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.

### `--plaintext_email`
Send plain-text email instead of HTML.

### `--monochrome_logs`
Do not use coloured log outputs.

### `--hook_url`
Incoming hook URL for messaging service. Currently, MS Teams and Slack are supported.

### `--multiqc_methods_description`
Custom MultiQC yaml file containing HTML including a methods description.

### `--ignore_failed_samples`
Whether to ignore samples that fail QC and taxonomy profiling steps and carry on with good samples. Default True. Set to False to force the pipeline to abort when any samples fails.

### `--report_name`
Name for report title and report file name.

### `--genome`
If using a reference genome configured in the pipeline using iGenomes, use this parameter to give the ID for the reference. This is then used to build the full paths for all required reference genome files e.g. `--genome GRCh38`.

See the [nf-core website docs](https://nf-co.re/usage/reference_genomes) for more details.

### `--igenomes_ignore`
Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.
