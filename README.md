## Introduction

This is a bioinformatics analysis pipeline used for shotgun metagenomic data developed at Zymo Research. This pipeline was adpated from community-developed [nf-core/taxprofiler](https://github.com/nf-core/taxprofiler) pipeline version 1.0.0. Many changes were made to the original pipeline. Some are based on our experience or preferences. But more importantly, we want to make the pipeline and its results easier to use/understand by people without bioinformatics experience. People can run the pipeline on the point-and-click bioinformatics platform [Aladdin Bioinformatics](https://www.aladdin101.org). Changes include but are not limited to:
* Changed the behavior of pipeline so that the user may choose one taxonomy profiler instead of running all available taxonomy profilers. We found that some of the profilers have worse performances or outdated databases. We have disabled those profilers temporarily, but kept the code to run them, with a plan to add them if necessary. This is a philosophical change, we believe this approach offers simplicity and avoids confusion for researchers.
* Added [sourmash](https://github.com/sourmash-bio/sourmash) as the preferred taxonomy profiler.
* Added a Zymo version of the sourmash database that includes common host genomes so that host removal step does not need to be run anymore.
* Upgraded MetaPhlAn3 to MetaPhlAn4.
* Added antimicrobial resistance analysis with [AMRplusplus](https://github.com/Microbial-Ecology-Group/AMRplusplus/tree/master).
* Added visualizations of sourmash and MetaPhlAn4 results to the report.
* Added diversity analysis using [Qiime2](https://qiime2.org/) and corresponding visualizations to the report.
* Fixed, simplified, and improved the report.
* Made the pipeline more resistant to bad samples, so that they don't stop the processing of others.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies.

## Pipeline summary

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) or [`falco`](https://github.com/smithlabcode/falco) as an alternative option)
2. Performs optional read pre-processing (*code for long-read inherited from nf-core/taxprofiler, but not separately tested by us yet*)
   - Adapter clipping and merging (short-read: [`fastp`](https://github.com/OpenGene/fastp), [`AdapterRemoval2`](https://github.com/MikkelSchubert/adapterremoval); long-read: [`porechop`](https://github.com/rrwick/Porechop))
   - Low complexity and quality filtering (short-read: [`bbduk`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/), [`PRINSEQ++`](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus); long-read: [`Filtlong`](https://github.com/rrwick/Filtlong))
3. Perform Host-read removal
   - Host-read removal (short-read: [`BowTie2`](http://bowtie-bio.sourceforge.net/bowtie2/); long-read: [`Minimap2`](https://github.com/lh3/minimap2)). This is not performed when `sourmash-zymo` is selected as the database because it already contains host sequences. 
   - Statistics for host-read removal ([`Samtools`](http://www.htslib.org/))
4. Run merging when applicable
5. Identifies antimicrobial resistance genes in samples from database [MEGARes version 3](https://academic.oup.com/nar/article/51/D1/D744/6830666)
   - Reads are aligned to MEGARes reference ([`bwa mem`](https://bio-bwa.sourceforge.net/bwa.shtml))
   - Resistome statistics are quantified and compiled for each sample ([`AMRplusplus`](https://github.com/Microbial-Ecology-Group/AMRplusplus/tree/master))
7. Performs taxonomic profiling using one of: (*nf-core/taxprofiler has more choices for this step, if there are tools you'd like for this step, please let us know.*)
   - [`sourmash`](https://github.com/sourmash-bio/sourmash)
   - [`MetaPhlAn4`](https://huttenhower.sph.harvard.edu/metaphlan/)
8. Merge all taxonomic profiling results into one table ([`Qiime2`](https://qiime2.org/)) and generate composition plot ([`Krona`](https://github.com/marbl/Krona))
9. Perform alpha/beta diversity analysis ([`Qiime2`](https://qiime2.org/)) and beta diversity differential expression ([`ANCOM-BC`](https://www.bioconductor.org/packages/release/bioc/html/ANCOMBC.html))
10. Present all results in above steps in a report ([`MultiQC`](http://multiqc.info/))

## Quick Start

We recommend you run this pipeline via the [Aladdin Bioinformatics platform](https://www.aladdin101.org). It is much easier to run without any requirement for coding. Also, because the Zymo sourmash database is private, public users would not be able to use it via the command line. If you would still like to run the pipeline via CLI, see instruction below.

### Prerequisites
* [Nextflow](https://www.nextflow.io) version 22.10.1 or later
* At least 8 CPU threads and 60GB memory. The default config file [nextflow.config](./nextflow.config) has higher settings, please modify to fit your device.
* [Docker](https://www.docker.com/) if you are running locally.
* Permissions to AWS S3 and Batch resources if you are running on AWS Batch.

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
1. The parameter `--design` is required. It must be a CSV file with the following format.
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
2. The parameter `--database` is used to change taxonomy profiler and database. It has a default value 'sourmash-zymo'. You can skip this if you don't want to change it.
3. The parameter `--run_amr` is used to run antimicrobial resistance analysis. This parameter is by default false. If you wish to skip this analysis, remove `--run_amr` from the command line or set it to "false".
4. The parameters `--awsregion`, `--awsqueue`, `-work-dir`, and `--outdir` are required when running on AWS Batch, the latter two must be directories on S3.
5. The parameter `-r` will run a specific release of the pipeline. If not specified, it will run the the `main` branch instead.
6. The parameter `-name` will define the title of the MultiQC report.

There are many other options built in the pipeline to customize your run and handle specific situations, please refer to the [Usage Documentation](docs/usage.md).

### Using Docker locally
```bash
nextflow run Zymo-Research/aladdin-shotgun \
    -profile docker,dev \
    --design "<path to design CSV file>" \
    --database sourmash-zymo-2024
```
Please see above for requirements of the design CSV file.

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
Using SLURM on ZymoCloud is only supported after release 0.0.13
## Credits

This pipeline was adapted from nf-core/taxprofiler version 1.0.0. Please refer to [credits](https://github.com/nf-core/taxprofiler#credits) for list of orginal contributors. Contributors from Zymo Research include:
- Zymo Research microbiomics team (source code, database, review)
- [Nora Sharp](https://github.com/nsharp2) (Pipeline coding)
- [Zhenfeng Liu](https://github.com/zxl124) (Pipeline coding)
