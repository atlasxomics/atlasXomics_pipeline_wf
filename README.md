# ATX epigenomic preprocessing
The pipeline composed of three steps:

##### Filtering 
In this step we use [bbmap](https://github.com/BioInfoTools/BBMap/blob/master/sh/bbduk.sh) to identify and subset the reads with primers.

##### Alignment 
The most recent and fastest alignment tools for epigenomics data [chromap](https://github.com/haowenz/chromap) is applied in this step. Chromap will process barcodes from fastq_R2.gz and  align reads of both fastq_R1.gz and fastq_R2.gz files. The result would be a fragment file which can be used for downstream analysis by [ArchR](https://www.archrproject.com/) or [Signac](https://stuartlab.org/signac/).

##### Statistics 
The quality control indices like FRIP, TSS and peak files prepared in this step using [pycisTopic](https://github.com/aertslab/pycisTopic) package.
