# aln

A nextflow pipeline for alignment of short WGS reads.

## Preparing index sequences

- prepare reference index for bwa-mem2 
  ```
   bwa-mem2 index ref.fa
  ```

- prepare elprep files using the folowing commands:
  ```
  #reference file
  elprep fasta-to-elfasta ref.fa ref.fa.elfasta
  #variats files GTAK4_Bunddle for variant calibration
  elprep vcf-to-elsites <vcf-file> <elsites-file>
  ``` 

## Conda Environment

We will build a micromamba environment with the needed software

```
micromamba create -f conda.yml 
```

## Param files

Here we provide an example to map the short reads to hg38, the context of the params files is the following


```
dbsnp: /mnt/beegfs/labs/DiGenomaLab/databases/references/human/GATK_Bundle/Homo_sapiens_assembly38.dbsnp138.elsites
dbindel: /mnt/beegfs/labs/DiGenomaLab/databases/references/human/GATK_Bundle/Mills_and_1000G_gold_standard.indels.hg38.elsites
ref: /mnt/beegfs/labs/DiGenomaLab/databases/references/human/bwa2/hs38DH.fa
elpre_ref: /mnt/beegfs/labs/DiGenomaLab/databases/references/human/hs38DH.fa.elfasta
alt_js: /mnt/beegfs/home/adigenova/micromamba/envs/aln/bin/bwa-postalt.js
bqsr: true
```
Save the above content in a file (i.e) : aln-params.yml 

## Read file
provide a csv file wiht the following information:

```
sampleId,read1,read2
test2,./test_reads/test2.R1.fq.gz,./test_reads/test2.R2.fq.gz
test3,./test_reads/test3.R1.fq.gz,./test_reads/test3.R2.fq.gz
test,./test_reads/test.R1.fq.gz,./test_reads/test.R2.fq.gz
```
Currently if a sample is split into several files is necesary to merge the reads before runing the pipeline.


Save the above content in a file (i.e) : reads.csv
## runnig the pipeline

```
nextflow run main.nf --csv reads.cvs -profile uoh -params-file aln-params.yml
```
in case of failure use:

```
nextflow run main.nf --csv reads.cvs -profile uoh -params-file aln-params.yml -resume
```

that will generate directory called ***results***

## Creating an agrregated result

To create an aggregated report across all the samples, is possible to run multiqc on the result directory:

- load the environment 

  ```
  micromamba activate aln
  ```
- run multiqc

 ```
  multiqc .
 ```

