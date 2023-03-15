.DELETE_ON_ERROR:

dbsnp:=/mnt/beegfs/labs/DiGenomaLab/databases/references/human/GATK_Bundle/Homo_sapiens_assembly38.dbsnp138.elsites
indel:=/mnt/beegfs/labs/DiGenomaLab/databases/references/human/GATK_Bundle/Mills_and_1000G_gold_standard.indels.hg38.elsites
ref:=/mnt/beegfs/labs/DiGenomaLab/databases/references/human/hs38DH.fa

test.aln.bam:
	/mnt/beegfs/home/adigenova/micromamba/envs/aln/bin/seqtk mergepe reads/B0873IM_R1_FKDN230059971.merge.1.fq.gz  reads/B0873IM_R1_FKDN230059971.merge.2.fq.gz  \
  	| /mnt/beegfs/home/adigenova/micromamba/envs/aln/bin/bwa mem -p -t 14 -R'@RG\tID:sampe\tSM:bar' /mnt/beegfs/labs/DiGenomaLab/databases/references/human/hs38DH.fa - 2> test.log.bwamem \
  	| /mnt/beegfs/home/adigenova/micromamba/envs/aln/bin/k8 /mnt/beegfs/home/adigenova/micromamba/envs/aln/bin/bwa-postalt.js -p test.hla /mnt/beegfs/labs/DiGenomaLab/databases/references/human/hs38DH.fa.alt \
  	| /mnt/beegfs/home/adigenova/micromamba/envs/aln/bin/samtools view -1 - > test.aln.bam;
	/mnt/beegfs/home/adigenova/micromamba/envs/aln/bin/run-HLA test.hla > test.hla.top 2> test.log.hla;
	touch test.hla.HLA-dummy.gt; cat test.hla.HLA*.gt | grep ^GT | cut -f2- > test.hla.all;
	rm -f test.hla.HLA*;
test.out.bam: test.aln.bam
	elprep sfm $< $@ --mark-duplicates --mark-optical-duplicates test.output.metrics \
	--sorting-order coordinate \
	--bqsr test.output.recal --known-sites ${dbsnp},${indel} \
	--reference ${ref}.elfasta  --nr-of-threads 14


test.mosdepth.summary.txt:test.out.bam
	samtools index test.out.bam
	mosdepth -t 14 test test.out.bam
#bam2cram
test.out.cram:test.out.bam
	samtools view -C -T ${ref} $< -o $@
	samtools index $@
	samtools flagstats  $@

#depth
#qualimap
test.qualimap:test.out.bam
	qualimap  bamqc  -bam test.out.bam  -outdir test.qualimap --java-mem-size=30G -nt 4
	
all : test.mosdepth.summary.txt test.out.cram test.qualimap
