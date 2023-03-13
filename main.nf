#!/usr/bin/env nextflow



nextflow.enable.dsl = 2


//check some variables before execution


process PRINT_VERSIONS {
    publishDir "$params.outdir/software", mode: "copy"

    output:
    path("versions.txt")
    """
    echo "bwa-mem2: 2.2.1" > versions.txt
    echo "samtools: 1.16.1" >> versions.txt
    echo "elprep: 5.1.3" >> versions.txt
    echo "fastqc: v0.12.1" >> versions.txt
    echo "qualimap: v.2.2.2-dev" >> versions.txt
    echo "mostdepth: 0.3.3" >> versions.txt
    """
}



process FASTQC {
    tag "$sampleId-mem"
    //label 'process_medium'
    publishDir "$params.outdir/QC/FASTQC", mode: "copy"

    input:
    tuple val(sampleId), file(read1), file(read2)

    output:
    path("${sampleId}.fastqc"), emit: fqc 

    script:
    if(params.debug == true){
    """
    echo fastqc -o ${sampleId}.fastqc $read1 $read2
    mkdir -p ${sampleId}.fastqc
    touch ${sampleId}.fastqc/report.fastqc
    
    """
    } else{
    """
     fastqc -t $task.cpus -o ${sampleId}.fastqc $read1 $read2
    """
    }
    
}

//we do run bwa-mem2 or bwa mem
process BWAMEM {
    
    tag "$sampleId-mem"
    //label 'process_high'
    publishDir "$params.outdir/BWA", mode: "copy", pattern: '*.log.*'
    publishDir "$params.outdir/BWA/HLA", mode: "copy", pattern: '*.hla.all'

    input:
    tuple val(sampleId), file(read1), file(read2)

    output:
    tuple val("${sampleId}"), file("${sampleId}.aln.bam"), emit: bams
    path("${sampleId}.log.bwamem") 
    path("${sampleId}.hla.all") , optional: true
    path("${sampleId}.log.hla") , optional: true
 
   script:
    def aln="bwa-mem2" 
    //we define the aln tool
    if(params.aligner=="bwa"){
	aln="bwa"
    }	
    if(params.debug == true){
    """
    echo "seqtk mergepe $read1 $read2 | ${aln} mem -p -t $task.cpus -R'@RG\tID:${sampleId}\tSM:ILL' ${params.ref} - 2 > ${sampleId}.log.bwamem | k8 bwa-postalt.js -p ${sampleId}.hla ${params.ref}.alt | samtools view -1 - > ${sampleId}.aln.bam; "
    echo "run-HLA ${sampleId}.hla > ${sampleId}.hla.top 2> ${sampleId}.log.hla;i"
    echo "touch ${sampleId}.hla.HLA-dummy.gt; cat ${sampleId}.hla.HLA*.gt | grep ^GT | cut -f2- > ${sampleId}.hla.all"
    echo "rm -f ${sampleId}.hla.HLA*;"
    touch ${sampleId}.aln.bam
    touch ${sampleId}.log.bwamem
    touch ${sampleId}.hla.all
    """
    }else{
    if(params.hla == "true"){
    """
	seqtk mergepe ${reads}  \
  	| ${aln} mem -p -t $task.cpus -R'@RG\tID:${sampleId}\tSM:ILL' ${ref} - 2 > ${sampleId}.log.bwamem \
  	| k8 bwa-postalt.js -p ${sampleId}.hla ${ref}.alt \
  	| samtools view -1 - > ${sampleId}.aln.bam;
	run-HLA ${sampleId}.hla > ${sampleId}.hla.top 2> ${sampleId}.log.hla;
	touch ${sampleId}.hla.HLA-dummy.gt; cat ${sampleId}.hla.HLA*.gt | grep ^GT | cut -f2- > ${sampleId}.hla.all;
	rm -f ${sampleId}.hla.HLA*;
    """
    }
    else if (params.alt == "true"){
     """
	seqtk mergepe ${reads}  \
  	| ${aln} mem -p -t $task.cpus  -R'@RG\tID:${sampleId}\tSM:ILL' ${ref} - 2 > ${sampleId}.log.bwamem \
  	| k8 bwa-postalt.js -p ${sampleId}.hla hs38DH.fa.alt \
  	| samtools view -1 - > ${sampleId}.aln.bam
     """	
    } else{
	//normal mapping mode
     """
	seqtk mergepe ${reads}  \
  	| ${aln} mem -p -t $task.cpus  -R'@RG\tID:${sampleId}\tSM:ILL' ${ref} - 2 > ${sampleId}.log.bwamem \
  	| samtools view -1 - > ${sampleId}.aln.bam
     """	
    }
    }

}

//we do preproces the bam file
process ELPREP {

    tag "$sampleId-elprep"
    //label 'process_high'
    //we save some stats from elprep
    publishDir "$params.outdir/ELPREP", mode: "copy", pattern: '*.output.*'

    input:
     tuple val(sampleId), file(bam)
    output:
    tuple val("${sampleId}"), file("${sampleId}.out.bam"), file("${sampleId}.out.bam.bai"), emit: bams
    //path("${sampleId}.out.bam.bai") , emit : bindex
    path("${sampleId}.output.metrics"), emit : metrics
    path("${sampleId}.output.recal"), optional: true

    script:
    
    if(params.debug == true){
    """
    echo elprep sfm ${bam} ${sampleId}.out.bam --mark-duplicates --mark-optical-duplicates ${sampleId}.output.metrics \
        --sorting-order coordinate \
        --bqsr  ${sampleId}.output.recal --known-sites ${params.dbsnp},${params.dbindel} \
        --reference ${params.elpre_ref}  --nr-of-threads $task.cpus
    echo samtools index ${sampleId}.out.bam
    #output files
    touch ${sampleId}.out.bam
    touch ${sampleId}.out.bam.bai
    touch ${sampleId}.output.metrics
    touch ${sampleId}.output.recal
    """
    }else{
    if(params.bqsr == true){
    """
    elprep sfm ${bam} ${sampleId}.out.bam --mark-duplicates --mark-optical-duplicates ${sampleId}.output.metrics \
        --sorting-order coordinate \
        --bqsr  ${sampleId}.output.recal --known-sites ${params.dbsnp},${params.dbindel} \
        --reference ${params.elpre_ref}  --nr-of-threads $task.cpus
    samtools index ${sampleId}.out.bam
    """
   }
   else{
    """
    elprep sfm ${bam} ${sampleId}.out.bam --mark-duplicates --mark-optical-duplicates ${sampleId}.output.metrics \
        --sorting-order coordinate \
        --reference ${params.elpre_ref}  --nr-of-threads $task.cpus
    #we index the resulting bam file
    samtools index ${sampleId}.out.bam
    """
    }
   }
}


process QUALIMAP{
    tag "$sampleId-qualimap"
    //label 'process_medium'
    
    publishDir "$params.outdir/QC/QUALIMAP", mode: "copy"

    input:
    tuple val(sampleId), file(bam), file(bai)
    

    output:
    path("${sampleId}.qualimap") , emit : qc
    
    script:
    if(params.debug == true){
    """
    echo qualimap  bamqc  -bam $bam  -outdir ${sampleId}.qualimap --java-mem-size=30G -nt $task.cpus
    mkdir ${sampleId}.qualimap
    touch ${sampleId}.qualimap/summaryQualimap.txt
    """
    }else{
    """
    qualimap  bamqc  -bam $bam  -outdir ${sampleId}.qualimap --java-mem-size=30G -nt $task.cpus
    """
    }
    
}

process B2C{
     tag "$sampleId-b2c"
    //label 'process_medium'
    
    publishDir "$params.outdir/CRAM", mode: "copy"

    input:
    tuple val(sampleId), file(bam), file(bai)
    

    output:
    tuple val("${sampleId}"), file("${sampleId}.out.cram"), file("${sampleId}.out.cram.crai"), emit: crams
    path("${sampleId}.out.flagstat"), emit: flags
    
    script:
    if(params.debug == true){
    """
    echo samtools view -C -T ${params.ref} $bam -o ${sampleId}.out.cram
    echo samtools index ${sampleId}.out.cram
    echo samtools flagstats  ${sampleId}.out.cram > ${sampleId}.out.flagstat
    touch ${sampleId}.out.cram 
    touch ${sampleId}.out.cram.crai
    touch ${sampleId}.out.flagstat
    """
    }else{
    """
    samtools view -C -T ${params.ref} $bam -o ${sampleId}.out.cram -@ $task.cpus
    samtools index ${sampleId}.out.cram
    samtools flagstats  ${sampleId}.out.cram > ${sampleId}.out.flagstat
    """
    }   
}

process DEPTH{
    tag "$sampleId-depth"
    //label 'process_medium'
    
    publishDir "$params.outdir/DEPTH", mode: "copy"

    input:
    tuple val(sampleId), file(cram), file(crai)
    

    output:
    path("${sampleId}.depth.*"), emit: depth
    
    script:
    if(params.debug == true){
    """
    echo mosdepth -t $task.cpus ${sampleId}.depth $cram
    touch ${sampleId}.depth.mosdepth.dist.txt
    touch ${sampleId}.depth.mosdepth.summary.txt
    """
    }else{
    """
    mosdepth -t $task.cpus ${sampleId}.depth $cram
    """
    }   
}



//we declare the workflow for index

workflow {
    // TODO do a parameter check
    PRINT_VERSIONS()
    //we read pairs from regex 
    if(params.reads != null){
    // --reads "./reads/B087*.merge.{1,2}.fq.gz"
    read_pairs_ch = channel.fromFilePairs(params.reads)
    }else if(params.csv != null){
    //we reads pairs from csv
    read_pairs_ch=Channel.fromPath(params.csv) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleId, file(row.read1), file(row.read2)) }
    }else{
        println "Error: reads regex or path"
    }

    read_pairs_ch.view()
    //ref = path(params.ref)
    //fastqc read quality
    FASTQC(read_pairs_ch)
    //read aligment alt/hla
    BWAMEM(read_pairs_ch)
    //bam procesisng sort/duplciates/bqrs
    ELPREP(BWAMEM.out.bams)
    //Quality of alignments
    QUALIMAP(ELPREP.out.bams)
    //BAM->CRAM conversion
    B2C(ELPREP.out.bams)
    //Coverage Stats
    DEPTH(B2C.out.crams)
}
