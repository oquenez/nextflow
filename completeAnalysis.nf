#!/usr/bin/env nextflow
reads_ch	= Channel.fromFilePairs(params.reads)
kit_file	= file(params.kit)
srcFolder	= file(params.srcFolder)
seqCenter	= params.seqCenter
project		= params.projectName


process 'align_sample' {
        tag "${sampleId}"
	cpus params.cpuBwa
	
	input :
		set sampleId, file(reads) from reads_ch
		
	output :
		set sampleId, file("${sampleId}_bwa-mem-P-M.sam")  into sam_ch
	script :
	"""
		bwa mem -t ${params.cpuBwa} -M $params.genome $reads -R "@RG\tID:bwa\tSM:${sampleId}\tCN:${seqCenter}\tPU:unknown\tLB:unknown\tPL:ILLUMINA" > ${sampleId}_bwa-mem-P-M.sam
	"""
}

process 'sortSam' {
	tag "${sampleId}"
	cpus params.cpuPicard
	
	input:
                set sampleId, file(sam) from sam_ch
        output:
                set sampleId, file("${sampleId}.sorted.bam")  into sorted_ch
        script:
        """
		java -XX:ParallelGCThreads=${params.cpuPicard} -Xmx${params.ramPicard}g -jar ${params.picardD}/SortSam.jar INPUT=$sam OUTPUT=${sampleId}.sorted.bam SORT_ORDER=coordinate TMP_DIR=TMP
	"""

}

process 'dedup_sample' {
	tag "${sampleId}"
	
	cpus params.cpuPicard
	
	input:
                set sampleId, file(bam) from sorted_ch
        output:
                set sampleId, file("${sampleId}.sorted.dedup.bam")  into dedup_ch
        script:
	"""
		java -XX:ParallelGCThreads=${params.cpuPicard} -Xmx${params.ramPicard}g -jar ${params.picardD}/MarkDuplicates.jar INPUT=$bam OUTPUT=${sampleId}.sorted.dedup.bam METRICS_FILE=${sampleId}.metrics ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT TMP_DIR=TMP
	"""
}

process 'buildIndex_sample' {
        tag "${sampleId}"

        cpus params.cpuPicard

	input:
                set sampleId, file(bam) from dedup_ch
	output:
		set sampleId, file(bam), file("${sampleId}.sorted.dedup.bai") into buildRecall_ch

	script:
	"""
		java -XX:ParallelGCThreads=${params.cpuPicard} -Xmx${params.ramPicard}g -jar ${params.picardD}/BuildBamIndex.jar INPUT=$bam

	"""
}

process 'realign_sample'{
	tag "${sampleId}"
	cpus params.cpuGATK

	input:
		set sampleId, file(bam), file(bai) from buildRecall_ch
	output:
		set sampleId, file("${sampleId}.sorted.dedup.withRG.real.bam"), file("${sampleId}.sorted.dedup.withRG.real.bai") into realign_ch
	script:
	"""
		java -Xmx4g -jar ${params.GATK} -T IndelRealigner -nt ${params.cpuGATK} -R $params.genome -I $bam -targetIntervals ${params.intervals} -known ${params.mills} -known ${params.dbsnp} -o ${sampleId}.sorted.dedup.withRG.real.bam
	"""
}

process 'buildBQSR_table' {
	tag "${sampleId}"

	cpus params.cpuGATK

	input:
                set sampleId, file(bam), file(bai) from realign_ch
	output:
		set sampleId, file(bam), file (bai), file("${sampleId}.BQSR.table") into recall_table_ch
	script:
	"""
		java -Xmx4g -jar ${params.GATK} -T BaseRecalibrator -nt 8 -R $params.genome -I $bam -knownSites ${params.dbsnp} -knownSites ${params.mills} -knownSites ${params.indels} -o ${sampleId}.BQSR.table
	""" 
}

