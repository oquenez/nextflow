profiles {
	standard {
		process.executor	= 'local'

		params.cpuBwa		= 8
		params.cpuGATK		= 16
		params.cpuPicard	= 8
		params.ramPicard	= 48

		params.picardD		= "/opt/picard-tools-1.108"
		params.GATK		= "/opt/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar"
		
		params.mills		= "/opt/REFERENCE/Mills_and_1000G_gold_standard.indels.b37.vcf"
		params.dbsnp		= "/opt/REFERENCE/dbsnp_141.vcf"
		params.indels		= "/opt/REFERENCE/1000G_phase1.indels.b37.vcf"
		params.intervals	= "/opt/REFERENCE/target_intervals.list"
		params.genome           = "/opt/REFERENCE/human_g1k_v37.fasta"
	}
}

