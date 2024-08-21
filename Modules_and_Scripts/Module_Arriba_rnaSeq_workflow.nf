#!/usr/bin/env nextflow

/*
##################################################################################################
###################################### Module Arriba #############################################
##################################################################################################
*/

nextflow.enable.dsl=2

/*

#################################### Module STAR_Arriba ##########################################

    This module performs the fusions detection with the Arriba tool starting as input from the subdirectory returned by the CreateDir process (Module_OrderFiles_FastQC_rnaSeq_workflow.nf) and returns the subdirectory containing the results of the process.		

*/

if (params.help) {
	log.info ''
	exit 0
}    

process Star_Arriba_module {
	
conda '/home/.conda/envs/ipx-arriba'

	errorStrategy 'ignore'
	
	input:
	val directory
	
	output:
	stdout
	
	label 'Heavy_Tools'
		
	script:
	"""
	R1=\$(ls \$(echo $directory)/*_R1_*.fastq.gz)
	R2=\$(ls \$(echo $directory)/*_R2_*.fastq.gz)
	path=\$(basename \$(echo \${R1}) | awk -F "_R1" '{print \$1}')
	
	mkdir -p \$(echo $directory)/Arriba_results
	
	STAR	--runThreadN 8\
		--genomeDir ${params.Ensembl}\
		--genomeLoad NoSharedMemory\
		--readFilesIn \${R1} \${R2}\
		--readFilesCommand zcat\
		--outStd BAM_Unsorted\
		--outSAMtype BAM Unsorted\
		--outSAMunmapped Within\
		--outBAMcompression 0\
		--outFilterMultimapNmax 50\
		--peOverlapNbasesMin 10\
		--alignSplicedMateMapLminOverLmate 0.5\
		--alignSJstitchMismatchNmax 5 -1 5 5\
		--chimSegmentMin 10\
		--chimOutType WithinBAM HardClip\
		--chimJunctionOverhangMin 10\
		--chimScoreDropMax 30\
		--chimScoreJunctionNonGTAG 0\
		--chimScoreSeparation 1\
		--chimSegmentReadGapMax 3\
		-chimMultimapNmax 50|
	arriba	-x /dev/stdin\
		-o "\$(echo $directory)/Arriba_results/\${path}_fusions.tsv"\
		-O "\$(echo $directory)/Arriba_results/\${path}_fusions_discarded.tsv"\
		-a "${params.Ensembl}/${params.Ensembl_fasta}"\
		-g "${params.Ensembl}/${params.Ensembl_gtf}"\
		-b "${params.Arriba_data}/${params.Arriba_blacklist}"\
		-k "${params.Arriba_data}/${params.Arriba_known_fusions}"\
		-t "${params.Arriba_data}/${params.Arriba_known_fusions}"\
		-p "${params.Arriba_data}/${params.Arriba_protein_domains}" > "\$(echo $directory)/Arriba_results/Arriba_stdout.txt"
		
	echo \${path}
	"""
}

