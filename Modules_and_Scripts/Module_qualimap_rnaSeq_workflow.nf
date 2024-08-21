#!/usr/bin/env nextflow

/*
###################################################################################################
######################################## Module qualimap ##########################################
###################################################################################################
*/

nextflow.enable.dsl = 2

/*

######################################## Module qualimap ##########################################

    This module contains a process that performs the quality control with the qualimap tool of the BAM file resulting from the STAR alignment.

*/

if (params.help) {
	log.info ''
	exit 0
}    

process qualimap_module {

conda '/home/.conda/envs/qualimap_conda'

	errorStrategy 'ignore'
	
	input:
	val bam	
	
	output:
	stdout
	
	script:
	"""
	echo $bam
	dir_to_files=\$(dirname \$(echo $bam))
	
	if [[ \$(basename \$(echo $bam)) != "pseudoalignments.bam" ]]; then
	        path=\$(basename \$(echo $bam) | awk -F "_split" '{print \$1}')
	        bam_al=\${dir_to_files}/\${path}_aligned_merged_sorted.bam
	else
        	bam_al=\$(echo $bam)
	fi
	
	qualimap rnaseq -bam \${bam_al} -gtf ${params.Ensembl}/${params.Ensembl_gtf} -outdir \${dir_to_files}/qualimap_result --java-mem-size=15G
	"""
}

