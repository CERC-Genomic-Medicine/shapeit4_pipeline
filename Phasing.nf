#!/usr/bin/env nextflow

/*
* AUTHOR: CERC Genomic Medicine,  Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2023
*/

process get_chunks {
	//executor "local"
	cache "lenient"
	cpus 1
	memory "4GB"

	input:
	tuple path(vcf), path(vcf_index)

        output:
        tuple path("*.chunk"), path(vcf), path(vcf_index)

	"""
        chrom=`bcftools index -s ${vcf} | cut -f1`
        start_bp=`bcftools view -HG ${vcf} | head -n1 | cut -f2`
	stop_bp=`bcftools index -s ${vcf} | cut -f2`
       
       	extend=0
	for i in `seq \${start_bp} ${params.window} \${stop_bp}`; do
		if [ \${extend} -eq 0 ]; then
			chunk_start=\$((${params.flank} > i ? 0 : i - ${params.flank}))
		fi
		chunk_stop=\$((i + ${params.window} + ${params.flank}))
		n=`bcftools view -HG ${vcf} \${chrom}:\${chunk_start}-\${chunk_stop} | wc -l`
		if [ \${n} -gt 0 ]; then
			printf "\${chrom}\t\${chunk_start}\t\${chunk_stop}\t\${n}\n" > \${chrom}_\${chunk_start}_\${chunk_stop}.chunk
			extend=0
		else
			extend=1
		fi
	done
	if [ \${extend} -eq 1 ]; then
		printf "\${chrom}\t\${chunk_start}\t\${chunk_stop}\t\${n}\n" > \${chrom}_\${chunk_start}_\${chunk_stop}.chunk
	fi
	"""
}


process phase_chunks {

	cache "lenient"
        scratch true
        //errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return "retry" }
	//maxRetries 3
	cpus 1
	memory "16 GB"
	time "48h"

	input:
	tuple path(chunk), path(vcf), path(vcf_index)

	output:
	tuple path("*.phased.bcf"), path("*.phased.log")

	publishDir "results/logs", pattern: "*.phased.log"

	"""
        chrom=`head -n1 ${chunk} | cut -f1`
        start_bp=`head -n1 ${chunk} | cut -f2`
	stop_bp=`head -n1 ${chunk} | cut -f3`
	${params.shapeit_exec} --sequencing --input ${vcf} --region \${chrom}:\${start_bp}-\${stop_bp} --map ${params.shapeit_maps}/\${chrom}.b38.gmap.gz --output \${chrom}_\${start_bp}_\${stop_bp}.phased.bcf --log \${chrom}_\${start_bp}_\${stop_bp}.phased.log
	"""
}


process ligate {
	
	cache "lenient"
	scratch true
	cpus 1
	memory "16G"
	time "24h"

	input:
	tuple val(chrom), path(phased_vcfs)

	output:
	path "${chrom}.phased.bcf*"

	publishDir "results"

	"""
	for f in ${phased_vcfs}; do bcftools index \${f}; done
	for f in ${phased_vcfs}; do echo \${f}; done | sort -V > files_list.txt
	bcftools concat -f files_list.txt -l -Ob -o ${chrom}.phased.bcf
	bcftools index ${chrom}.phased.bcf
	"""
}

workflow {
	vcfs = Channel.fromPath(params.unphased_vcfs).map { file -> [file, file + (file.getExtension() == "bcf" ? ".csi" : ".tbi")] }
	vcfs_chunks = get_chunks(vcfs)
	phased_vcfs_chunks = phase_chunks(vcfs_chunks.transpose())
	ligate(phased_vcfs_chunks.map { it -> [ it[0].getSimpleName().split('_')[0], it[0]]  }.groupTuple())
}
