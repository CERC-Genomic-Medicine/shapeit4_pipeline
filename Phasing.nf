

process chunk {
	scratch "~/scratch/"

	input:
	set file(vcf), file(index) from Channel.fromPath(params.unphased_vcfs).map{ vcf -> [ vcf, vcf + ".csi" ] }

        output:
        tuple file("*.chunk"), file(vcf), file(index) into chunks

	"""
        chrom=`bcftools view -HG ${vcf} | head -n1 | cut -f1`
        start_bp=`bcftools view -HG ${vcf} | head -n1 | cut -f2`
	stop_bp=`bcftools view -HG ${vcf} | tail -n1 | cut -f2`
       
	for i in `seq \${start_bp} ${params.window} \${stop_bp}`; do
		chunk_start=\$((${params.flank} > i ? 0 : i - ${params.flank}))
		chunk_stop=\$((i + ${params.window} + ${params.flank}))
		printf "\${chrom}\t\${chunk_start}\t\${chunk_stop}\n" > \${chrom}_\${chunk_start}_\${chunk_stop}.chunk
	done	
	"""
}


process phase {
	scratch "~/scratch/"

	input:
	set file(chunk), file(vcf), file(index) from chunks.transpose()

	output:
	tuple stdout, file("*.phased.bcf") into phased_chunks
	file "*.phased.log" into phased_logs

	publishDir "results/logs", pattern: "*.phased.log"


	"""
        chrom=`head -n1 ${chunk} | cut -f1`
        start_bp=`head -n1 ${chunk} | cut -f2`
	stop_bp=`head -n1 ${chunk} | cut -f3`
	${params.shapeit_exec} --input ${vcf} --region \${chrom}:\${start_bp}-\${stop_bp} --map ${params.shapeit_maps}/\${chrom}.b38.gmap.gz --output \${chrom}_\${start_bp}_\${stop_bp}.phased.bcf --log \${chrom}_\${start_bp}_\${stop_bp}.phased.log >/dev/null
	printf "\${chrom}"
	"""
}


process ligate {
	scratch "~/scratch/"

	input:
	set val(chrom), file(phased_bcfs) from phased_chunks.groupTuple()

	output:
	file "${chrom}.phased.bcf*" into phased

	publishDir "results"

	"""
	for f in ${phased_bcfs}; do bcftools index \${f}; done
	for f in ${phased_bcfs}; do echo \${f}; done | sort -V > files_list.txt
	bcftools concat -f files_list.txt -l -Ob -o ${chrom}.phased.bcf
	bcftools index ${chrom}.phased.bcf
	"""
}
