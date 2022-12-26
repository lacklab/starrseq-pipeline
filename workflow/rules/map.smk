rule map_bwa:
	input:
		get_fqs
	output:
		temp("results_{ref}/mapping/{raw}.raw.bam")
	params:
		idx = config["REF"]["BWA_IDX"]
	threads:
		32
	shell:
		"""
		bwa mem -t {threads} {params.idx} {input} \
		| samtools view -bS - > {output}
		"""

rule bam_process:
	input:
		"results_{ref}/mapping/{raw}.raw.bam"
	output:
		temp("results_{ref}/mapping/{raw}.coorsorted.bam")
	threads:
		32
	params:
		config["OUTPUT"]["BAMPROCESS_PARAMS"]
	shell:
		"""
		samtools view -h {params} {input} \
		| samtools fixmate -m -@ {threads} - - \
		| samtools sort -@ {threads} -m 10G - \
		| samtools markdup -@ {threads} - {output}
		"""

rule bam_filter:
	input:
		"results_{ref}/mapping/{raw}.coorsorted.bam"
	output:
		temp("results_{ref}/mapping/{raw}.filtered.bam")
	threads:
		32
	params:
		config["REF"]["FA"],
		get_filter_p
	shell:
		"""
		samtools view {input} | egrep -v "chrM" | \
		samtools view -b -@ {threads} -T {params} > {output}
		"""

rule bam_merge:
	input:
		get_reps
	output:
		"results_{ref}/mapping/{raw}.final.bam"
	threads:
		32
	run:
		if str(input).find(' ') != -1:
			shell("""
			samtools merge -@ {threads} -o {output} {input}
			samtools index {output}
			""")
		else:
			shell("""
			mv {input} {output}
			samtools index {output}
			""")

rule bam_nsort:
	input:
		"results_{ref}/mapping/{raw}.final.bam"
	output:
		"results_{ref}/mapping/{raw}.ns.bam"
	threads:
		16
	shell:
		"""
		samtools sort -n -@ {threads} -m 10G {input} -o {output}
		"""