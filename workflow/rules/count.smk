rule multicov:
	input:
		get_multicov_i
	output:
		"results_{ref}/capture_counts.tsv"
	params:
		config["CAPTURE_REGIONS"]
	shell:
		"""
		bedtools multicov -bams {input} -bed {params} > {output}
		"""

rule coverage:
	input:
		"results_{ref}/mapping/{raw}.ns.bam"
	output:
		"results_{ref}/coverage/{raw}.cov"
	params:
		config["CAPTURE_REGIONS"]
	shell:
		"""
		bedtools bamtobed -bedpe -i {input} \
		| cut -f 1,2,6,7 \
		| bedtools intersect -a {params} -b - -c > {output}
		"""
