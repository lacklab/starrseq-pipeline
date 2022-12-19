rule fastqc:
	input:
		get_fastqc
	output:
		"qc/{raw}_fastqc.zip"
	shell:
		"""
		fastqc {input} -o qc
		"""

rule stats:
	input:
		"results_{ref}/mapping/{raw}.bam"
	output:
		"qc/{ref}:{raw}.bam_stats"
	shell:
		"""
		samtools stats {input} > {output}
		"""

rule flagstat:
	input:
		"results_{ref}/mapping/{raw}.bam"
	output:
		"qc/{ref}:{raw}.bam_flagstat"
	shell:
		"""
		samtools flagstat {input} > {output}
		"""

from scipy.stats import gaussian_kde
rule coverage_qc:
	input:
		"results_{ref}/coverage/{raw}.cov"
	output:
		"qc/{ref}:{raw}_coverage_mqc.txt"
	run:
		cov = pd.read_table(input[0], names=["Chr", "Start","End", "Cov"])
		model = gaussian_kde(np.log10(cov["Cov"]+1))
		x = np.arange(1,20, step=0.1)
		with open(output[0], "w") as f:
			f.write(assets["coverage_qc"]+"\n")
			for a, b in zip(x, model(x)):
				f.write(f"{a}\t{b}\n")

				
# TODO: This can be written in a script.	

rule multiqc:
	input:
		get_multiqc
	output:
		"qc/multiqc_report.html"
	shell:	
		"""
		cd qc/ && multiqc .
		"""

# TODO: get FRiP values and add to 'multiqc_data/multiqc_general_stats.txt'
# TODO: work on custom-content (https://multiqc.info/docs/#custom-content)
# TODO: QC output can include outputs; such as https://nf-co.re/chipseq/1.2.2/output