rule fastqc:
	input:
		get_fastqc
	output:
		"qc/fastqc/{raw}_fastqc.zip"
	shell:
		"""
		fastqc {input} -o qc/fastqc
		"""

rule stats:
	input:
		"results_{ref}/mapping/{raw}.bam"
	output:
		"qc/stats/{ref}:{raw}"
	shell:
		"""
		samtools stats {input} > {output}
		"""

rule flagstat:
	input:
		"results_{ref}/mapping/{raw}.bam"
	output:
		"qc/flagstat/{ref}:{raw}"
	shell:
		"""
		samtools flagstat {input} > {output}
		"""

from scipy.stats import gaussian_kde
rule coverage_qc:
	input:
		"results_{ref}/coverage/{raw}.cov"
	output:
		"qc/coverage/{ref}:{raw}_mqc.txt"
	run:
		cov = pd.read_table(input[0], names=["Chr", "Start","End", "Name", "Cov"])
		with open(output[0], "w") as f:
			model = gaussian_kde(np.log10(cov["Cov"]+1))
			x = np.arange(1,6, step=0.1)
			f.write(assets["coverage_qc"]+"\n")
			for a, b in zip(x, model(x)):
				f.write(f"{a}\t{b}\n")
		
rule ontarget_qc:
	input:
		"results_{ref}/coverage/{raw}.cov",
		"results_{ref}/coverage/{raw}.bed"
	output:
		"qc/ontarget/{ref}:{raw}_mqc.txt"
	run:
		cov = pd.read_table(input[0], names=["Chr", "Start","End", "Name", "Cov"])
		mapped = pd.read_table(input[1], names=["Chr", "Start","End", "Name"]).shape[0]
		ontarg = sum(cov['Cov'])
		with open(output[0], "w") as f:
			f.write(assets["ontarget_qc"]+"\n")
			f.write(f"On-target\t{ontarg}\n")
			f.write(f"Off-target\t{mapped - ontarg}\n")

rule length_qc:
	input:
		"results_{ref}/coverage/{raw}.bed"
	output:
		"qc/length/{ref}:{raw}_mqc.txt"
	run:
		bed = pd.read_table(input[0], names=["Chr", "Start","End", "Name"])
		with open(output[0], "w") as f:
			f.write(assets["length_qc"]+"\n")
			ys, xs = np.histogram(bed["End"] - bed["Start"], bins=20)
			for x, y in zip(xs, ys):
				f.write(f"{x}\t{y}\n")

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