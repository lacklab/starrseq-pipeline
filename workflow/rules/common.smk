import pandas as pd
import numpy as np
#import yaml
#from pathlib import Path
#config = yaml.safe_load(Path("config/config.yaml").read_text())

# `samples` includes 'IP/CONTROL' pairs 
samples = pd.read_table(config["SAMPLES"])

# `units` includes all files that need to be preprocessed
units = pd.read_table(config["UNITS"])
units["Raw"] = units["Name"] + "_" + units["Unit"].astype(str)

ref = config["REF"]["NAME"]
# >>> utils >>>
def get_lib(wildcards):
	return units.loc[units["Raw"] == wildcards.raw, "Library"].unique()[0]

def get_units(wildcards):
	return units.loc[units["Name"] == wildcards.raw, "Raw"].unique()

def get_fq1(wildcards):
	return units.loc[units["Raw"] == wildcards.raw, "Fastq1"].unique()[0]

def get_fq2(wildcards):
	return units.loc[units["Raw"] == wildcards.raw, "Fastq2"].unique()[0]

def get_contol(wildcards):
	return samples.loc[samples["Name"] == wildcards.raw, "Control"].unique()[0]
# <<< utils <<<

# >>> `qc.smk` >>>
fastqc_map = {}
for u in units["Fastq1"].tolist() + units["Fastq2"].tolist():	
	if (u.find("gz") != -1) or (u.find("zip") != -1):
		fastqc_map[u.rsplit(".", 2)[0].split("/")[-1]] = u
	elif u.find("SRR") != -1:
		fastqc_map[u] = f"sra-data/{u}_1.fastq.gz"
	elif u != "-":
		fastqc_map[u.rsplit(".", 1)[0].split("/")[-1]] = u

def get_fastqc(wildcards):
	return fastqc_map[wildcards.raw]

def get_multiqc(wildcards):
	out = []
	for i, row in units.iterrows():
		lib = row["Library"]
		fq1 = row["Fastq1"].split("/")[-1]
		fq2 = row["Fastq2"].split("/")[-1]
		if fq1.find("gz") != -1:
			fq1 = fq1.rsplit(".", 2)[0]
			fq2 = fq2.rsplit(".", 2)[0]
		if lib == "Single":
			out.append(f"qc/{fq1}_fastqc.zip")
		elif lib == "Paired":
			out.append(f"qc/{fq1}_fastqc.zip")
			out.append(f"qc/{fq2}_fastqc.zip")
		out.append(f"qc/{ref}:{row['Raw']}.raw.bam_flagstat")
		out.append(f"qc/{ref}:{row['Name']}.final.bam_flagstat")
		out.append(f"qc/{ref}:{row['Name']}.final.bam_stats")
		out.append(f"qc/{ref}:{row['Name']}_coverage_mqc.txt")
		out.append(f"qc/{ref}:{row['Name']}_ontarget_mqc.txt")
		out.append(f"qc/{ref}:{row['Name']}_length_mqc.txt")
	return expand(out)
# <<< `qc.smk` <<<

# >>> `map.smk` functions >>>
def get_fqs(wildcards):
	name, unit = wildcards.raw.rsplit("_",1)
	fq1 = units.loc[units["Name"] == name,"Fastq1"].unique()[0]
	source = str(fq1).find("SRR") != -1
	lib = get_lib(wildcards)
	if source:
		srr = fq1
		if lib == "Single":
			return f"sra-data/{srr}_1.fastq.gz"
		elif lib == "Paired":
			return f"sra-data/{srr}_1.fastq.gz", f"sra-data/{srr}_2.fastq.gz"
	else:
		fq1 = get_fq1(wildcards)
		fq2 = get_fq2(wildcards)
		if lib == "Single":
			return fq1
		elif lib == "Paired":
			return fq1, fq2

def get_filter_p(wildcards):
	param = "-F 2048 "
	lib = get_lib(wildcards)
	if lib == "Single":
		param += ""
	elif lib == "Paired":
		param += "-f 2"
	return param
def get_reps(wildcards):
	reps = get_units(wildcards)
	return expand("results_{{ref}}/mapping/{rep}.filtered.bam", rep=reps)
# <<< `map.smk` functions <<<

# >>> `peak.smk` functions >>>
def get_macs_p(wildcards):
	lib = units.loc[units["Name"] == wildcards.raw, "Library"].unique()[0]
	input_c = get_contol(wildcards)
	param = f"-t results_{wildcards.ref}/mapping/{wildcards.raw}.final.bam "
	if not input_c == "-":
		param += f" -c results_{wildcards.ref}/mapping/{input_c}.final.bam "
	if lib == "Single":
		param += " -f BAM"
	elif lib == "Paired":
		param += " -f BAMPE"
	return param

def get_macs_i(wildcards):
	inputs = []
	input_c = get_contol(wildcards)
	if not input_c == "-":
		inputs.append(f"results_{wildcards.ref}/mapping/{input_c}.final.bam")
	inputs.append(f"results_{wildcards.ref}/mapping/{wildcards.raw}.final.bam")
	return inputs
# >>> `peak.smk` functions <<<

# >>> `count.smk` >>>
def get_multicov_i(wildcards):
	out = [] 
	for i, row in units.iterrows():
		out.append(f"results_{ref}/mapping/{row['Name']}.final.bam")
	return out
# <<< `count.smk` <<<

bwNorm = config["OUTPUT"]["BW_NORMALIZATIONS"]

outputs = []

if config["OUTPUT"]["RUN"]["QC"]:
	outputs += ["qc/multiqc_report.html"]

q = config['OUTPUT']['MACS_THRESHOLD']
if config["OUTPUT"]["RUN"]["PEAKS"]:
	outputs += [
		f"results_{ref}/peaks/{raw}_{q}_peaks.narrowPeak"
		for raw in samples["Name"]
	]

if config["OUTPUT"]["RUN"]["BWS"]:
	outputs += [
		f"results_{ref}/bigwig/{raw}.genomecov.{norm}.bw"
		for raw in samples["Name"]
		for norm in bwNorm
	]

if config["OUTPUT"]["RUN"]["COUNTTABLE"]:
	outputs += [f"results_{ref}/capture_counts.tsv"]

if config["OUTPUT"]["RUN"]["COVERAGE"]:
	outputs += [
		f"results_{ref}/coverage/{raw}.cov"
		for raw in samples["Name"]
]


assets = {}
with open("assets/coverage_multiqc.asset", "r") as f:
	assets["coverage_qc"]	= ""
	for line in f.readlines():
		assets["coverage_qc"] += line

with open("assets/ontarget_multiqc.asset", "r") as f:
	assets["ontarget_qc"]	= ""
	for line in f.readlines():
		assets["ontarget_qc"] += line
		
with open("assets/length_multiqc.asset", "r") as f:
	assets["length_qc"]	= ""
	for line in f.readlines():
		assets["length_qc"] += line
		