# Snakefile

# CONFIGURATION
# Set the input directory name - change this as needed
INPUT_DIR = "fastq_files"

# Auto-detect samples from the input directory - change this as needed
# note: Snakemake doesn't use shell wildcards like '*' for inputs
# check rule fastqc, remove_human, summarize_human_read_loss
SAMPLES,_ = glob_wildcards(f"{INPUT_DIR}/{{sample}}_S{{unit}}_R1_001.fastq.gz")
# DEBUG: Print detected samples
print(f"DEBUG: Detected {len(SAMPLES)} samples: {SAMPLES}")

# Set the project name for output files
PROJECT_NAME = "fl119_igaseq_NOVA1284P"

# Define GB for downstream resource allocation. Do not change this
GB = 1024

# Snakemake works backwards. Specify the final files you want to check for
# if you don't specify all of the required output files, some rules will not run
rule all:
    input:
        # This will trigger fastqc + multiqc
        "fastqc_summary/multiqc_report.html",
    
########################################
# Step 0: fastqc
########################################
rule fastqc:
    input: 
        r1=f"{INPUT_DIR}/{{sample}}_S{{unit}}_R1_001.fastq.gz",
        r2=f"{INPUT_DIR}/{{sample}}_S{{unit}}_R2_001.fastq.gz"
    output:
        html1="fastqc_output/{sample}_R1_fastqc.html",
        zipped1="fastqc_output/{sample}_R1_fastqc.zip",
        html2="fastqc_output/{sample}_R2_fastqc.html",
        zipped2="fastqc_output/{sample}_R2_fastqc.zip"
    threads: 2
    resources:
        mem_mb=8*GB,
        runtime=60
    shell:
        """ 
        source /etc/profile.d/modules.sh
        module load fastqc/0.12.1
        fastqc -o fastqc_output -f fastq -t 2 {input.r1} {input.r2}
        """
########################################
# Step 0.1: multiqc
########################################
rule multiqc:
    input:
        expand("fastqc_output/{sample}_R1_fastqc.zip", sample=SAMPLES),
        expand("fastqc_output/{sample}_R2_fastqc.zip", sample=SAMPLES)
    output:
        "fastqc_summary/multiqc_report.html"
    shell:
        """
        source /etc/profile.d/modules.sh
        module load multiqc/1.23
        multiqc fastqc_output -o fastqc_summary
        """