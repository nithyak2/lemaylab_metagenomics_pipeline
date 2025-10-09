# Snakefile
# for FL119 processing metagenomes 
###############################
# USAGE INSTRUCTIONS
# ------------------
# 1. Navigate to project root directory:
#    cd /path/to/project_root
# 
# 2. Get required files from MetaPhlAn github (this may need to be updated if databases change)
#    wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/sgb_to_gtdb_profile.py -O scripts/sgb_to_gtdb_profile.py
#    wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/util_fun.py -O scripts/util_fun.py
#    wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/mpa_vJan25_CHOCOPhlAnSGB_202503_SGB2GTDB.tsv -O scripts/mpa_vJan25_CHOCOPhlAnSGB_202503_SGB2GTDB.tsv
#    wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/merge_metaphlan_tables.py -O scripts/merge_metaphlan_tables.py  
#
# 3. Load snakemake v9.11.4 into a conda environment (if necessary):
#    eval "$(mamba shell hook --shell bash)"
#    mamba create -n snakemake_env -c conda-forge -c bioconda snakemake=9.11.4
#    conda activate snakemake_env
#
# 4. Install slurm executor plugin for snakemake v8+ (only needs to be done once):
#    pip install snakemake-executor-plugin-slurm
#
# 5. Quick check to make sure there are no errors (dry run):
#    snakemake -s scripts/Snakefile -n
#
# 6. Run the pipeline using one of these methods:
#
#    METHOD A - Submit via sbatch script (recommended):
#    sbatch scripts/submit_snakefile.sh
#
#    METHOD B - Run directly with snakemake:
#    snakemake -s scripts/Snakefile --executor slurm --jobs 20 --use-conda \
#        --default-resources slurm_account=dglemaygrp mem_mb=4096 runtime=600
#
# 7. Monitor progress:
#    tail -f logs/snakemake_<jobid>.out
###############################
# REQUIRED DIRECTORY STRUCTURE:
# Make sure the .py scripts from the MetaPhlAn github have execute permissions
# ls -l scripts/sgb_to_gtdb_profile.py 
# chmod +x scripts/sgb_to_gtdb_profile.py 
# 
# project_root/
# │
# ├── scripts/
# │   ├── Snakefile.py                    [this file]
# │   ├── sgb_to_gtdb_profile.py      [required script from MetaPhlAn github]
# │   ├── util_fun.py                 [required script from MetaPhlAn github]
# │   ├── mpa_vJan25_CHOCOPhlAnSGB_202503_SGB2GTDB.tsv [required for sgb to gtdb script, from MetaPhlAn github]
# │   ├── merge_metaphlan_tables.py   [required script from MetaPhlAn github]
# │   └── submit_snakefile.sh         [submit this file with sbatch]
# │     └──source ~/.bashrc
# │     └──conda activate snakemake_env 
# │     └──snakemake -s scripts/Snakefile --executor slurm --job 20 --use-conda --default-resources slurm_account=dglemaygrp mem_mb=4096 runtime=600 
# │
# └── {INPUT_DIR}/                 [REQUIRED - place input files here]
#     ├── sample1_R1.fastq.gz
#     ├── sample1_R2.fastq.gz
#     ├── sample2_R1.fastq.gz
#     └── sample2_R2.fastq.gz
#
# All output directories will be created automatically by Snakemake
#
###############################
#
# CONFIGURATION
# Set the input directory name - change this as needed
INPUT_DIR = "fastq_files"

# Auto-detect samples from the input directory - change this as needed
# note: Snakemake doesn't use shell wildcards like '*' for inputs
# check rule fastqc, remove_human, summarize_human_read_loss
SAMPLES, = glob_wildcards(f"{INPUT_DIR}/{{sample}}_S[^_]+_L[^_]+_R1_001.fastq.gz")
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
        # This will trigger the entire fastp pipeline
        expand("step2_fastp/merged/{sample}.merged.fastq.gz", sample=SAMPLES),
        
        # This will trigger fastqc + multiqc
        "fastqc_summary/multiqc_report.html",
        
        # This will trigger the summary reports
        "step2_fastp/pairedend/reports/pairedend_summary.txt",
        "step2_fastp/merged/reports/merged_summary.txt",
        "step1_bowtie2/human_read_loss_summary.txt",

        # This will trigger MetaPhlAn stepts
        f"step3_metaphlan/{PROJECT_NAME}_sgb.metaphlan4-2-2.txt",

########################################
# Step 0: fastqc
########################################
rule fastqc:
    input:
        r1=f"{INPUT_DIR}/{{sample}}_S*_L*_R1_001.fastq.gz",
        r2=f"{INPUT_DIR}/{{sample}}_S*_L*_R2_001.fastq.gz"
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
########################################
# Step 1: Remove human reads
# run bowtie2-build command if needed, only need to run once
# bowtie2-build /share/lemaylab-backedup/milklab/database/human_GRCh38_p13/GCF_000001405.39_GRCh38.p13_genomic.fna /share/lemaylab-backedup/nkumar/database/human_genome
########################################
rule remove_human:
    input:
        r1=f"{INPUT_DIR}/{{sample}}_S*_L*_R1_001.fastq.gz",
        r2=f"{INPUT_DIR}/{{sample}}_S*_L*_R2_001.fastq.gz"
    output:
        pe1="step1_bowtie2/{sample}_nohuman_pe.1.fastq",
        pe2="step1_bowtie2/{sample}_nohuman_pe.2.fastq",
        se="step1_bowtie2/{sample}_nohuman_se.fastq",
        stats="step1_bowtie2/{sample}_alignment_stats.txt"
    log:
        "logs/{sample}.remove_human.log"
    threads: 8
    params:
        index="/quobyte/dglemaygrp/BACKED-UP/nkumar/database/human_genome",
        basename="step1_bowtie2/{sample}_nohuman_pe.fastq"
    resources:
        mem_mb=64*GB,
        runtime=24*60
    shell:
        """
        source /etc/profile.d/modules.sh
        module load bowtie2/2.5.2
        # run bowtie2
        bowtie2 -q -p {threads} -x {params.index} \
            -1 {input.r1} -2 {input.r2} \
            --un-conc {params.basename} \
            --un {output.se} \
            1>/dev/null 2>{output.stats}
        """
########################################
# Step 1.1: Zip reads to save space
########################################
rule zip_reads:
    input:
        pe1="step1_bowtie2/{sample}_nohuman_pe.1.fastq",
        pe2="step1_bowtie2/{sample}_nohuman_pe.2.fastq",
        se="step1_bowtie2/{sample}_nohuman_se.fastq"
    output:
        pe1="step1_bowtie2/{sample}_nohuman_pe.1.fastq.gz",
        pe2="step1_bowtie2/{sample}_nohuman_pe.2.fastq.gz",
        se="step1_bowtie2/{sample}_nohuman_se.fastq.gz"
    shell:
        """ 
        source /etc/profile.d/modules.sh
        module load pigz 

        pigz {input.pe1} {input.pe2} {input.se}

        """
########################################
# Step 1.2: Summarize reads lost
########################################
rule summarize_human_read_loss:
    input:
        # Files after human read removal
        expand("step1_bowtie2/{sample}_nohuman_pe.1.fastq.gz", sample=SAMPLES),
        expand("step1_bowtie2/{sample}_nohuman_pe.2.fastq.gz", sample=SAMPLES),
        # Original input files
        expand(f"{INPUT_DIR}/{{sample}}_S*_L*_R1_001.fastq.gz", sample=SAMPLES),
        expand(f"{INPUT_DIR}/{{sample}}_S*_L*_R2_001.fastq.gz", sample=SAMPLES)
    output:
        "step1_bowtie2/human_read_loss_summary.txt"
    params:
        samples=" ".join(SAMPLES),  # Pass samples as a parameter
        input_dir=INPUT_DIR
    resources:
        mem_mb=4*GB,
        runtime=6*60
    run:
        import gzip
        import glob
        
        with open(output[0], "w") as out:
            out.write("sample\ttotal_reads\tremaining_reads\tpercent_lost\n")
            
            for sample in SAMPLES:
                # Use glob to find the actual files with wildcards
                r1_files = glob.glob(f"{params.input_dir}/{sample}_S*_L*_R1_001.fastq.gz")
                r2_files = glob.glob(f"{params.input_dir}/{sample}_S*_L*_R2_001.fastq.gz")
                
                if not r1_files or not r2_files:
                    continue
                    
                r1_in = r1_files[0]
                r2_in = r2_files[0]
                pe1 = f"step1_bowtie2/{sample}_nohuman_pe.1.fastq.gz"
                pe2 = f"step1_bowtie2/{sample}_nohuman_pe.2.fastq.gz"
                
                # Count total reads before filtering
                total_reads = 0
                for file in [r1_in, r2_in]:
                    with gzip.open(file, 'rt') as f:
                        total_reads += sum(1 for line in f) // 4
                
                # Count remaining reads after filtering
                remaining_reads = 0
                for file in [pe1, pe2]:
                    with gzip.open(file, 'rt') as f:
                        remaining_reads += sum(1 for line in f) // 4
                
                # Calculate percentage
                percent_lost = 100 * (total_reads - remaining_reads) / total_reads if total_reads > 0 else 0.0
                
                # Write results
                out.write(f"{sample}\t{total_reads}\t{remaining_reads}\t{percent_lost:.2f}\n")
########################################
# Step 2: Quality trimming with fastp
# This first pass saves paired end reads 
########################################
rule trim_quality:
    input:
        pe1="step1_bowtie2/{sample}_nohuman_pe.1.fastq.gz",
        pe2="step1_bowtie2/{sample}_nohuman_pe.2.fastq.gz"
    output:
        r1="step2_fastp/pairedend/{sample}_R1_paired.fastq.gz",
        r2="step2_fastp/pairedend/{sample}_R2_paired.fastq.gz",
        ur1="step2_fastp/pairedend/{sample}_R1_unpaired.fastq.gz",
        ur2="step2_fastp/pairedend/{sample}_R2_unpaired.fastq.gz",
        html="step2_fastp/pairedend/reports/{sample}.fastp.html",
        json="step2_fastp/pairedend/reports/{sample}.fastp.json"
    log:
        "logs/{sample}.fastp_pairedend.log"
    conda:
        "/quobyte/dglemaygrp/nkumar/conda_envs/fastp.yaml" 
    threads: 4
    resources:
        mem_mb=16*GB,
        runtime=6*60
    shell:
        """
        fastp --in1 {input.pe1} --in2 {input.pe2} \
            --thread {threads} --detect_adapter_for_pe \
            --trim_poly_g --dedup \
            --length_required 99 \
            --qualified_quality_phred 15 \
            --unpaired1 {output.ur1} --unpaired2 {output.ur2} \
            --out1 {output.r1} --out2 {output.r2} \
            --json {output.json} \
            --html {output.html}
        """
########################################
# Step 2.1: Merge paired reads with fastp
########################################
rule merge_reads:
    input:
        r1="step2_fastp/pairedend/{sample}_R1_paired.fastq.gz",
        r2="step2_fastp/pairedend/{sample}_R2_paired.fastq.gz",
    output:
        merged="step2_fastp/merged/{sample}.merged.fastq.gz",
        um1="step2_fastp/merged/{sample}_R1_unmerged.fastq.gz",
        um2="step2_fastp/merged/{sample}_R2_unmerged.fastq.gz",
        html="step2_fastp/merged/reports/{sample}.fastp.html",
        json="step2_fastp/merged/reports/{sample}.fastp.json"
    log:
        "logs/{sample}.fastp_merge.log"
    conda:
        "/quobyte/dglemaygrp/nkumar/conda_envs/fastp.yaml"
    resources:
        mem_mb=16*GB,
        runtime=6*60
    threads: 4
    shell:
        """
        fastp --in1 {input.r1} --in2 {input.r2} \
        --thread {threads} \
        --merge \
        --merged_out {output.merged} \
        --out1 {output.um1} --out2 {output.um2} \
        --overlap_len_require 10 \
        --overlap_diff_percent_limit 10 \
        --disable_trim_poly_g \
        --disable_adapter_trimming \
        --json {output.json} \
        --html {output.html}
        """
rule summarize_fastp_reports:
    input:
        pe_reports=expand("step2_fastp/pairedend/reports/{sample}.fastp.json", sample=SAMPLES),
        merged_reports=expand("step2_fastp/merged/reports/{sample}.fastp.json", sample=SAMPLES)
    output:
        pe_summary="step2_fastp/pairedend/reports/pairedend_summary.txt",
        merged_summary="step2_fastp/merged/reports/merged_summary.txt"
    run:
        import json
        # -----------------------------
        # Paired-end summary
        # -----------------------------
        with open(output.pe_summary, "w") as out_pe:
            # write header
            out_pe.write("sample_id\tbefore_filtering_total_reads\tafter_filtering_total_reads\tlow_quality_reads\ttoo_many_N_reads\ttoo_short_reads\tduplication_rate\tinsert_peak\tpercent_lost\n")
            
            for sample in SAMPLES:
                report_file = f"step2_fastp/pairedend/reports/{sample}.fastp.json"
                with open(report_file) as f:
                    data = json.load(f)

                # Extract fields from the JSON
                total_before = data["summary"]["before_filtering"]["total_reads"]
                total_after  = data["summary"]["after_filtering"]["total_reads"]
                low_qual     = data["filtering_result"]["low_quality_reads"]
                too_many_N   = data["filtering_result"]["too_many_N_reads"]
                too_short    = data["filtering_result"]["too_short_reads"]
                dup_rate     = data["duplication"]["rate"]
                insert_peak  = data["insert_size"]["peak"]
                percent_lost = 100 * (total_before - total_after) / total_before if total_before > 0 else 0

                # write line for this sample
                out_pe.write(f"{sample}\t{total_before}\t{total_after}\t{low_qual}\t{too_many_N}\t{too_short}\t{dup_rate}\t{insert_peak}\t{percent_lost}\n")

        # -----------------------------
        # Merged summary
        # -----------------------------
        with open(output.merged_summary, "w") as out_merged:
            out_merged.write("sample_id\tbefore_filtering_total_reads\tafter_filtering_total_reads\tlow_quality_reads\ttoo_many_N_reads\ttoo_short_reads\tduplication_rate\tinsert_peak\tpercent_lost\n")

            for sample in SAMPLES:
                report_file = f"step2_fastp/merged/reports/{sample}.fastp.json"
                with open(report_file) as f:
                    data = json.load(f)

                total_before = data["summary"]["before_filtering"]["total_reads"] / 2
                total_after  = data["summary"]["after_filtering"]["total_reads"]
                low_qual     = data["filtering_result"]["low_quality_reads"]
                too_many_N   = data["filtering_result"]["too_many_N_reads"]
                too_short    = data["filtering_result"]["too_short_reads"]
                dup_rate     = data["duplication"]["rate"]
                insert_peak  = data["insert_size"]["peak"]
                percent_lost = 100 * (total_before - total_after) / total_before if total_before > 0 else 0

                out_merged.write(f"{sample}\t{total_before}\t{total_after}\t{low_qual}\t{too_many_N}\t{too_short}\t{dup_rate}\t{insert_peak}\t{percent_lost}\n")

########################################
# Step 3: MetaPhlAn
# written for MetaPhlAn 4.2.2
# make sure to specify the index and check if any new databases are available 
########################################
rule metaphlan:
    input:
        r1="step2_fastp/pairedend/{sample}_R1_paired.fastq.gz",
        r2="step2_fastp/pairedend/{sample}_R2_paired.fastq.gz"
    output:
        bt="step3_metaphlan/bowtie_outs/{sample}.metaphlan.bz2",
        sams="step3_metaphlan/sams/{sample}.metaphlan.sam",
        profile="step3_metaphlan/sgb_profile/{sample}.sgb.metaphlan4-2-2.txt"
    conda:
        "/quobyte/dglemaygrp/nkumar/conda_envs/metaphlan.yaml"
    threads: 4
    resources:
           mem_mb=64*GB,
           runtime=12*60
    shell:
        """
        metaphlan {input.r1},{input.r2} --input_type fastq \
        --db_dir /quobyte/dglemaygrp/BACKED-UP/nkumar/database/ \
        --nproc 4 \
        --index mpa_vJan25_CHOCOPhlAnSGB_202503 \
        --mapout {output.bt} -s {output.sams} \
        -o {output.profile} --verbose --offline 
        """
########################################
# Step 3.1: GTDB taxonomy
# python script from MetaPhlAn github
# this script is specific to the index used in rule metaphlan,
# so there will probably be an updated script if a new index is released 
# wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/sgb_to_gtdb_profile.py -O scripts/sgb_to_gtdb_profile.py
# wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/util_fun.py -O scripts/util_fun.py
# wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/mpa_vJan25_CHOCOPhlAnSGB_202503_SGB2GTDB.tsv -O scripts/mpa_vJan25_CHOCOPhlAnSGB_202503_SGB2GTDB.tsv
########################################
rule gtdb:
    input: 
        profile="step3_metaphlan/sgb_profile/{sample}.sgb.metaphlan4-2-2.txt"
    output: 
        profile2="step3_metaphlan/gtdb_profile/{sample}.gtdb.metaphlan4-2-2.txt"
    conda:
        "/quobyte/dglemaygrp/nkumar/conda_envs/metaphlan.yaml"
    shell:
        """
        scripts/sgb_to_gtdb_profile.py -i {input.profile} -o {output.profile2}
        """
########################################
# Step 3.2: Combine metaphlan
# python script from MetaPhlAn github
# this script is specific to the index used in rule metaphlan,
# so there will probably be an updated script if a new index is released 
# wget https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/merge_metaphlan_tables.py -O scripts/merge_metaphlan_tables.py
########################################
rule combine_metaphlan_sgb:
    input:
        expand("step3_metaphlan/sgb_profile/{sample}.sgb.metaphlan4-2-2.txt", sample=SAMPLES)
    output:
        f"step3_metaphlan/{PROJECT_NAME}_sgb.metaphlan4-2-2.txt"
    conda:
        "/quobyte/dglemaygrp/nkumar/conda_envs/metaphlan.yaml"
    shell:
        """
        scripts/merge_metaphlan_tables.py {input} > {output}
        """