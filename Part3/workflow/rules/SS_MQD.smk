# NOTE: Laser-capture microdissection analysis ---------------------------------------
# Preprocess fastq with fastp
rule preprocess_paired_reads_with_fastp:
    input:
        fastq_read1 = lambda wildcards: config["MQD_FASTQ"][wildcards.sample][0].replace("RAWDATA_PATH", config["RAWDATA_PATH"]),
        fastq_read2 = lambda wildcards: config["MQD_FASTQ"][wildcards.sample][1].replace("RAWDATA_PATH", config["RAWDATA_PATH"])
    output:
        fastq_read1 = "results/MQD_analysis/fastq_preprocessed/{sample}_afterQC_1.fq.gz",
        fastq_read2 = "results/MQD_analysis/fastq_preprocessed/{sample}_afterQC_2.fq.gz"
    log:
        "logs/" + DATE + "_{sample}_fastp.log"
    threads: 20
    shell:
        """
        bash workflow/scripts/fastp.sh \
            -t {threads} \
            -a {input.fastq_read1} \
            -b {input.fastq_read2} \
            -x {output.fastq_read1} \
            -y {output.fastq_read2} \
            -l {log}
        """

rule preprocess_single_reads_with_fastp:
    input:
        fastq_read = lambda wildcards: config["MQD_FASTQ"][wildcards.sample][0].replace("RAWDATA_PATH", config["RAWDATA_PATH"])
    output:
        fastq_read = "results/MQD_analysis/fastq_preprocessed/{sample}_afterQC.fq.gz"
    log:
        "logs/" + DATE + "_{sample}_fastp.log"
    threads: 20
    shell:
        """
        bash workflow/scripts/fastp.sh \
            -t {threads} \
            -c {input.fastq_read} \
            -z {output.fastq_read} \
            -l {log}
        """

# Build index file for read mapping with hisat2
rule build_mapping_index_for_hisat2:
    input:
        genome_fasta = lambda wildcards: config["REF_GENOME_FASTA"][wildcards.species],
        genome_gtf = lambda wildcards: config["REF_GENOME_GTF"][wildcards.species]
    output:
        ref_index_1 = "results/MQD_analysis/mapping_index/{species}_index.1.ht2"
    params:
        ref_ss = "results/MQD_analysis/mapping_index/{species}_ss",
        ref_exon = "results/MQD_analysis/mapping_index/{species}_exon",
        ref_index = "results/MQD_analysis/mapping_index/{species}_index"
    threads:
        workflow.cores
    shell:
        """
        echo '=====Build the index file for reads alignment with HISAT2====='
        hisat2_extract_splice_sites.py {input.genome_gtf} > {params.ref_ss}
        hisat2_extract_exons.py {input.genome_gtf} > {params.ref_exon}

        hisat2-build \
            -f \
            -p {threads} \
            --exon {params.ref_exon} \
            --ss {params.ref_ss} \
            {input.genome_fasta} \
            {params.ref_index}
        """

# Map reads and turn sam into sorted bam
rule map_paired_reads_to_genome_with_hisat2:
    input:
        ref_index_1 = lambda wildcards: "results/MQD_analysis/mapping_index/" + config["REF_SPECIES"][wildcards.sample] + "_index.1.ht2",
        fastq_read1 = "results/MQD_analysis/fastq_preprocessed/{sample}_afterQC_1.fq.gz",
        fastq_read2 = "results/MQD_analysis/fastq_preprocessed/{sample}_afterQC_2.fq.gz"
    output:
        sorted_bam = "results/MQD_analysis/sorted_bam_paired/{sample}.sorted.bam"
    params:
        ref_index = lambda wildcards: "results/MQD_analysis/mapping_index/" + config["REF_SPECIES"][wildcards.sample] + "_index",
        max_intron_length = lambda wildcards: config["REF_MAX_INTRON_LENGTH"][config["REF_SPECIES"][wildcards.sample]]
    log:
        "logs/" + DATE + "_{sample}_summary.log"
    threads:
        workflow.cores
    shell:
        """
        echo '=====Align the reads with HISAT2====='
        date +%Y-%m-%d_%H:%M:%S
        hisat2 \
            --threads {threads} \
            -q --phred33 \
            --max-intronlen {params.max_intron_length} \
            --secondary \
            --fr \
            --rna-strandness FR \
            --new-summary \
            --summary-file {log} \
            -x {params.ref_index} \
            -1 {input.fastq_read1} \
            -2 {input.fastq_read2} \
            -S results/MQD_analysis/sorted_bam_paired/{wildcards.sample}.sam \
            --dta
        # --rna-strandness FR is for stranded library
        # --dta is for StringTie

        echo '=====Convert the sam file into bam file (sam -> bam -> sorted.bam)====='
        date +%Y-%m-%d_%H:%M:%S
        samtools view -b --threads {threads} \
        results/MQD_analysis/sorted_bam_paired/{wildcards.sample}.sam |
        samtools sort -m 2G --threads {threads} \
        -o results/MQD_analysis/sorted_bam_paired/{wildcards.sample}.sorted.bam
        rm results/MQD_analysis/sorted_bam_paired/{wildcards.sample}.sam
        
        echo '=====Build the bam index file [sorted.bam +> bai]====='
        date +%Y-%m-%d_%H:%M:%S
        samtools index -@ {threads} \
        results/MQD_analysis/sorted_bam_paired/{wildcards.sample}.sorted.bam \
        results/MQD_analysis/sorted_bam_paired/{wildcards.sample}.sorted.bam.bai
        """

rule map_single_reads_to_genome_with_hisat2:
    input:
        ref_index_1 = lambda wildcards: "results/MQD_analysis/mapping_index/" + config["REF_SPECIES"][wildcards.sample] + "_index.1.ht2",
        fastq_read = "results/MQD_analysis/fastq_preprocessed/{sample}_afterQC.fq.gz",
    output:
        sorted_bam = "results/MQD_analysis/sorted_bam_single/{sample}.sorted.bam"
    params:
        ref_index = lambda wildcards: "results/MQD_analysis/mapping_index/" + config["REF_SPECIES"][wildcards.sample] + "_index",
        max_intron_length = lambda wildcards: config["REF_MAX_INTRON_LENGTH"][config["REF_SPECIES"][wildcards.sample]]
    log:
        "logs/" + DATE + "_{sample}_summary.log"
    threads:
        workflow.cores
    shell:
        """
        echo '=====Align the reads with HISAT2====='
        date +%Y-%m-%d_%H:%M:%S
        hisat2 \
            --threads {threads} \
            -q --phred33 \
            --max-intronlen {params.max_intron_length} \
            --secondary \
            --fr \
            --new-summary \
            --summary-file {log} \
            -x {params.ref_index} \
            -U {input.fastq_read} \
            -S results/MQD_analysis/sorted_bam_single/{wildcards.sample}.sam \
            --dta
        # --rna-strandness FR is for stranded library
        # --dta is for StringTie

        echo '=====Convert the sam file into bam file (sam -> bam -> sorted.bam)====='
        date +%Y-%m-%d_%H:%M:%S
        samtools view -b --threads {threads} \
        results/MQD_analysis/sorted_bam_single/{wildcards.sample}.sam |
        samtools sort -m 2G --threads {threads} \
        -o results/MQD_analysis/sorted_bam_single/{wildcards.sample}.sorted.bam
        rm results/MQD_analysis/sorted_bam_single/{wildcards.sample}.sam
        
        echo '=====Build the bam index file [sorted.bam +> bai]====='
        date +%Y-%m-%d_%H:%M:%S
        samtools index -@ {threads} \
        results/MQD_analysis/sorted_bam_single/{wildcards.sample}.sorted.bam \
        results/MQD_analysis/sorted_bam_single/{wildcards.sample}.sorted.bam.bai
        """

# Calcualte transcript abundance with stringtie
rule quantify_with_stringtie:
    input:
        genome_gtf = lambda wildcards: config["REF_GENOME_GTF"][config["REF_SPECIES"][wildcards.sample]],
        sorted_bam = lambda wildcards: "results/MQD_analysis/sorted_bam_single/{sample}.sorted.bam" if len(config["MQD_FASTQ"][wildcards.sample]) == 1 else "results/MQD_analysis/sorted_bam_paired/{sample}.sorted.bam"
    output:
        assembled_transcripts_gtf = "results/MQD_analysis/quantification/{sample}_assembled_transcripts.gtf",
        gene_abundance_tsv = "results/MQD_analysis/quantification/{sample}_gene_abundances.tsv"
    log:
        "logs/" + DATE + "_{sample}_stringtie.log"
    threads:
        workflow.cores
    shell:
        """
        echo '=====Calculate the read counts of each transcript with StringTie====='
        stringtie {input.sorted_bam} \
            --fr \
            -G {input.genome_gtf} \
            -o {output.assembled_transcripts_gtf} \
            -e -A {output.gene_abundance_tsv} \
            -p {threads} \
            -v > {log}
        """

# Turn stringtie quantification results into csv
rule turn_stringtie_results_into_csv:
    input:
        gene_abundance_tsv_list = lambda wildcards: expand(
            "results/MQD_analysis/quantification/{sample}_gene_abundances.tsv",
            sample = config["MQD_SAMPLE"][wildcards.group]
        )
    output:
        gene_abundance_csv = "results/MQD_analysis/quantification/{group}_gene_abundances.csv"
    script:
        "../scripts/turn_stringtie_results_into_csv.R"

# Differential expresssion analysis
def get_all_required_samples(DEAplan: str) -> list:
    all_groups = config["DEA_RULE"][DEAplan]["EXP_GROUP"] + config["DEA_RULE"][DEAplan]["CTRL_GROUP"]
    all_samples = []
    for group in all_groups:
        all_samples += config["MQD_SAMPLE"][group]
    return all_samples

rule turn_stringtie_gtf_into_read_counts:
    input:
        assembled_transcripts_gtf = "results/MQD_analysis/quantification/{sample}_assembled_transcripts.gtf"
    output:
        gene_count_matrix_csv = "results/MQD_analysis/quantification/hypothetical_read_counts/{sample}_gene_count_matrix.csv",
        transcript_count_matrix_csv = "results/MQD_analysis/quantification/hypothetical_read_counts/{sample}_transcript_count_matrix.csv"
    params:
        read_length = lambda wildcards: config["SEQ_LENGTH"][wildcards.sample]
    shell:
        """
        echo '=====Extract raw counts from the StringTie-output GTFs====='
        mkdir --parents temp/{wildcards.sample}/{wildcards.sample}
        cp {input.assembled_transcripts_gtf} temp/{wildcards.sample}/{wildcards.sample}
        prepDE.py3 -i temp/{wildcards.sample} \
                -l {params.read_length} \
                -g {output.gene_count_matrix_csv} \
                -t {output.transcript_count_matrix_csv}
        rm -r temp/{wildcards.sample}
        """

rule differential_expression_analysis_with_DEseq2:
    input:
        gene_count_matrix_csvs = lambda wildcards: expand(
            "results/MQD_analysis/quantification/hypothetical_read_counts/{sample}_gene_count_matrix.csv",
            sample = get_all_required_samples(wildcards.DEAplan)
        )
    output:
        DEA_folder = directory("results/MQD_analysis/differential_expression_analysis/{DEAplan}")
    params:
        DEAplan = lambda wildcards: config["DEA_RULE"][wildcards.DEAplan],
        Group2Sample = lambda wildcards: config["MQD_SAMPLE"]
    script:
        "../scripts/do_DEA_with_DESeq2.R"
    # shell:
    #     """
    #     echo '=====Conduct diffetential epxression analysis with DESeq2====='
    #     date +%Y-%m-%d_%H:%M:%S
    #     Rscript process_DESeq2.R $inputFileList outputStringtie_gene_count_matrix.csv $isPaired $sigSide $FDR $foldChange
    #     """
