configfile: "bsa-pipeline-config.yaml"

rule all:
    input:
        [f"results/variant_to_table/{group}.table" for group in config["samples"]]

rule fastp:
    input:
        R1 = lambda w: config["samples"][w.group][w.bulk]["R1"],
        R2 = lambda w: config["samples"][w.group][w.bulk]["R2"]
    output:
        R1_fastp = "resources/fastp/{group}_{bulk}_R1_fastp.fastq.gz",
        R2_fastp = "resources/fastp/{group}_{bulk}_R2_fastp.fastq.gz",
        json = "results/fastp-reports/{group}_{bulk}_fastp.json",
        html = "results/fastp-reports/{group}_{bulk}_fastp.html"
    log:
        "logs/fastp/{group}_{bulk}.log"
    shell:
        '''
        fastp \
        -i {input.R1} \
        -I {input.R2} \
        -o {output.R1_fastp} \
        -O {output.R2_fastp} \
        -j {output.json} \
        -h {output.html} \
        -R {wildcards.group}_{wildcards.bulk}
        '''

rule make_bwa_index:
    input:
        ref_genome_fa = config["ref_genome"]
    output:
        bwa_index = expand("resources/bwa_index.{suffix}", suffix = ["amb", "ann", "pac", "bwt", "sa"])
    log:
        "logs/make_bwa_index.log"
    resources:
        mem_gb = 12
    shell:
        "bwa index -p resources/bwa_index {input.ref_genome_fa}"

rule bwa_alignment:
    input:
        R1_fastp = "resources/fastp/{group}_{bulk}_R1_fastp.fastq.gz",
        R2_fastp = "resources/fastp/{group}_{bulk}_R2_fastp.fastq.gz",
        bwa_index = expand("resources/bwa_index.{suffix}", suffix = ["amb", "ann", "pac", "bwt", "sa"])
    output:
        sorted_bam = "resources/bwa_alignment/{group}_{bulk}_sorted.bam"
    log:
        "logs/bwa_aligment/{group}_{bulk}.log"
    params:
        RG_ID = "{group}-{bulk}"
    resources:
        core_num = 12,
        mem_gb = 12
    shell:
        '''
        bwa mem -M -t {resources.core_num} \
            -R '@RG\\tID:{params.RG_ID}\\tSM:{params.RG_ID}' \
            resources/bwa_index \
            {input.R1_fastp} \
            {input.R2_fastp} | \
        samtools sort \
            -o {output.sorted_bam}
        '''

rule filter:  # Remove unmapped and duplicated reads
    input:
        sorted_bam = "resources/bwa_alignment/{group}_{bulk}_sorted.bam"
    output:
        filtered_bam = "resources/filter/{group}_{bulk}_sorted_filtered.bam",
        dup_metrics = "resources/filter/{group}_{bulk}_dup_metrics.txt",
        filtered_index = "resources/filter/{group}_{bulk}_sorted_filtered.bai",
        temp_bam = temp("temp/{group}_{bulk}.bam")
    log:
        "logs/filter/{group}_{bulk}.log"
    resources:
        mem_gb = 12
    shell:  # MarkDuplicates can NOT read input from STDIN
        '''
        samtools view -b -h -F 4 {input.sorted_bam} > {output.temp_bam}
        gatk --java-options "-Xmx{resources.mem_gb}G"  MarkDuplicates \
            --REMOVE_DUPLICATES true \
            --INPUT {output.temp_bam} \
            --OUTPUT {output.filtered_bam} \
            --METRICS_FILE {output.dup_metrics} \
            --CREATE_INDEX true
        '''

rule make_gatk_fa_index:  # Make fai and dict index file for ref genome FASTA, both files are required by GATK pipeline
    input:  # fai index also required by VEP
        ref_genome_fa = config["ref_genome"],
    output:
        ref_genome_fa_copy = "resources/ref_genome.fa",
        ref_genome_fa_index = "resources/ref_genome.fa.fai",
        ref_genome_fa_dict = "resources/ref_genome.dict"
    log:
        "logs/make_gatk_fa_index.log"
    run:
        if input.ref_genome_fa[-3:] == ".gz":
            shell("zcat < {input.ref_genome_fa} > {output.ref_genome_fa_copy}")
        else:
            shell("cp {input.ref_genome_fa} {output.ref_genome_fa_copy}")
        shell("samtools faidx {output.ref_genome_fa_copy}")
        shell("gatk CreateSequenceDictionary -R {output.ref_genome_fa_copy}")

rule haplotype_caller:
    input:
        ref_genome_fa_copy = "resources/ref_genome.fa",
        ref_genome_fa_index = "resources/ref_genome.fa.fai",
        ref_genome_fa_dict = "resources/ref_genome.dict",
        filtered_bam = "resources/filter/{group}_{bulk}_sorted_filtered.bam",
        filtered_index = "resources/filter/{group}_{bulk}_sorted_filtered.bai",
    output:
        gvcf = "resources/haplotype_caller/{group}_{bulk}.g.vcf.gz",
        gvcf_index = "resources/haplotype_caller/{group}_{bulk}.g.vcf.gz.tbi"
    log:
        "logs/haplotype_caller/{group}_{bulk}.log"
    resources:
        mem_gb = 12
    shell:
        '''
        gatk --java-options "-Xmx{resources.mem_gb}G" HaplotypeCaller \
            -R {input.ref_genome_fa_copy} \
            -I {input.filtered_bam} \
            -ERC GVCF \
            -O {output.gvcf}
        '''

rule combine_gvcfs:
    input:
        gvcf = lambda w: [f"resources/haplotype_caller/{w.group}_{bulk}.g.vcf.gz" for bulk in config["samples"][w.group]],
        ref_genome_fa_copy = "resources/ref_genome.fa",
        ref_genome_fa_index = "resources/ref_genome.fa.fai",
        ref_genome_fa_dict = "resources/ref_genome.dict"
    output:
        combined_gvcfs = "resources/combine_gvcfs/{group}_combined.g.vcf.gz",
        combined_gvcfs_index = "resources/combine_gvcfs/{group}_combined.g.vcf.gz.tbi"
    log:
        "logs/combine_gvcfs/{group}.log"
    resources:
        mem_gb = 4
    run:
        command_fragments = [
            "gatk --java-options '-Xmx{resources.mem_gb}G' CombineGVCFs",
            "-R {input.ref_genome_fa_copy}",
            "-O {output.combined_gvcfs}"
        ]
        command_fragments.extend([f"--variant {gvcf}" for gvcf in input.gvcf])
        shell(" ".join(command_fragments))

rule genotype_gvcfs:
    input:
        combined_gvcfs = "resources/combine_gvcfs/{group}_combined.g.vcf.gz",
        ref_genome_fa_copy = "resources/ref_genome.fa",
        ref_genome_fa_index = "resources/ref_genome.fa.fai",
        ref_genome_fa_dict = "resources/ref_genome.dict"
    output:
        genotyped_vcf = "resources/genotype_gvcfs/{group}_combined_genotyped.vcf.gz",
        genotyped_vcf_index = "resources/genotype_gvcfs/{group}_combined_genotyped.vcf.gz.tbi"
    log:
        "logs/genotype_gvcfs/{group}.log"
    resources:
        mem_gb = 4
    shell:
        '''
        gatk --java-options "-Xmx{resources.mem_gb}G" GenotypeGVCFs \
            -R {input.ref_genome_fa_copy} \
            -V {input.combined_gvcfs} \
            -O {output.genotyped_vcf}
        '''

rule tabix_gff:
    input:
        anno_gff = config["anno_gff"]
    output:
        bgzip_gff= "resources/annotation.gff.gz",
        bgzip_gff_index = "resources/annotation.gff.gz.tbi"
    log:
        "logs/tabix_gff.log"
    run:
        pipeline = [
            "cat < {input.anno_gff}",
            "grep -v '^#'",
            "sort -k1,1 -k4,4n -k5,5n -t$'\\t'",  # sort gff (from vep document)
            "bgzip -c > {output.bgzip_gff}"
        ]
        if input.anno_gff[-3:] == ".gz":
            pipeline[0] = "zcat < {input.anno_gff}"
        shell(" | ".join(pipeline))
        shell("tabix -p gff {output.bgzip_gff}")

rule vep_annotation:
    input:
        genotyped_vcf = "resources/genotype_gvcfs/{group}_combined_genotyped.vcf.gz",
        ref_genome_fa_copy = "resources/ref_genome.fa",
        ref_genome_fa_index = "resources/ref_genome.fa.fai",
        bgzip_gff= "resources/annotation.gff.gz",
        bgzip_gff_index = "resources/annotation.gff.gz.tbi"
    output:
        annotated_vcf = "resources/vep_annotation/{group}_combined_genotyped_annotated.vcf.gz",
        annotated_vcf_index = "resources/vep_annotation/{group}_combined_genotyped_annotated.vcf.gz.tbi",
        annotation_summary = "results/vep_annotation/{group}_annotation_summary.html"
    log:
        "logs/vep_annotation/{group}.log"
    shell:
        '''
        vep --vcf \
            --input_file {input.genotyped_vcf} \
            --gff {input.bgzip_gff} \
            --fasta {input.ref_genome_fa_copy} \
            --compress_output bgzip \
            --fields "Consequence,Gene,Amino_acids" \
            --output_file {output.annotated_vcf} \
            --stats_file {output.annotation_summary}
        tabix -p vcf {output.annotated_vcf}
        '''

rule variant_to_table:
    input:
        annotated_vcf = "resources/vep_annotation/{group}_combined_genotyped_annotated.vcf.gz",
        ref_genome_fa_copy = "resources/ref_genome.fa",
        ref_genome_fa_index = "resources/ref_genome.fa.fai",
        ref_genome_fa_dict = "resources/ref_genome.dict"
    output:
        table = "results/variant_to_table/{group}.table"
    log:
        "logs/variant_to_table/{group}.log"
    shell:
        '''
        gatk VariantsToTable \
            -R {input.ref_genome_fa_copy} \
            -V {input.annotated_vcf} \
            -F CHROM -F POS -F REF -F ALT \
            -GF AD -GF DP -GF GQ -GF PL \
            -O {output.table}
        '''
