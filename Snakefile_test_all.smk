configfile: "config.yaml"
import os

def preprocess_vamb_clusters(cluster_in, cluster_out, fasta):
    fasta_contigs = set()
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                contig = line[1:].strip().split()[0]
                fasta_contigs.add(contig)
    with open(cluster_in) as fin, open(cluster_out, "w") as fout:
        next(fin)  # Пропускаем заголовок
        for line in fin:
            parts = line.strip().split("\t")
            if len(parts) < 2: continue
            cluster, contig = parts[0], parts[1]
            if contig in fasta_contigs:
                fout.write(f"{cluster}\t{contig}\n")

rule all:
    input:
        os.path.join(config["results_dir"], "phamb/viral_bins")

rule create_annotations_dir:
    output:
        touch(os.path.join(config["results_dir"], "annotations/.dir_created"))
    shell:
        "mkdir -p $(dirname {output}) && touch {output}"

rule fastp:
    input:
        r1 = os.path.join(config["raw_reads_dir"], "{sample}_R1.fastq"),
        r2 = os.path.join(config["raw_reads_dir"], "{sample}_R2.fastq")
    output:
        r1 = os.path.join(config["results_dir"], "fastp/{sample}_R1.trimmed.fq.gz"),
        r2 = os.path.join(config["results_dir"], "fastp/{sample}_R2.trimmed.fq.gz")
    conda: "../envs/fastp"
    log: "../workflow/logs/fastp_{sample}.log"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} 2> {log}"

rule megahit:
    input:
        r1 = expand(os.path.join(config["results_dir"], "fastp/{sample}_R1.trimmed.fq.gz"), sample=config["samples"]),
        r2 = expand(os.path.join(config["results_dir"], "fastp/{sample}_R2.trimmed.fq.gz"), sample=config["samples"])
    output:
        contigs = os.path.join(config["results_dir"], "megahit/contigs.fa")
    conda: "../envs/megahit"
    log: "../workflow/logs/megahit.log"
    threads: config["threads"]
    shell:
        """
        rm -rf {config[results_dir]}/megahit_tmp
        megahit -1 {input.r1} -2 {input.r2} -o {config[results_dir]}/megahit_tmp --min-contig-len {config[contig_min_length]} -t {threads} 2> {log}
        mv {config[results_dir]}/megahit_tmp/final.contigs.fa {output.contigs}
        """

rule strobealign:
    input:
        ref = rules.megahit.output.contigs,
        r1 = expand(os.path.join(config["results_dir"], "fastp/{sample}_R1.trimmed.fq.gz"), sample=config["samples"]),
        r2 = expand(os.path.join(config["results_dir"], "fastp/{sample}_R2.trimmed.fq.gz"), sample=config["samples"])
    output:
        bam = os.path.join(config["results_dir"], "strobealign/aligned.bam")
    conda: "../envs/strobealign"
    log: "../workflow/logs/strobealign.log"
    threads: config["threads"]
    shell:
        "strobealign -t {threads} {input.ref} {input.r1} {input.r2} | samtools sort -@ {threads} -o {output.bam} - 2> {log}"

rule normalize_headers:
    input: rules.megahit.output.contigs
    output: contigs = os.path.join(config["results_dir"], "megahit/contigs_clean.fa")
    shell: "sed -E 's/^>([^ ]+).*/>\\1/' {input} > {output.contigs}"

rule vamb:
    input:
        contigs = rules.megahit.output.contigs,
        bam = rules.strobealign.output.bam
    output:
        clusters = os.path.join(config["results_dir"], "vamb/vae_clusters_unsplit.tsv"),
        metadata = os.path.join(config["results_dir"], "vamb/vae_clusters_metadata.tsv")
    conda: "./envs/vamb"
    log: "../workflow/logs/vamb.log"
    shell:
        """
        rm -rf {config[results_dir]}/vamb
        vamb bin default --outdir {config[results_dir]}/vamb --fasta {input.contigs} --bamfiles {input.bam} 2> {log}
        """

rule preprocess_clusters:
    input:
        clusters = rules.vamb.output.clusters,
        contigs = rules.megahit.output.contigs
    output:
        filtered_clusters = os.path.join(config["results_dir"], "vamb/vae_clusters_unsplit.filtered.tsv")
    run: preprocess_vamb_clusters(input.clusters, output.filtered_clusters, input.contigs)

rule prodigal:
    input:
        contigs = rules.megahit.output.contigs,
        dir_marker = rules.create_annotations_dir.output
    output:
        proteins = os.path.join(config["results_dir"], "annotations/proteins.faa")
    conda: "../envs/prodigal"
    log: "../workflow/logs/prodigal.log"
    shell:
        "prodigal -i {input.contigs} -a {output.proteins} -p meta 2> {log}"

rule hmmer_vog:
    input:
        proteins = rules.prodigal.output.proteins,
        dir_marker = rules.create_annotations_dir.output
    output:
        annotations = os.path.join(config["results_dir"], "annotations/all.hmmVOG.tbl")
    params: db = config["vog_db"]
    conda: "../envs/hmmer"
    log: "../workflow/logs/hmmer_vog.log"
    threads: 4
    shell:
        "hmmscan --cpu {threads} --tblout {output.annotations} {params.db} {input.proteins} 2> {log}"

rule hmmer_micomplete:
    input:
        proteins = rules.prodigal.output.proteins,
        dir_marker = rules.create_annotations_dir.output
    output:
        annotations = os.path.join(config["results_dir"], "annotations/all.hmmMiComplete105.tbl")
    params: db = config["miComplete_db"]
    conda: "../envs/hmmer"
    log: "../workflow/logs/hmmer_micomplete.log"
    threads: 4
    shell:
        "hmmscan --cpu {threads} --tblout {output.annotations} {params.db} {input.proteins} 2> {log}"

rule deepvirfinder:
    input:
        contigs = rules.normalize_headers.output.contigs,
        dir_marker = rules.create_annotations_dir.output
    output:
        predictions = os.path.join(config["results_dir"], "annotations/all.DVF.predictions.txt")
    conda: "../envs/deepvirfinder"
    log: "../workflow/logs/deepvirfinder.log"
    shell:
        r"""
        cd ../DeepVirFinder
        ./dvf.py -i "{input.contigs}" -o "$(dirname {output.predictions})" -l 1000 2> {log}
        # Сохраняем оригинальный формат DVF (4 столбца)
        mv "$(dirname {output.predictions})/$(basename {input.contigs})_gt1000bp_dvfpred.txt" {output.predictions}
        """

rule phamb:
    input:
        contigs = rules.normalize_headers.output.contigs,
        clusters = rules.preprocess_clusters.output.filtered_clusters,
        vog = rules.hmmer_vog.output.annotations,
        micomplete = rules.hmmer_micomplete.output.annotations,
        dvf = rules.deepvirfinder.output.predictions
    output:
        viral_bins = directory(os.path.join(config["results_dir"], "phamb/viral_bins"))
    conda: "../envs/phamb"
    log: "../workflow/logs/phamb.log"
    params:
        script = config["phamb_script"],
        annotations_path = os.path.join(config["results_dir"], "annotations")
    shell:
        "python {params.script} {input.contigs} {input.clusters} {params.annotations_path} {output.viral_bins} 2> {log}"
