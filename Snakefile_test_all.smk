configfile: "config.yaml"

import os

# --- Функция для предобработки кластеров VAMB ---
def preprocess_vamb_clusters(cluster_in, cluster_out, fasta):
    """
    Удаляет заголовок из файла кластеров VAMB и фильтрует строки,
    оставляя только те контиги, которые реально присутствуют в FASTA.
    Это предотвращает KeyError в PHAMB и другие downstream-ошибки.
    """
    # Собираем имена всех контигов из FASTA
    fasta_contigs = set()
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                contig = line[1:].strip().split()[0]
                fasta_contigs.add(contig)
    # Фильтруем кластеры, удаляя заголовок и отсутствующие контиги
    with open(cluster_in) as fin, open(cluster_out, "w") as fout:
        next(fin)  # Пропускаем заголовок
        for line in fin:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            cluster, contig = parts[0], parts[1]
            if contig in fasta_contigs:
                fout.write(f"{cluster}\t{contig}\n")

# --- Главная цель пайплайна ---
rule all:
    input:
        os.path.join(config["results_dir"], "phamb/viral_bins")

# --- Создание директории для аннотаций ---
rule create_annotations_dir:
    output:
        directory(os.path.join(config["results_dir"], "annotations"))
    shell:
        "mkdir -p {output}"

# --- Обрезка и фильтрация ридов ---
rule fastp:
    input:
        r1 = os.path.join(config["raw_reads_dir"], "{sample}_R1.fastq"),
        r2 = os.path.join(config["raw_reads_dir"], "{sample}_R2.fastq")
    output:
        r1 = os.path.join(config["results_dir"], "fastp/{sample}_R1.trimmed.fq.gz"),
        r2 = os.path.join(config["results_dir"], "fastp/{sample}_R2.trimmed.fq.gz")
    conda:
        "../envs/fastp"
    log:
        "../workflow/logs/fastp_{sample}.log"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} 2> {log}"

# --- Сборка контигов ---
rule megahit:
    input:
        r1 = lambda wildcards: expand(
            os.path.join(config["results_dir"], "fastp/{sample}_R1.trimmed.fq.gz"),
            sample=config["samples"]
        ),
        r2 = lambda wildcards: expand(
            os.path.join(config["results_dir"], "fastp/{sample}_R2.trimmed.fq.gz"),
            sample=config["samples"]
        )
    output:
        contigs = os.path.join(config["results_dir"], "megahit/contigs.fa")
    conda:
        "../envs/megahit"
    log:
        "../workflow/logs/megahit.log"
    threads: config["threads"]
    shell:
        """
        rm -rf {config[results_dir]}/megahit_tmp
        megahit -1 {input.r1} -2 {input.r2} -o {config[results_dir]}/megahit_tmp --min-contig-len {config[contig_min_length]} -t {threads} 2> {log}
        mv {config[results_dir]}/megahit_tmp/final.contigs.fa {output.contigs}
        """

# --- Выравнивание ридов на сборку ---
rule strobealign:
    input:
        ref = rules.megahit.output.contigs,
        r1 = lambda wildcards: expand(
            os.path.join(config["results_dir"], "fastp/{sample}_R1.trimmed.fq.gz"),
            sample=config["samples"]
        ),
        r2 = lambda wildcards: expand(
            os.path.join(config["results_dir"], "fastp/{sample}_R2.trimmed.fq.gz"),
            sample=config["samples"]
        )
    output:
        bam = temp(os.path.join(config["results_dir"], "strobealign/aligned.bam"))
    conda:
        "../envs/strobealign"
    log:
        "../workflow/logs/strobealign.log"
    threads: config["threads"]
    shell:
        "strobealign -t {threads} {input.ref} {input.r1} {input.r2} | samtools sort -@ {threads} -o {output.bam} - 2> {log}"

# --- Бининг контигов с помощью VAMB ---
rule vamb:
    input:
        contigs = rules.megahit.output.contigs,
        bam = rules.strobealign.output.bam
    output:
        clusters = os.path.join(config["results_dir"], "vamb/vae_clusters_unsplit.tsv"),
        metadata = os.path.join(config["results_dir"], "vamb/vae_clusters_metadata.tsv")
    conda:
        "../envs/vamb"
    log:
        "../workflow/logs/vamb.log"
    shell:
        """
        rm -rf {config[results_dir]}/vamb
        vamb bin default --outdir {config[results_dir]}/vamb --fasta {input.contigs} --bamfiles {input.bam} 2> {log}
        """

# --- Фильтрация кластеров VAMB (удаление заголовка и несуществующих контигов) ---
rule preprocess_clusters:
    input:
        clusters = rules.vamb.output.clusters,
        contigs = rules.megahit.output.contigs
    output:
        filtered_clusters = os.path.join(config["results_dir"], "vamb/vae_clusters_unsplit.filtered.tsv")
    run:
        preprocess_vamb_clusters(
            input.clusters,
            output.filtered_clusters,
            input.contigs
        )

# --- Генерация белковых последовательностей с помощью Prodigal ---
rule prodigal:
    input:
        contigs = rules.megahit.output.contigs,
        dir = rules.create_annotations_dir.output
    output:
        proteins = os.path.join(config["results_dir"], "annotations/proteins.faa")
    conda:
        "../envs/prodigal"
    log:
        "../workflow/logs/prodigal.log"
    shell:
        "prodigal -i {input.contigs} -a {output.proteins} -p meta 2> {log}"

# --- Поиск вирусных последовательностей с помощью DeepVirFinder ---
rule deepvirfinder:
    input:
        contigs = rules.megahit.output.contigs,
        dir = rules.create_annotations_dir.output
    output:
        predictions = os.path.join(config["results_dir"], "annotations/dvf_predictions.tsv")
    conda:
        "../envs/deepvirfinder"
    log:
        "../workflow/logs/deepvirfinder.log"
    shell:
        r"""
        set -exuo pipefail
        cd ../DeepVirFinder
        ./dvf.py -i "$(realpath {input.contigs})" \
                 -o "$(realpath {config[results_dir]}/annotations)" \
                 -l 1000 2> "{log}"
        mv "$(realpath {config[results_dir]}/annotations)/$(basename {input.contigs})_gt1000bp_dvfpred.txt" "{output.predictions}"
        """

# --- Финальная агрегация вирусных бинов с помощью PHAMB ---
rule phamb:
    input:
        contigs = rules.megahit.output.contigs,
        clusters = rules.preprocess_clusters.output.filtered_clusters,
        annotations = rules.create_annotations_dir.output
    output:
        viral_bins = directory(os.path.join(config["results_dir"], "phamb/viral_bins"))
    conda:
        "../envs/phamb"
    log:
        "../workflow/logs/phamb.log"
    params:
        script = config["phamb_script"]
    shell:
        "python {params.script} {input.contigs} {input.clusters} {input.annotations} {output.viral_bins} 2> {log}"
