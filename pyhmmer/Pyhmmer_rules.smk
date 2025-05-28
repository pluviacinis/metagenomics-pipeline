rule pyhmmer_vog:
    input:
        proteins = rules.prodigal.output.proteins,
        dir_marker = rules.create_annotations_dir.output
    output:
        annotations = os.path.join(config["results_dir"], "annotations/all.hmmVOG.tbl")
    params: db = config["vog_db"] + ".h3m"
    conda: "../envs/pyhmmer"
    log: "../workflow/logs/pyhmmer_vog.log"
    threads: 4
    shell:
        "./pyhmmer_scan.py {params.db} {input.proteins} "
        "-o {output.annotations} "
        "--cpu {threads} "
        "-E 1e-5 "
        "2> {log}"

rule pyhmmer_micomplete:
    input:
        proteins = rules.prodigal.output.proteins,
        dir_marker = rules.create_annotations_dir.output
    output:
        annotations = os.path.join(config["results_dir"], "annotations/all.hmmMiComplete105.tbl")
    params: db = config["miComplete_db"] + ".h3m"
    conda: "../envs/pyhmmer"
    log: "../workflow/logs/pyhmmer_micomplete.log"
    threads: 4
    shell:
        "./pyhmmer_scan.py {params.db} {input.proteins} "
        "-o {output.annotations} "
        "--cpu {threads} "
        "-E 1e-5 "
        "2> {log}"