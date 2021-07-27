SEMIBIN_DATA_PATH = os.path.join(DBDIR, "SemiBin_GTDB")

localrules: semibin_download_gtdb
rule semibin_download_gtdb:
    output:
        directory(SEMIBIN_DATA_PATH),
    log:
        "log/download/Semibin.txt",
    conda:
        "../envs/semibin.yaml"
    threads: 1
    shell:
        "SemiBin download_GTDB --reference-db {output}/GTDB 2> {log}"

        # Semibin 0.2 has the following error https://github.com/BigDataBiology/SemiBin/issues/31


rule semibin_predict_taxonomy:
    input:
        fasta=rules.combine_contigs.output,
        db=SEMIBIN_DATA_PATH,
    output:
        "Crossbinning/SemiBin/inconsitent_taxonomy.tsv",
    conda:
        "../envs/semibin.yaml"
    threads: 1
    resources:
        mem=config["large_mem"],
        time=config["runtime"]["default"],
    log:
        "log/semibin/predict_taxonomy.log",
    benchmark:
        "log/benchmarks/semibin/predict_taxonomy.tsv"
    params:
        output_dir=lambda wc, output: os.path.dirname(output[0]),
        name=lambda wc, output: os.path.basename(output[0]),
    shell:
        "SemiBin predict_taxonomy "
        " --input-fasta {input.fasta} "
        " --output {params.output_dir} "
        " --cannot-name {params.name} "
        " --reference-db {input.db}/GTDB "
        " > {log} 2> {log}"


rule semibin_generate_data_multi:
    input:
        fasta=rules.combine_contigs.output,
        bams=expand(rules.sort_bam.output, sample=SAMPLES),
    output:
        expand("Crossbinning/SemiBin/samples/{sample}/{files}",
        sample=SAMPLES,
        files=  ["data.csv","data_split.csv"]
        )
    conda:
        "../envs/semibin.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
        time=config["runtime"]["default"],
    log:
        "log/semibin/generate_data_multi.log",
    benchmark:
        "log/benchmarks/semibin/generate_data_multi.tsv"
    params:
        output_dir="Crossbinning/SemiBin",
        separator="_",
    shell:
        "SemiBin generate_data_multi "
        " --input-fasta {input.fasta} "
        " --input-bam {input.bams} "
        " --output {params.output_dir} "
        " --threads {threads} "
        " --separator {params.separator} "
        " 2> {log}"


rule semibin_train:
    input:
        fasta=rules.combine_contigs.output,
        bams=expand(rules.sort_bam.output, sample=SAMPLES),
        data="Crossbinning/SemiBin/samples/{sample}/data.csv",
        data_split="Crossbinning/SemiBin/samples/{sample}/data_split.csv",
        cannot_link=rules.semibin_predict_taxonomy.output[0],
    output:
        "Crossbinning/SemiBin/samples/{sample}/model.h5",
    conda:
        "../envs/semibin.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
        time=config["runtime"]["default"],
    log:
        "log/semibin/train/{sample}.log",
    benchmark:
        "log/benchmarks/semibin/train/{sample}.tsv"
    params:
        output_dir=lambda wc, output: os.path.dirname(output[0]),
        extra=" --epoches 20",
    shell:
        "SemiBin train "
        " --input-fasta {input.fasta} "
        " --input-bam {input.bams} "
        " --output {params.output_dir} "
        " --threads {threads} "
        " --data {input.data} "
        " --data_split {input.data_split} "
        " --cannot-link {input.cannot_link} "
        " {params.extra} "
        " 2> {log}"


rule run_semibin:
    input:
        fasta=rules.combine_contigs.output,
        bams=expand(rules.sort_bam.output, sample=SAMPLES),
        data="Crossbinning/SemiBin/samples/{sample}/data.csv",
        model=rules.semibin_train.output[0],
    output:
        touch("Crossbinning/SemiBin/samples/{sample}/finished"),
    conda:
        "../envs/semibin.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
        time=config["runtime"]["default"],
    log:
        "log/semibin/bin/{sample}.log",
    benchmark:
        "log/benchmarks/semibin/bin/{sample}.tsv"
    params:
        output_dir="Crossbinning/SemiBin",
        extra=" --minfasta-kbs 200 --recluster --max-node 1 --max-edges 200 ",
    shell:
        "SemiBin train "
        " --input-fasta {input.fasta} "
        " --input-bam {input.bams} "
        " --output {params.output_dir} "
        " --threads {threads} "
        " --data {input.data} "
        " --model {input.model} "
        " {params.extra} "
        " 2> {log}"

rule semibin:
    input:
        expand("Crossbinning/SemiBin/samples/{sample}/finished", sample=SAMPLES)
# alternative to pretrained model --environment: Environment for the built-in model(human_gut/dog_gut/ocean).”