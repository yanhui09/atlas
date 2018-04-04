import os
import re
import sys
import tempfile
import warnings

from default_values import *

combined_contigs_folder='contigs'

#### combine contigs
config['perform_genome_binning']= False
MAGs=['metagenome']


rule combine_contigs_report:
    input:
        combined_contigs= COMBINED_CONTIGS,
        combined_contigs_stats="contigs/combined_contigs_stats.txt",
        median_coverage="contigs/combined_median_coverage.tsv",
        gc_stats = "contigs/combined_contigs_stats_gc.tsv",
        binned_coverage = "contigs/combined_coverage_binned.tsv.gz",
        gene_counts= expand('annotations/{MAG}/gene_counts.tsv',MAG=MAGs),
        eggNOG= expand('annotations/{MAG}/eggNOG_annotation.tsv',MAG=MAGs),
        taxonomy = expand("annotations/{MAG}/refseq/tax_assignments.tsv",MAG=MAGs),
        gene_info = expand("annotations/{MAG}/feature_counts/gene_info.tsv",MAG=MAGs),
        # genes = expand("annotations/{MAG}/predicted_genes/genes_plus.tsv",MAG=MAGs), # error in prokka
        # instead of annotations= expand("annotations/{MAG}/annotations.txt",MAG=MAGs)
        # concoct="{folder}/binning/{file}".format(folder=combined_contigs_folder,file='means_gt2500.csv')
    output:
        touch("Combined_contigs_done")



rule combine_contigs:
    input:
        expand("{sample}/assembly/{sample}_prefilter_contigs.fasta",sample=SAMPLES)
    output:
        combined_contigs=temp("{folder}/combined_contigs_oldnames.fasta"),
        cluster_stats="{folder}/combined_contigs_kmerfreq.txt",
        dot="{folder}/combined_contigs_graph.dot"
    benchmark:
        "logs/benchmarks/combine_contigs.txt"
    log:
        "logs/combine_contigs.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("java_mem", JAVA_MEM)
    params:
       input=lambda wc,input: ','.join(input),
       min_length=config.get("minimum_contig_length", MINIMUM_CONTIG_LENGTH),
       min_overlap=config['combine_contigs_params']['min_overlap'],
       max_substitutions=config['combine_contigs_params']['max_substitutions'],
       dont_allow_N= 't' if config['combine_contigs_params']['dont_allow_N'] else 'f',
       remove_cycles='t' if config['combine_contigs_params']['remove_cycles'] else 'f',
       trim_contradictions='t' if config['combine_contigs_params']['trim_contradictions'] else 'f'

    shell:
        """
            dedupe.sh in={params.input} findoverlaps cluster processclusters \
            out={output.combined_contigs} \
            csf={output.cluster_stats} \
            dot={output.dot} \
            minoverlap={params.min_overlap}\
            minscaf={params.min_length} \
            maxsubs={params.max_substitutions} \
            threads={threads} \
            sort=length \
            maxspanningtree={params.remove_cycles} \
            exact={params.dont_allow_N}\
            fixcanoncontradictions={params.trim_contradictions}\
            -Xmx{resources.mem}G 2> >(tee {log})
        """

# vizualize dot, takes enormous times
# dot -Tpdf combined_contigs_graph.dot -o combined_clusters.pdf

localrules: rename_combined_contigs

rule rename_combined_contigs:
    # standardizes header labels within contig FASTAs
    input:
        rules.combine_contigs.output.combined_contigs
    output:
        "{folder}/combined_contigs.fasta"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    params:
        prefix='C'
    shell:
        """rename.sh in={input} out={output} ow=t prefix={params.prefix}"""


rule combined_contigs_stats:
    input:
        rules.rename_combined_contigs.output
    output:
        "{folder}/combined_contigs_stats.txt"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        1
    resources:
        mem = config.get("java_mem", JAVA_MEM)
    shell:
        "stats.sh in={input} format=3 -Xmx{resources.mem}G > {output}"

# TODO: this my interfere with a rule in assemble
rule align_reads_to_combined_contigs:
    input:
        unpack(get_quality_controlled_reads),
        fasta = "{folder}/{Reference}.fasta",
    output:
        sam = temp("{folder}/sequence_alignment_{Reference}/{sample}.sam"),
        unmapped= expand("{{folder}}/sequence_alignment_{{Reference}}/unmapped/{{sample}}_unmapped_{fraction}.fastq.gz",fraction= MULTIFILE_FRACTIONS)
    benchmark:
        "logs/benchmarks/sequence_alignment_{Reference}/{sample}.txt"
    params:
        input= lambda wc,input : input_params_for_bbwrap(wc,input),
        maxsites = config.get("maximum_counted_map_sites", MAXIMUM_COUNTED_MAP_SITES),
        unmapped= lambda wc,output: "outu1={0},{2} outu2={1},null".format(*output.unmapped) if PAIRED_END else "outu={0}".format(*output.unmapped),
        max_distance_between_pairs=config.get('contig_max_distance_between_pairs',CONTIG_MAX_DISTANCE_BETWEEN_PAIRS),
        paired_only = 't' if config.get("contig_map_paired_only", CONTIG_MAP_PAIRED_ONLY) else 'f',
        ambiguous = 'all' if CONTIG_COUNT_MULTI_MAPPED_READS else'best',
        min_id= config.get('contig_min_id',CONTIG_MIN_ID),
    log:
        "{folder}/logs/sequence_alignment_{Reference}/{sample}.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("java_mem", JAVA_MEM)
    shell:
        """    bbwrap.sh nodisk=t \
               ref={input.fasta} \
               {params.input} \
               trimreaddescriptions=t \
               outm={output.sam} \
               {params.unmapped} \
               threads={threads} \
               pairlen={params.max_distance_between_pairs} \
               pairedonly={params.paired_only} \
               mdtag=t \
               xstag=fs \
               nmtag=t \
               sam=1.3 \
               local=t \
               ambiguous={params.ambiguous} \
               secondary=t \
               ssao=t \
               maxsites={params.maxsites} \
               -Xmx{resources.mem}G \
               2> {log}

               #max_distance_between_pairs : pairlen=32000           Set max allowed distance between paired reads.
               #(insert size)=(pairlen)+(read1 length)+(read2 length)
        """
rule pileup_combined_contigs:
    input:
        fasta = "{folder}/{Reference}.fasta",
        sam=temp("{folder}/sequence_alignment_{Reference}/{sample}.sam"),
    output:
        covstats = "{folder}/sequence_alignment_{Reference}/{sample}/{sample}_coverage.txt",
        basecov=temp("{folder}/sequence_alignment_{Reference}/{sample}/{sample}_base_coverage.txt.gz"),
        covhist= "{folder}/sequence_alignment_{Reference}/{sample}/{sample}_coverage_histogram.txt.gz",
        bincov ="{folder}/sequence_alignment_{Reference}/{sample}/{sample}_coverage_binned.txt.gz",
    params:
        pileup_secondary='t' if config.get("count_multi_mapped_reads",True) else 'f'
    benchmark:
        "logs/benchmarks/pileup_{Reference}/{sample}.txt"
    log:
        "{folder}/logs/pileup_{Reference}/{sample}.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("java_mem", JAVA_MEM)
    shell:
        """
            pileup.sh ref={input.fasta} in={input.sam} threads={threads} \
            -Xmx{resources.mem}G covstats={output.covstats} \
            hist={output.covhist} basecov={output.basecov} physcov secondary={params.pileup_secondary} bincov={output.bincov} 2>> {log}
        """

bam_combined_contigs_alignemnt= "contigs/sequence_alignment_combined_contigs/{sample}.bam"


localrules: combine_coverages_of_combined_contigs, combine_bined_coverages_of_combined_contigs

rule combine_coverages_of_combined_contigs:
    input:
        covstats = expand("{folder}/sequence_alignment_{Reference}/{sample}/{sample}_coverage.txt",
            sample=SAMPLES,Reference='combined_contigs',folder=combined_contigs_folder)
    output:
        "contigs/combined_median_coverage.tsv",
        "contigs/combined_readcounts.tsv",
        gc_stats = "contigs/combined_contigs_stats_gc.tsv"
    run:

        import pandas as pd
        import os

        combined_cov={}
        combined_N_reads={}
        for cov_file in input:

            sample= os.path.split(cov_file)[-1].split('_')[0]
            data= pd.read_table(cov_file,index_col=0)

            if cov_file == input[0]:
                data[['Length','Ref_GC']].to_csv(output.gc_stats,sep='\t')

            data.loc[data.Median_fold<0,'Median_fold']=0
            combined_cov[sample]= data.Median_fold
            combined_N_reads[sample] = data.Plus_reads+data.Minus_reads

        pd.DataFrame(combined_cov).to_csv(output[0],sep='\t')
        pd.DataFrame(combined_N_reads).to_csv(output[1],sep='\t')


rule combine_bined_coverages_of_combined_contigs:
    input:
        expand("{folder}/sequence_alignment_{Reference}/{sample}/{sample}_coverage_binned.txt.gz",
            sample=SAMPLES,Reference='combined_contigs',folder=combined_contigs_folder)
    output:
        "contigs/combined_coverage_binned.tsv.gz",
    run:

        import pandas as pd
        import os

        binCov={}
        for cov_file in input:

            sample= os.path.split(cov_file)[-1].split('_')[0]

            binCov[sample] = pd.read_table(cov_file,compression='gzip',comment='#',header=None,index_col=[0,2],usecols=[0,1,2],squeeze=True)

        binCov = pd.DataFrame(binCov)
        binCov.index.names=['Contig','Position']
        binCov.to_csv(output[0],sep='\t',compression='gzip')


#TODO detect read length automatically.
if config.get("perform_genome_binning", True):
    if config['combine_contigs_params']['binner']=='concoct':

        rule run_concoct:
          input:
              coverage= "{folder}/sequence_alignment_{Reference}/combined_median_coverage.tsv".format(Reference='combined_contigs',folder=combined_contigs_folder),
              fasta= "{folder}/{Reference}.fasta".format(Reference='combined_contigs',folder=combined_contigs_folder)
          output:
              expand("{folder}/binning/{file}",folder=combined_contigs_folder,file=['means_gt2500.csv','PCA_components_data_gt2500.csv','original_data_gt2500.csv','PCA_transformed_data_gt2500.csv','pca_means_gt2500.csv','args.txt','responsibilities.csv']),
          params:
              basename= lambda wc,output: os.path.dirname(output[0]),
              Nexpected_clusters= config['concoct']['Nexpected_clusters'],
              read_length= config['concoct']['read_length'],
              min_length=config["minimum_contig_length"],
              niterations=config["concoct"]["Niterations"]
          benchmark:
              "logs/benchmarks/binning/concoct.txt"
          log:
              "{folder}/binning/log.txt".format(folder=combined_contigs_folder)
          conda:
              "%s/concoct.yaml" % CONDAENV
          threads:
              10 # concoct uses 10 threads by default, wit for update: https://github.com/BinPro/CONCOCT/issues/177
          resources:
              mem = config.get("java_mem", JAVA_MEM)
          shell:
              """
                  concoct -c {params.Nexpected_clusters}\
                  --coverage_file {input.coverage}\
                  --composition_file {input.fasta}\
                  --basename {params.basename}\
                  --read_length {params.read_length} \
                  --length_threshold {params.min_length}\
                  --converge_out \
                  --iterations {params.niterations}
              """
# TODO: metabat conda recipie :-(
#     elif config['combine_contigs_params']['binner']=='metabat':
#         rule get_metabat_deph_file:
#               input:
#                     bam= expand(bam_combined_contigs_alignemnt,sample=SAMPLES)
#               output:
#                   expand("{folder}/binning/metabat_depth.txt",folder=combined_contigs_folder)
#               params:
#               log:
#                   "{folder}/binning/metabat.log".format(folder=combined_contigs_folder)
#               conda:
#                   "%s/metabat.yaml" % CONDAENV
#               threads:
#                   config['threads']
#               resources:
#                   mem = config.get("java_mem", JAVA_MEM)
#               shell:
#                     """
#                     jgi_summarize_bam_contig_depths --outputDepth {output} {input.bam} &> >(tee {log})
#                     """
#         rule run_metabat:
#             input:
#                 coverage= "{folder}/sequence_alignment_{Reference}/combined_median_coverage.tsv".format(Reference='combined_contigs',folder=combined_contigs_folder),
#                 fasta= "{folder}/{Reference}.fasta".format(Reference='combined_contigs',folder=combined_contigs_folder)
#             output:
#                 cluster_membership="{folder}/binning_cluster_membership.txt".format(folder=combined_contigs_folder),
#             params:
#                   sensitivity = 500 if config['binning_sensitivity']=='sensitive' else 200
#                   output_dir="annotations"
#             benchmark:
#                 "logs/benchmarks/binning/metabat.txt"
#             log:
#                 "logs/binning/metabat.txt"
#             conda:
#                 "%s/metabat.yaml" % CONDAENV
#             threads:
#                 config["threads"]
#             resources:
#                 mem = config.get("java_mem", JAVA_MEM)
#             shell:
#                   """
#                   metabat2 -i {input.contigs} \
#                   --abdFile {input.depth_file} \
#                   --minContig {params.min_contig_len} \
#                   --numThreads {threads} \
#                   --saveCls {output.cluster_membership} \
#                   --unbinned \
#                   --maxEdges {params.sensitivity} \
#                   -o {params.output_dir}/bin \
#                   &> >(tee {log})
#                   """
#

# https://bitbucket.org/berkeleylab/metabat/wiki/Best%20Binning%20Practices


    else:
        raise NotImplementedError("We don't have implemented the binning method: {}\ntry 'concoct'".format(config['combine_contigs_params']['binner']))
else:
    rule analyse_whole_metagenome:
        input:
            COMBINED_CONTIGS
        output:
            "annotations/metagenome/metagenome_contigs.fasta"
        shell:
            """
                ln -fs ../../{input} {output}
            """
# TODO: get absolute path
### GENE prediction

# if config['gene_predicter']=='prokka':
#
#     rule predict_genes:
#         input:
#             "annotations/{MAG}/{MAG}_contigs.fasta"
#         output:
#             faa = "annotations/{MAG}/predicted_genes/genes.faa",
#             ffn = "annotations/{MAG}/predicted_genes/genes.ffn",
#             fna = "annotations/{MAG}/predicted_genes/genes.fna",
#             fsa = "annotations/{MAG}/predicted_genes/genes.fsa",
#             gff = "annotations/{MAG}/predicted_genes/genes.gff",
#             log = "annotations/{MAG}/predicted_genes/genes.log",
#             tbl = "annotations/{MAG}/predicted_genes/genes.tbl",
#             tsv = "annotations/{MAG}/predicted_genes/genes.tsv",
#             txt = "annotations/{MAG}/predicted_genes/genes.txt"
# #            discrepancy = "annotations/{MAG}/prokka/genes.err",
#         benchmark:
#             "logs/benchmarks/prokka/{MAG}.txt"
#         params:
#             outdir = lambda wc, output: os.path.dirname(output.faa),
#             kingdom = config.get("prokka_kingdom", PROKKA_KINGDOM)
#         conda:
#             "%s/required_packages.yaml" % CONDAENV
#         threads:
#             config.get("threads", 1)
#         shell:
#             """prokka --outdir {params.outdir} \
#                    --force \
#                    --prefix genes \
#                    --locustag {wildcards.MAG} \
#                    --kingdom {params.kingdom} \
#                    --metagenome \
#                    --cpus {threads} \
#                    {input}
#             """

# elif config['gene_predicter']=='prodigal':
warnings.warn("gene_predicter=prodigal creates a gff which is not compatibl ewith downstream workflow.\ngenes don't follow naming rule C_0_1 C_0_2")
rule predict_genes:
    input:
        "annotations/{MAG}/{MAG}_contigs.fasta"
    output:
        fna="annotations/{MAG}/predicted_genes/genes.fna",
        faa="annotations/{MAG}/predicted_genes/genes.faa",
        gff="annotations/{MAG}/predicted_genes/genes.gff"

    conda:
        "%s/gene_catalog.yaml" % CONDAENV
    log:
        "logs/benchmarks/prodigal/{MAG}.txt"
    threads:
        1
    shell:
        """
            prodigal -i {input} -o {output.gff} -d {output.fna} -a {output.faa} -p meta -f gff 2> >(tee {log})
        """
# else:
#     raise NotADirectoryError("There is no genepredicter iplemented for: {} try one of ['prodigal','prokka']".format(config['gene_predicter']))


rule update_gene_table:
    input:
        "annotations/{MAG}/predicted_genes/genes.gff"
    output:
        "annotations/{MAG}/predicted_genes/genes_plus.tsv"
    shell:
        """atlas gff2tsv {input} {output}"""

localrules: renameeggNOG_annotation
rule renameeggNOG_annotation:
    input:
        "annotations/{MAG}/predicted_genes/genes.emapper.annotations"
    output:
        "annotations/{MAG}/eggNOG_annotation.tsv"
    shell:
        "cp {input} {output}"

include: "gene_annotation.snakefile"

# Taxonomy

rule MAG_run_diamond_blastp:
    input:
        fasta = "annotations/{MAG}/predicted_genes/genes.faa",
        db = config["diamond_db"]
    output:
        "annotations/{MAG}/refseq/hits.tsv"
    benchmark:
        "logs/benchmarks/run_diamond_blastp/{MAG}.txt"
    params:
        tmpdir = "--tmpdir %s" % TMPDIR if TMPDIR else "",
        top_seqs = config.get("diamond_top_seqs", DIAMOND_TOP_SEQS),
        e_value = config.get("diamond_e_value", DIAMOND_E_VALUE),
        min_identity = config.get("diamond_min_identity", DIAMOND_MIN_IDENTITY),
        query_cover = config.get("diamond_query_coverage", DIAMOND_QUERY_COVERAGE),
        gap_open = config.get("diamond_gap_open", DIAMOND_GAP_OPEN),
        gap_extend = config.get("diamond_gap_extend", DIAMOND_GAP_EXTEND),
        block_size = config.get("diamond_block_size", DIAMOND_BLOCK_SIZE),
        index_chunks = config.get("diamond_index_chunks", DIAMOND_INDEX_CHUNKS),
        run_mode = "--more-sensitive" if not config.get("diamond_run_mode", "") == "fast" else ""
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """diamond blastp \
               --threads {threads} \
               --outfmt 6 \
               --out {output} \
               --query {input.fasta} \
               --db {input.db} \
               --top {params.top_seqs} \
               --evalue {params.e_value} \
               --id {params.min_identity} \
               --query-cover {params.query_cover} \
               {params.run_mode} \
               --gapopen {params.gap_open} \
               --gapextend {params.gap_extend} \
               {params.tmpdir} \
               --block-size {params.block_size} \
               --index-chunks {params.index_chunks}"""


rule MAG_add_contig_metadata:
    input:
        hits = "annotations/{MAG}/refseq/hits.tsv",
        gff = "annotations/{MAG}/predicted_genes/genes.gff"
    output:
        temp("annotations/{MAG}/refseq/hits_plus.tsv")
    shell:
        "atlas munge-blast {input.hits} {input.gff} {output}"


rule MAG_sort_munged_blast_hits:
    # ensure blast hits are grouped by contig, ORF, and then decreasing by bitscore
    input:
        "annotations/{MAG}/refseq/hits_plus.tsv"
    output:
        "annotations/{MAG}/refseq/hits_plus_sorted.tsv"
    shell:
        "sort -k1,1 -k2,2 -k13,13rn {input} > {output}"


rule MAG_parse_blastp:
    # assign a taxonomy to contigs using the consensus of the ORF assignments
    input:
        "annotations/{MAG}/refseq/hits_plus_sorted.tsv"
    output:
        "annotations/{MAG}/refseq/tax_assignments.tsv"
    params:
        namemap = config["refseq_namemap"],
        treefile = config["refseq_tree"],
        summary_method = config.get("summary_method", SUMMARY_METHOD),
        aggregation_method = config.get("aggregation_method", AGGREGATION_METHOD),
        majority_threshold = config.get("majority_threshold", MAJORITY_THRESHOLD),
        min_identity = config.get("diamond_min_identity", DIAMOND_MIN_IDENTITY),
        min_bitscore = config.get("min_bitscore", MIN_BITSCORE),
        min_length = config.get("min_length", MIN_LENGTH),
        max_evalue = config.get("diamond_e_value", DIAMOND_E_VALUE),
        max_hits = config.get("max_hits", MAX_HITS),
        top_fraction = (100 - config.get("diamond_top_seqs", 5)) * 0.01
    shell:
        """atlas refseq \
               --summary-method {params.summary_method} \
               --aggregation-method {params.aggregation_method} \
               --majority-threshold {params.majority_threshold} \
               --min-identity {params.min_identity} \
               --min-bitscore {params.min_bitscore} \
               --min-length {params.min_length} \
               --max-evalue {params.max_evalue} \
               --max-hits {params.max_hits} \
               --top-fraction {params.top_fraction} \
               {input} \
               {params.namemap} \
               {params.treefile} \
               {output}"""

# if config.get("perform_genome_binning", True):
#     rule MAG_merge_sample_tables:
#         input:
#             prokka = "{sample}/annotation/prokka/{sample}_plus.tsv",
#             refseq = "annotations/{MAG}/refseq/tax_assignments.tsv",#
#             counts = "{sample}/annotation/feature_counts/{sample}_counts.txt",
#             completeness = "{sample}/genomic_bins/checkm/completeness.tsv",
#             taxonomy = "{sample}/genomic_bins/checkm/taxonomy.tsv"
#         output:
#             "{sample}/{sample}_annotations.txt"
#         params:
#             fastas = lambda wc: " --fasta ".join(glob("{sample}/genomic_bins/{sample}.*.fasta".format(sample=wc.sample)))
#         shell:
#             "atlas merge-tables \
#                  --counts {input.counts} \
#                  --completeness {input.completeness} \
#                  --taxonomy {input.taxonomy} \
#                  --fasta {params.fastas} \
#                  {input.prokka} \
#                  {input.refseq} \
#                  {output}"
#
#
# else:
rule MAG_merge_sample_tables:
    input:
        genes = "annotations/{MAG}/predicted_genes/genes_plus.tsv",
        refseq = "annotations/{MAG}/refseq/tax_assignments.tsv",
        counts = "annotations/{MAG}/feature_counts/gene_info.tsv"
    output:
        "annotations/{MAG}/annotations.txt"
    shell:
        "atlas merge-tables \
             --counts {input.counts} \
             {input.genes} \
             {input.refseq} \
             {output}"

rule counts_genes:
    input:
        gtf = "annotations/{MAG}/predicted_genes/genes.gtf",
        bam = bam_combined_contigs_alignemnt
    output:
        summary = "annotations/{MAG}/feature_counts/{sample}_counts.txt.summary",
        counts = "annotations/{MAG}/feature_counts/{sample}_counts.txt"
    params:
        min_read_overlap = config.get("minimum_region_overlap", MINIMUM_REGION_OVERLAP),
        paired_only= "-B" if config.get('contig_map_paired_only',CONTIG_MAP_PAIRED_ONLY) else "",
        paired_mode = "-p" if PAIRED_END else "",
        multi_mapping = "-M --fraction" if config.get("contig_count_multi_mapped_reads",CONTIG_COUNT_MULTI_MAPPED_READS) else "--primary",
        feature_counts_allow_overlap = "-O --fraction" if config.get("feature_counts_allow_overlap", FEATURE_COUNTS_ALLOW_OVERLAP) else ""
    log:
        "{sample}/logs/counts_per_region.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """featureCounts \
                --minOverlap {params.min_read_overlap} \
                {params.paired_mode} \
                {params.paired_only} \
               -F GTF \
               -T {threads} \
               {params.multi_mapping} \
               {params.feature_counts_allow_overlap} \
               -t CDS \
               -g ID \
               -a {input.gtf} \
               -o {output.counts} \
               {input.bam} 2> {log}"""

rule combine_gene_counts:
    input:
        expand("annotations/{{MAG}}/feature_counts/{sample}_counts.txt",
            sample=SAMPLES)
    output:
        'annotations/{MAG}/gene_counts.tsv',
        'annotations/{MAG}/feature_counts/gene_info.tsv'
    run:
        import pandas as pd
        import os
        C= {}

        for file in input:
            D= pd.read_table(file,index_col=0,comment='#')
            # contigs/sequence_alignment_combined_contigs/S1/S1.bam
            sample= D.columns[-1].split('/')[-2]
            C[sample]= D.iloc[:,-1]
        C= pd.DataFrame(C)
        C.to_csv(output[0],sep='\t')

        D.iloc[:,:-1].to_csv(output[1],sep='\t')
