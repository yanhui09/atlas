#! /bin/bash
set -euo pipefail


atlas --version


databaseDir=".test/databases"


ressource_args=" --config java_mem=10 assembly_mem=10"


#git clone https://github.com/metagenome-atlas/example_data.git


WD="example_data/binning"
reads_dir="example_data/binning_reads"

rm -f $WD/samples.tsv
#
atlas init --db-dir $databaseDir --threads 3  -w $WD $reads_dir --skip-qc

touch $WD/finished_assembly

atlas run binning -w $WD $ressource_args assembler=spades final_binner=metabat interleaved_fastq=true $@

# genomes need databases
# atlas run genomes -w $WD $ressource_args assembler=spades final_binner=metabat --omit-from get_genome_for_cat $@

#atlas run all -w $WD $ressource_args assembler=spades final_binner=metabat $@
