#!/usr/bin/bash
# encoding: utf-8
# author: Jinxin Meng
# created date: 2023-11-10, 11:33:11
# modified date: 2024-03-10, 11:33:11

[ $# -ne 2 ] && echo "Usage: $0 [gtdbtk_classify_directory] [output_file]" && exit -1

in=$1; out=$2
[ ! $(which gtdbtk 2>/dev/null) ] && source activate /share/data1/software/miniconda3/envs/gtdbtk >/dev/null
/share/data1/software/miniconda3/envs/gtdbtk/bin/gtdb_to_ncbi_majority_vote.py \
    --gtdbtk_output_dir $in \
    --ar53_metadata_file /share/data1/database/gtdb/taxonomy/ar53_metadata.tsv \
    --bac120_metadata_file /share/data1/database/gtdb/taxonomy/bac120_metadata.tsv \
    --output_file $out
conda deactivate
