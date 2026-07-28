#!/bin/bash

cd /home/yoyerush/yo/methylome_pipeline/MethylC-analyzer_test/

#for CNTX in 'CG' 'CHG' 'CHH'; do
#
#if [ "$CNTX" = "CG" ]; then
#    min_proportion=0.4
#elif [ "$CNTX" = "CHG" ]; then
#    min_proportion=0.2
#elif [ "$CNTX" = "CHH" ]; then
#    min_proportion=0.1
#fi

python MethylC-analyzer/scripts/MethylC_yo_edited.py all samples_table_ddm1_test.txt TAIR10_GFF3_genes.gff /home/yoyerush/yo/methylome_pipeline/MethylC-analyzer_test/ -a ddm1 -b wt -r 100 # -context "$CNTX" -dmrc "$min_proportion"
#done