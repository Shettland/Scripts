#!/bin/bash

## RUNNING GET HOMOLOGUES ##


nohup ~/get_homologues/get_homologues-x86_64-20210217/get_homologues.pl -d fasta_files/ -n 4  &> output.log &
waiter="$!"

nohup ~/get_homologues/get_homologues-x86_64-20210217/get_homologues.pl -d fasta_filesM/ -M -A -P -n 4  &> outputM.log &
waiter="$waiter $!"

nohup ~/get_homologues/get_homologues-x86_64-20210217/get_homologues.pl -d fasta_filesG/ -G -A -P -n 4  &> outputG.log &
waiter="$waiter $!"

wait $waiter
echo "FINISHED GET_HOMOLOGUES"

nohup ~/get_homologues/get_homologues-x86_64-20210217/compare_clusters.pl -o trichinella_intersection -m -d \
fasta_filesM_homologues/trichinellapapuae_f0_alltaxa_algOMCL_e0_,\
fasta_files_homologues/trichinellapapuae_f0_alltaxa_algBDBH_e0_,\
fasta_filesG_homologues/trichinellapapuae_f0_alltaxa_algCOG_e0_ \
&> compare.log &


: '
~/get_homologues/get_homologues-x86_64-20210217/plot_matrix_heatmap.sh -i trichinella_intersection/pangenome_matrix_t0.tab -o pdf \
 -r -H 8 -W 14 -m 28 -t "sample pangenome (clusters=180)" -k "genes per cluster"

'
