# Scripts
This is a place where all my scripts are saved.

## R Scripts 

### Merge_results_fin.R

__*Merge_results_fin.R*__ is a script that starts with a list of samples and a raw table with gene homologues. 

**The script does various things:**
1. Edits the raw **gene_table** to delete redundant information and present a clean table to work with.
2. Runs ballgown on all the samples given, creating one new table for each of those results my merging them with the gene_table.
3. Merges all the recently made tables and merges all of them into one final_table.
4. Adds valuable information to the final_table and order to ease its comprehension.

For this task it uses ***for*** loops and ***if*** conditionals to iterate over the data, generating the results automatically.

## Bash scripts

### pipeline.sh

A pipeline to run trimmomatic and trinity over all the paired-end found in given directory or subdirectories.

### lengther.sh

A script that selects **n** number of genes from a gene_name table and finds their sequence searching through fasta files in a given directory, printing its length and the sum of total aminoacids of the selected gene-sequences. Its useful to know if your request to other software like secretomeP will surpass the permitted limit or not.
