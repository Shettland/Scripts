# R_Scripts
This is a place where R scripts are saved for later
## Merge_results_fin.R

**Merge_results_fin.R** is a script that starts with a list of samples and a raw table with gene homologues. 

**The script does various things:**
1. Edits the raw **gene_table** to delete redundant information and present a clean table to work with.
2. Runs ballgown on all the samples given, creating one new table for each of those results my merging them with the gene_table.
3. Merges all the recently made tables and merges all of them into one final_table.
4. Adds valuable information to the final_table and order to ease its comprehension.


