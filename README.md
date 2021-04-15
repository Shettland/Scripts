# Scripts
This is a place where all my scripts are saved.
- R scripts:
    - [merge_results_fin.R](#Merge_results_finr)
- Bash scripts:
    - [pipeline.sh](#pipelinesh)
    - [lengther.sh](#lengthersh)
## R Scripts 

### Merge_results_fin.R 

__*Merge_results_fin.R*__ is a script that starts with a list of samples and a raw table with gene homologues. 

**The script does various things:**
1. Edits the raw **gene_table** to delete redundant information and present a clean table to work with.
2. Runs **ballgown** on all the samples given, creating one new table for each of those results my merging them with the gene_table.
3. Takes the recently made tables and merges all of them into one final_table.
4. Adds valuable information to the final_table and order to ease its comprehension.

For this task it uses ***for*** loops and ***if*** conditionals to iterate over the data, generating the results automatically.

**Warning:** *The script doesn't include previous RNAseq analysis which are necessary to generate .ctab files that the script needs to work.
These are:*
- ***Quality analysis:** Using FastQC.*
- ***Trimming:** Using trimmomatic.*
- ***Alignment:** Using HISAT2 or STAR.*
- ***Post Processing:***
    - ***Marking and removing duplicates:** Using Preseq and dupRadar.*
    - ***Generating FPKM tables:** Using featureCounts and Stringtie (which includes ballgown).*

*Later analysis could be: differential expression analysis (Deseq2), plotting (EdgeR or ggplots) and MultiQC for quality analysis.*

**Required libraries:** ***tidyverse, ballgown, reshape, stringr***

## Bash scripts

### pipeline.sh

A pipeline to run trimmomatic and trinity over all the paired-end samples (2 *.fastq* files each) found in a given directory or subdirectories. Single-end and stranded are not supported, you might need to change the trinity parameters if you wish to use those kind of samples.
Example:
~~~
bash ./pipeline.sh -p home/Analysis/RNAseq_samples/ -t home/share/trimmomatic/trimmomatic.jar -a trimmomatic_adapters.file -y home/share/trinity/Trinity  
~~~

**Warning:** *Trimmomatic and Trinity parameters are predefined to work for a wide variety of species from different reigns. If you wish you can change them directly from the script. Also keep in mind that Trinity's De novo assembly is a high consuming process: it might take 1 hour for every million reads* 

### lengther.sh

A script that selects **-n** number of genes from a gene_name table and finds their sequence searching through fasta files in a given directory, printing its length and the sum of total aminoacids of the selected gene-sequences. Its useful to know if your request to other software like secretomeP will surpass the permitted limit or not.

~~~ 
bash ./lengther.sh -f path/to/gene_table.file -p path/to/fastas_directory/ -n num -m path/to/pangenome_matrix.file
~~~

**Warning:**  *the pathing used in the script is not global, the script might not work if you run it through different locations.*
