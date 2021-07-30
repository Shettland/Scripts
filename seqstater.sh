#!/bin/bash

helpFunction()
{
   echo ""
   echo "How to run script: $0 -p path/to/fasta_folder -g path/to/clean_gene.table -t /path/to/FPKM_full_table.tab"
   echo -e "\t -p: Fasta_files folder, containing the original sequences from the genes in .faa files"
   echo -e "\t -g: Clean table containing the genes needed to search for its sequences"
   echo -e "\t -t: Full Reference Table with every homologue + FPKM. Used to calculate average FPKM"
   echo -e "\t Tables given by -g & -t could be the same in case you are extracting directly from the main table"
   exit 1 # Exit script after printing help if no parameter is indicated
}

while getopts "p:g:t:" opt
do
   case "$opt" in
        p ) path="$OPTARG" ;;
        f ) gtable="$OPTARG" ;;
        t ) fulltab="$OPTARG" ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

BLUE='\033[1;34m'
NC='\033[0m'
YELLOW='\033[1;33m'
RED='\033[1;31m'

if [ -z "$path" ] || [ -z "$gtable" ] || [[ -z "$fulltab" ]];
then
   echo -e "${RED} Some or all of the parameters are empty ${NC} \n"
   echo "directories: example ${list[@]}"
   helpFunction
fi

###gtable=goodgenes_table.txt
###fulltab=gfin_FPKM_merged_table.txt

taxas=$(head -1 $gtable | tr "\t" "\n" | nl -nln | grep -v "FPKM" | cut -f 1 | tr "\n" "," | xargs | tr -d " " )

cols=$(tail -n +2 $gtable | cut -f ${taxas%?} | tr "\t" " " | tr -d '"' > trygenes.txt)

fpkms=$(head -1 $gtable | tr "\t" "\n" | nl -nln | grep "FPKM" | cut -f 1 | tr "\n" "," | xargs | tr -d " " )

fpkmheaders=$(head -1 $gtable | cut -f ${fpkms%?} | sed 's/_ballgown//g' | tr -d '"')

if [[ -f goodseqs.txt ]];
then
	rm goodseqs.txt
fi

if [[ -f prove.log ]];
then
        rm prove.log
fi


count=0;
total=0;

prep=$(echo "$fpkms" | tr "," " ")
colarray=($prep)

for col in ${colarray[@]};
do

	for i in $( awk -v col=$col -F "\t" 'NR>=2 { print $col; }' $fulltab )
   	do
     		total=$(echo $total+$i | bc )
     		((count++))
   	done

	avg=$(echo "scale=6; $total / $count" | bc)
	meanrow="${meanrow}\t${avg}"

done

avgcols=$(echo "$meanrow" | sed 's/^[[:space:]]*//')

for cluster in $path*.faa;
do

	while read line;
	do

	greppats=""
	awkpats=""
	leng=$(echo "$line" | tr " " "\n" | wc -l)
        for (( i=1; i<=$leng; i++ ));
        do
                subst[i]=$(echo "$line" | cut -d " " -f $i | sed 's/^/>/')
                string[i]="${subst[i]}"
		greppats="${greppats}|$string[i]"
		awkpats="/$awkpats/ && /$string[i]/"
        done

        a=$(echo "$subst1" | echo "$subst2" | echo "$subst3" >> prove.log)


	awkpats="${awkpats#?????}"
        greper=$(grep "$greppats" $cluster | tr "\n" " " > greptest.txt)

	if test $(awk -v patterns=$awkpats '$patterns' greptest.txt | wc -l) -eq 1;
        then

        	stats=$(grep "${string}" $gtable | cut -f ${fpkms%?})

        	echo -e "\tCluster number $n stats: \n" >> goodseqs.txt
        	echo -e "Samples:\t$fpkmheaders" >> goodseqs.txt
        	echo -e "Stats:\t$stats" >> goodseqs.txt
	        echo -e "Average:$avgcols \n" >> goodseqs.txt
        	echo -e "\tSequences:" >> goodseqs.txt

		for (( i=1; i<=$leng; i++ ));
                do
                        maxa=$(grep ">$subst[i]" $cluster | wc -l)
                        for (( n=1; n<=$maxa; n++ ));
                        do
                                partsub=(grep ">$subst[i]" $cluster | cut -d " " -f 1 | tr "\n" "\t" | cut -f $n)
                                greper=$(sed -n -e "/$partsub/,/gene/{/$partsub/p;/gene/!p;}" $cluster >> goodseqs.txt)
                        done
                done
	fi

	done < trygenes.txt
done


###trim=$(sed -e '$!N;/^>.*\n>/D' -e 'P;D' goodseqs.txt > goodfseq.txt) # In case duplicated transcripts were found, they could be removed with this command
###trim2=$(sed -e '$!N;/^>.*\n-/D' -e 'P;D' goodfseq.txt > goodfinseq.txt) #In case script failed and printed several '-' lines they'd be removed with this command

echo "done"


