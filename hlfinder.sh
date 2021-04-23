#!/bin/bash

: '
helpFunction()
{
   echo ""
   echo "How to run script: $0 -p path/to/samples_directory/ -t path/to/trimmomatic.jar -a trimmomatic_adapters.file -y path/to/Trinity"
   echo -e "\t -p: folder where the samples are located. Also works for samples within subdirectories"
   echo -e "\t -t: Path to trimmomatic.jar executable"
   echo -e "\t -a: Adapters file for trimmomatic"
   echo -e "\t -y: Path to trinity executable"
   echo -e "\t Paired End samples only. Single end is not supported"
   exit 1 # Exit script after printing help if no parameter is indicated
}

while getopts "p:a:t:y:" opt
do
   case "$opt" in
        p ) path="$OPTARG" ;;
        a ) adapter="$OPTARG" ;;
        t ) trimmomatic="$OPTARG" ;;
        y ) trinity="$OPTARG" ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

if [ -z "$path" ] || [ -z "$trimmomatic" ] || [ -z "$adapter" ] || [ -z "$trinity" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

'

FOLDER=$1

list=(00-proteins_homologues/schistosomahaematobium_f0_alltaxa_algBDBH_e0_/
00-proteins_M_homologues/schistosomahaematobium_f0_alltaxa_algOMCL_e0_/
00-proteins_G_homologues/schistosomahaematobium_f0_alltaxa_algCOG_e0_/)


if [[ -z $FOLDER ]];
then
echo "directories: ${list[@]}"
echo "script usage: $0 folder_with_alltaxa_directories"
echo "works for subdirectories too (maxdepth = 1)"
else

names=(T01_ T03_ T4A_ T4E_ T4C_ T4B_)
folders=()

fold=$(find ${1}* -maxdepth 1 -type d | grep "alltaxa" > tempfold.txt)
readarray -t folders < tempfold.txt

touch badgenes.txt

for (( i=0; i<${#names[@]}; i++ ));
do

	for folder in ${folders[@]};
	do

		for file in $folder/*.faa;
		do
#gene=$(grep "gene=${names[c]}" $file > badgenes.txt)#
#sed -n "s|matches found in folder ${FOLDER}|&|p" >> badgenes.txt

			if [[ -f ./badgenes.txt ]];
			then
				gene=$(sed -n "/gene=${names[i]}/p" $file >> badgenes.txt)
			else

				gene=$(sed -n "/gene=${names[i]}/p" $file >> badgenes.txt)
			fi

		done
	done

done


echo -e "\n" >> badgenes.txt
filter=$(awk '!seen[$3]++' badgenes.txt > badgenf.txt | echo " " >> badgenf.txt) ###removing duplicates based on column 3 (gene=...) and adding last empty column for later process

### sed -i "s/.*gene=.*/file found = file &/{p;n;p}" badgenes.txt
### sed -i "/gene/ s/^/file-found=${file}/" ###
### sed -i "s/.*gene=.*/file-found='"${file}"' &/" badgenes.txt

cutting=$(awk 'NF' badgenf.txt | cut -d " " -f 3 | sed 's/^gene=//' | uniq > tempfile.txt) #echo " " >> tempfile.txt)   ### removing empty lines and substring *gene=* to get a clear geneID.file
rm badgenes_table.txt
touch badgenes_table.txt

### FINDING MATCHES BETWEEN GENES FROM GENEID-FILE AND REFERENCE TABLE ###

while read line;
do

finder=$(grep "$line" ult_FPKM_merged_table.txt | cut -f 1-7,25 >> badgenes_table.txt)

done < tempfile.txt

sorting=$(sort -t$'\t' -k 8,8 -rg badgenes_table.txt | uniq > badgenes2_table.txt) ### Sorting by 8th column and deleting duplicates

badgenes=$(cat badgenes2_table.txt | wc -l)
echo "Number of genes to be filtered (matching) = $badgenes"

rm badgenes_table.txt


### CREATING FILTERED TABLE FROM MATCHING GENES ###

string=""

while read line;
do

string="$string|$line"

done < tempfile.txt

string=${string#?}

goodgenes=$(grep -vE "$string" ult_FPKM_merged_table.txt | cut -f 1-7,25 > goodgen_table.txt) ## grep -E needed so '\|' is not needed

sorting=$(sort -t$'\t' -k 8,8 -rg goodgen_table.txt | uniq > goodgenes_table.txt) ### Sorting by 8th column and deleting duplicates
headers=$(head -1 ult_FPKM_merged_table.txt | cut -f 1-7,25)
sed -i "1s/^/${headers}\n/" goodgenes_table.txt

rm goodgen_table.txt

num=$(cat goodgenes_table.txt | wc -l)
sum=$(( $num - 1 ))

echo " "
echo "Filtered gene table saved as goodgenes_table.txt"
echo " "
echo "Total genes saved: $sum"


fi


