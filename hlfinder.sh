#!/bin/bash

helpFunction()
{
   echo ""
   echo "How to run script: $0 -p path/to/main_folder -f path/to/geneIDs.file -t /path/to/FPKM_table.tab"
   echo -e "\t -p: Folder where get_homologues was executed so the script can find the *taxa* subfolder with clustered .faa files."
   echo -e "\t -f: File with the reference IDs to be found, one per line. Example file.txt:"
   echo -e "\t -t: Reference Table with homologues + FPKM"
   echo -e "\t T01_"
   echo -e "\t T03_"
   echo -e "\t T4A_"
   echo -e "\t ... with this file, $0 will look for matches based on geneIDs T01_* T03_* or T4A_*"
   exit 1 # Exit script after printing help if no parameter is indicated
}

while getopts "p:f:t:" opt
do
   case "$opt" in
        p ) path="$OPTARG" ;;
        f ) idfile="$OPTARG" ;;
	t ) gtable="$OPTARG" ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

list=(00-proteins_homologues/schistosomahaematobium_f0_alltaxa_algBDBH_e0_/
00-proteins_M_homologues/schistosomahaematobium_f0_alltaxa_algOMCL_e0_/
00-proteins_G_homologues/schistosomahaematobium_f0_alltaxa_algCOG_e0_/)

BLUE='\033[1;34m'
NC='\033[0m'
YELLOW='\033[1;33m'
RED='\033[1;31m'

if [ -z "$path" ] || [ -z "$idfile" ] || [[ -z "$gtable" ]];
then
   echo -e "${RED} Some or all of the parameters are empty ${NC} \n"
   echo "directories: example ${list[@]}"
   helpFunction
fi

if [[ ! -f $idfile ]] || [[ ! -s $idfile ]];
then
	echo -e "${RED} geneID.file given ($idfile) is empty or does not exist ${NC}"
	exit 1
fi

if test $(find ${path}* -maxdepth 1 -type d | grep "taxa_alg" | wc -l) -eq 0;
then
	echo -e "${RED} No cluster /*taxa* folders found in main_folder's subdirectories ${NC}"
	exit 1
fi
### Another possible way to exit if find doesn't returns 0: [[ ! -z `find ${path}] -maxdepth 1 -type d | grep "alltaxa"` ]] && helpfunction

DIR=$(dirname "$(readlink -f "$0")")

: '
if test $(find $DIR/ -maxdepth 1 -type f | grep "merged_table" | wc -l) -eq 0;
then
	echo -e "${RED} No reference gene table found in $DIR ${NC}"
	exit 1
fi

'

names=()
readarray -t names < $idfile


folders=()
folders=$(find ${path}* -maxdepth 2 -type d | grep "taxa_alg" > tempfold.txt)
readarray -t folders < tempfold.txt
rm tempfold.txt



if [[ -f $DIR/badgenes.txt ]];
then
	rm $DIR/badgenes.txt
fi

touch badgenes.txt

for (( i=0; i<${#names[@]}; i++ ));
do

	for folder in ${folders[@]};
	do

		for file in $folder/*.faa;
		do
#gene=$(grep "gene=${names[c]}" $file > badgenes.txt)#
#sed -n "s|matches found in folder ${FOLDER}|&|p" >> badgenes.txt

			if [[ -f $DIR/badgenes.txt ]];
			then
				gene=$(sed -n "/gene=${names[i]}/p" $file >> badgenes.txt)
				[[ ! -z $gene ]] && loge=$(echo "$gene found in file: $file" >> badgenlocs.log)
			else

				gene=$(sed -n "/gene=${names[i]}/p" $file >> badgenes.txt)
				[[ ! -z $gene ]] && loge=$(echo "$gene found in file: $file" >> badgenlocs.log)
			fi

		done
	done

done

echo -e "\n" >> badgenes.txt
filter=$(awk '!seen[$3]++' badgenes.txt > badgenf.txt | echo " " >> badgenf.txt) ###removing duplicates based on column 3 (gene=...) and adding an empty line for later awk process

rm badgenes.txt

### sed -i "s/.*gene=.*/file found = file &/{p;n;p}" badgenes.txt
### sed -i "/gene/ s/^/file-found=${file}/" ###
### sed -i "s/.*gene=.*/file-found='"${file}"' &/" badgenes.txt

cutting=$(awk 'NF' badgenf.txt | cut -d " " -f 3 | sed 's/^gene=//' | uniq > badgeneids.txt)  ### removing empty lines and substring *gene=* to get a clear geneID.file

if [[ -f $DIR/badgenes_table.txt ]];
then
	rm $DIR/badgenes_table.txt
fi

touch badgenes_table.txt

### FINDING MATCHES BETWEEN GENES FROM GENEID-FILE AND REFERENCE TABLE ###
colnums=$(head -1 $gtable | tr -s '\t' '\n' | nl -nln | cut -f1 > colfile.txt)
ballcols=$(head -1 $gtable | tr -s '\t' '\n' | nl -nln |  grep "ballgown" | cut -f1 | tr "\n" "|" | tr -d " ")
truecols=$(grep -vE ${ballcols%?} colfile.txt | xargs | tr " " ",") 	### Taking every column except those with 'ballgown' on it ###
maxcols=$(head -1 $gtable | tr "\t" "\n" | wc -l)

while read line;
do
	finder=$(grep "$line" $gtable | cut -f ${truecols} >> badgenes_table.txt)
done < badgeneids.txt

sorting=$(awk -F '\t' '{print $NF,$0}' badgenes_table.txt | sort -gr | cut -f2- -d' ' |  uniq > badgenes2_table.txt) ### Sorting by 8th column and deleting duplicates
### sort -t$'\t' -k 8,8 -rg goodgen_table.txt | uniq ###

badgenes=$(cat badgenes2_table.txt | wc -l)
echo "Number of genes to be filtered (matching) = $badgenes"

rm badgenes_table.txt

### CREATING FILTERED TABLE FROM MATCHING GENES ###

string=""

while read line;
do
	string="$string|$line"
done < badgeneids.txt

goodgenes=$(tail -n +2 $gtable | grep -vE "${string#?}" | cut -f 1-${maxcols} > goodgen_table.txt) ## grep -E needed to escape '\|' | tail -n +2 to grep from everything but the first line (headers)

refcol=$(head -1 $gtable | tr "\t" "\n" | nl -nln | cut -f 1 | tr "\n" "\t" | rev | cut -f 2 | rev)

sorting=$(awk -F '\t' '{print $NF,$0}' goodgen_table.txt | sort -rg | cut -f2- -d' ' | uniq > $DIR/goodgenes_table.txt) ### Sorting by last column and deleting duplicates

###sorting=$(sort -t$'\t' -k$refcol -nr goodgen_table.txt | uniq > $DIR/goodgenes_table.txt)
headers=$(head -1 $gtable | cut -f 1-${maxcols})
addhead=$(sed -i "1s/^/${headers}\n/" goodgenes_table.txt)

rm goodgen_table.txt
rm colfile.txt


### sort -t$'\t' -k 8,8 -rg goodgen_table.txt | uniq ###


num=$(cat goodgenes_table.txt | wc -l)
sum=$(( $num - 1 ))

fpkms=$(head -1 $gtable | tr "\t" "\n" | nl -nln | grep "FPKM" | cut -f 1 | tr "\n" "," | xargs | tr -d " " )

count=0;
total=0;

prep=$(echo "$fpkms" | tr "," " ")
colarray=($prep)

for col in ${colarray[@]};
do

        for i in $( awk -v col=$col -F "\t" 'NR>=2 { print $col; }' $gtable )
        do
                total=$(echo $total+$i | bc )
                ((count++))
        done

        avg=$(echo "scale=6; $total / $count" | bc)
        meanrow="${meanrow}\t${avg}"

done

avgcols=$(echo "$meanrow" | sed 's/^[[:space:]]*//')
truenums=$(head -1 $gtable | tr "\t" "\n" | nl -nln | grep -v "FPKM" | wc -l)

for (( i=1 ; i<${truenums} ; i++ ))
do
	tabnum="${tabnum}\t-"
done
avgwrite=$(echo -e "Averages:${tabnum}${meanrow}"  >> $DIR/goodgenes_table.txt)


echo " "
echo "Filtered gene table saved as goodgenes_table.txt"
echo " "
echo "Total genes saved: $sum"


