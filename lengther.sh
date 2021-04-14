#!/bin/bash

## SETTING HELP FUNCTION ##

helpFunction()
{
   echo ""
   echo "How to run script: $0 -f path/to/gene_table.file -p path/to/clusters_directory/ -n num -m path/to/pangenome_matrix.file"
   echo -e "\t -p: folder where the sequence clusters are located."
   echo -e "\t -f: Path to table with genenames/homologues"
   echo -e "\t -n: Number of top expresed sequences to measure"
   echo -e "\t -m: Path to pangenome matrix with cluster names added to genenames"
   exit 1 # Exit script after printing help if no parameter is indicated
}

while getopts "p:f:n:m:" opt
do
   case "$opt" in
        p ) path="$OPTARG" ;;
        f ) table="$OPTARG" ;;
	n ) num="$OPTARG" ;;
	m ) pangemat="$OPTARG" ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

if [ -z "$path" ] || [ -z "$table" ] || [ -z "$num" ] || [ -z "$pangemat" ];
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

NC='\033[0m'
YELLOW='\033[1;33m'
RED='\033[1;31m'

max=$(tail -n +2 $table | wc -l)

if [[ $num -gt $max ]];
then
echo -e "${RED}Number of genes requested excedes the number of sequences in file ${NC}"
echo -e "Max value =${YELLOW} $max ${NC}"
exit 1
fi

names=()
sed -n "2,${num}p" $table | cut -f 4 > filetemp.txt
readarray -t names < filetemp.txt

rm filetemp.txt


clusters=()
c=0

for name in ${names[@]}
do
	gene="${name%\"}"
	gene="${gene#\"}"

	greping=$(sed -n "s/${gene}\b/${gene}\b/p" $pangemat | cut -f 2 | sed -e 's/^"//' -e 's/"$//' >> filetemp2.txt)

	opt[c]=${gene}
	((c++))
done

readarray -t clusters < filetemp2.txt
rm filetemp2.txt

c=0
b=0

echo "selected genes: ${opt[@]}"

sequence=()
number=()

for cluster in ${clusters[@]}
do


	file=$(find $path | grep ${cluster})

	echo -e "${YELLOW} length of gene ${opt[c]} in cluster $cluster: ${NC}"
	length=$(sed -n "/${opt[c]}/{n;p;}" $file | awk "{print length}")
	echo "$length"

	line=(${length[@]})
	line=(${length[0]})


	if [[ $line -gt 250 ]];
	then
		echo -e "${RED} ${opt[c]}'s sequence length is >250!! ${NC}"
		sequence[${#sequence[@]}]=" Gene ${opt[c]}:length=${line} "
	fi


	values=(${length[@]})

	a=1
	e=0
	p=0
	for (( i= 0; i< ${#values[@]}; ++i ))
	do
		numbers[${#numbers[@]}]=${values[p]}

		if [[ $a -gt 1 && $a -lt 3 ]];
		then
			echo -e "${RED} Gene ${opt[c]} has more than one sequence!! ${NC}"
			duplicate[${#duplicate[@]}]=${opt[c]}
		fi
	((p++))
	((a++))
	((e++))
	done

	if [[ $e -gt 1 ]];
	then
		count[b]=$e
		((b++))
	fi

((c++))
done

echo -e "${YELLOW}Duplicated genes: ${duplicate[@]} ${NC}"
echo "${count[@]}"
### printf "%s\n" "${duplicate[@]}"  > dup_genes.txt ###

rm -f dup_genes.txt
touch dup_genes.txt

b=0
for dup in ${duplicate[@]};
do
printf "%s\n" "$dup is duplicated ${count[b]} times" >> dup_genes.txt
#sed -e "s/$/repeats ${count[b]} times/" -i dup_genes.txt
((b++))
done



echo ""
echo -e "${YELLOW} list of proteins with sequences greater than 250: ${NC}"
echo "${sequence[@]}"
echo ""

echo "all lengths: ${numbers[@]}"


for (( i = 0; i < ${#numbers[@]}; ++i ));
do
	sum=`expr $sum + ${numbers[i]}`
	### también sirve sum=$((sum + numbers[i])) ###
done

echo ""
echo -e "${YELLOW} total aminoacids = $sum ${NC}"
echo ""
echo "FIN"



#head -n 50 pangenome_matrix_edited.tab | cut -d "\t" -f 4
#sed -n '/${opt[c]}/{n;p;}'
#grep "T4E_6892" trichinella_intersection/79192_T4D_11040.1.faa -A 2 | awk "{print length}"


