#!/bin/bash

## SETTING HELP FUNCTION ##

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

#if [[ $2 == "PE" ]]
#then


BLUE='\033[1;34m'
NC='\033[0m'
YELLOW='\033[1;33m'

folders=()
names=()

folder=$(find $path | grep -F ".fastq" | sort | tr "/" " " | rev | cut -d " " -f 2-100 | rev | tr " " "/" | uniq > filetemp.txt)
readarray -t folders < filetemp.txt

echo " ${folders[@]}"
rm filetemp.txt

filter=$(find $path | grep -F ".fastq" | sort | tr "/" " " | rev | cut -d " " -f 1 | rev  > filtertemp.txt)
#echo "--temp 1--"
#cat filtertemp.txt

filter2=$(cat filtertemp.txt | tr "_" " " | rev | cut -d " " -f 2-100 | rev | tr " " "_" | uniq > filtertemp2.txt)
readarray -t names < filtertemp2.txt
#echo "--temp 2--"
#cat filtertemp2.txt

rm filtertemp.txt
rm filtertemp2.txt

echo -e "${BLUE} "
echo "Selected samples for trimming: ${names[@]}"
echo -e "${NC} "


c=1
p=1

## Setting variables to each sample's name ##
for m in ${names[@]}
do

opt[c]=${m}
((c++))
done

: '
## COMPRESS FILES WITH GZIP BEFORE SETTING FILENAMES ##
for folder in ${folders[@]}
do

gzip -r $folder

done
'

## SETTING ARRAY TO CONTAIN PATH+FILENAMES FOR TRIMMOMATIC ##
files1=()
files2=()

search=$(find $path | grep -F "_1.fastq.gz" | sort | uniq > file1temp.txt)
readarray -t files1 < file1temp.txt

rm file1temp.txt
#echo "${files1[@]}"

search2=$(find $path | grep -F "_2.fastq.gz" | sort | uniq > file2temp.txt)
readarray -t files2 < file2temp.txt

rm file2temp.txt
#echo "${files2[@]}"


## CREATING A SYMBOLIC LINK FOR EACH SAMPLE ##

p=1

mkdir ./symbolic_links

for file in ${files1[@]}
do

result1=$(ln -s $file ./symbolic_links/${opt[p]}_R1.fastq.gz)
((p++))
done

p=1
for file in ${files2[@]}
do

result2=$(ln -s $file ./symbolic_links/${opt[p]}_R2.fastq.gz)
((p++))

done

echo -e "${YELLOW} FINISHED CREATING SYMBOLIC LINKS TO EACH SAMPLE FILE ${NC}"

: '
p=1

for file in ${files1[@]}
do

## This function removes the last 11(???...) characters of the string ##
result=${file%???????????}

echo "trimmomatic input: ${result1} ${result2}"
echo "trimmomatic output: ${opt[p]}_paired_R1...."

done
'

mkdir ./trimmomatic_results

p=1

for name in ${names[@]}
do
	## TRIMMOMATIC SCRIPT ##

	java -jar ${trimmomatic} PE -threads 4 -phred33 ./symbolic_links/${opt[p]}_R1.fastq.gz  ./symbolic_links/${opt[p]}_R2.fastq.gz \
	./trimmomatic_results/${opt[p]}_paired_R1.fastq ./trimmomatic_results/${opt[p]}_unpaired_R1.fastq ./trimmomatic_results/${opt[p]}_paired_R2.fastq ./trimmomatic_results/${opt[p]}_unpaired_R2.fastq \
	ILLUMINACLIP:${adapters}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50 2> trimmomatic_results/trimmomatic.log

	echo -e "${YELLOW} FINISHED TRIMMING SAMPLE ${opt[p]} ${NC}"
	((p++))
	#wait $!
done

echo -e "${BLUE}"
echo "----------------------------------------------"
echo "      Finished trimming for all samples       "
echo "----------------------------------------------"
echo " Trimming results are in trimmomatic_results/ "
echo "----------------------------------------------"
echo -e "${NC}"


## Running trinity ##
## LOCATION echo "$PATH" ##

## TAKING THE FULL PATH OF UNIQUE TRIMMED SAMPLES AND APPEND TO ARRAY IN CASE SAMPLES ARE ALREADY TRIMMED ##

samples=()

archivos=$(find ./trimmomatic_results | grep -v "readcount" | grep _paired_ | sort > filetemp.txt)
filtrado=$(cat filetemp.txt | tr "_" " " | cut -d " " -f 1 | uniq > filetemp2.txt)
readarray -t samples < filetemp2.txt

rm filetemp.txt
rm filetemp2.txt

## STARTING TRINITY SCRIPT  ##

p=1

echo -e "${BLUE} "
echo "SELECTED SAMPLES FOR TRINITY ASSEMBLY: ${names[@]}"
echo -e "${NC} "


for name in ${names[@]}
do

echo -e "${YELLOW}"
echo "------------------------------------------------------------------------------------------------------------------------------"
echo "-----------------------------Running trinity for sample ${opt[p]}-----------------------------------------"
echo "------------------------------------------------------------------------------------------------------------------------------"
echo -e "${NC}"

	${trinity} --seqType fq \
	--left ./trimmomatic_results/${opt[p]}_paired_R1.fastq.gz \
	--right ./trimmomatic_results/${opt[p]}_paired_R2.fastq.gz \
	--max_memory 10G --CPU 2 \
	--output ${opt[p]}_trinity

echo -e "${BLUE} "
echo "------------------------------------------------------------------------------------------------------------------------------"
echo "----------------------------Finished assembly for sample ${opt[p]}----------------------------------------"
echo "------------------------------------------------------------------------------------------------------------------------------"
echo "------------------------------Output folder: ${opt[p]}_trinity/-------------------------------------------"
echo "------------------------------------------------------------------------------------------------------------------------------"
echo -e "${NC}"

((p++))

#wait $!

done


echo "FIN"
