#!/bin/bash

fastadir=~/TFG/Analysis/all_parasites/get_homologues/fasta_files  #folder containing every fasta file to be used: change it at your disposal

DIR=$(dirname "$(readlink -f "$0")")
echo "Folders created in: $DIR"

mkdir $DIR/ref_fasta

ln -s $fastadir/{trichinella_pseudospiralis*,trichinella_spiralis*,trichinella_britovi*} $DIR/ref_fasta/  #trichinella reference folder

filter=()	###  cut -d " " -f 5 <- To get specie's name (example: DIR/trichinella_spiralis.PRJNA12603.WBPS15.protein.faa --> DIR (trichinella_spiralis) PRJNA12603 WBPS15 protein faa ###
filt=$(find $fastadir/*.faa | grep -vE "trichinella_pseudospiralis|trichinella_spiralis|trichinella_britovi" | tr "/" " " | tr "." " " | rev | cut -d " " -f 5 | rev | uniq > $DIR/tempnames.txt)
readarray -t filter < $DIR/tempnames.txt

rm $DIR/tempnames.txt

for name in ${filter[@]};
do

mkdir $DIR/comp-$name
mkdir $DIR/comp-$name/fasta_files

ln -s $DIR/ref_fasta/* $DIR/comp-$name/fasta_files/
ln -s $fastadir/$name* $DIR/comp-$name/fasta_files/

ln -s $DIR/comp-$name/fasta_files/ $DIR/comp-$name/fasta_filesM
ln -s $DIR/comp-$name/fasta_files/ $DIR/comp-$name/fasta_filesG

done

folders=()
search=$(find $DIR/comp-* -maxdepth 0 -type d > tempfile.txt)
readarray -t folders < tempfile.txt

rm tempfile.txt

for folder in ${folders[@]};
do

pushd $folder 	### changing working directory temporarily to each new folder

nohup ~/get_homologues/get_homologues-x86_64-20210217/get_homologues.pl -d $folder/fasta_files -n 4  &> output.log &
waiter="$!"

nohup ~/get_homologues/get_homologues-x86_64-20210217/get_homologues.pl -d $folder/fasta_filesM -M -A -P -n 4  &> outputM.log &
waiter="$waiter $!"

nohup ~/get_homologues/get_homologues-x86_64-20210217/get_homologues.pl -d $folder/fasta_filesG -G -A -P -n 4  &> outputG.log &
waiter="$waiter $!"

wait $waiter


namef=$(echo "$folder" | tr "/" " " | rev | cut -d " " -f 1 | rev)

echo ""
echo "-----------------------------------------------------------------------------------------------------"
echo "            FINISHED GET_HOMOLOGUES FOR $namef"
echo "-----------------------------------------------------------------------------------------------------"
echo ""

foldOMCL=$(find $folder/fasta_filesM_homologues/*alltaxa* -maxdepth 0 -type d)
foldBDBH=$(find $folder/fasta_files_homologues/*alltaxa* -maxdepth 0 -type d)
foldCOG=$(find $folder/fasta_filesG_homologues/*alltaxa* -maxdepth 0 -type d)

nohup ~/get_homologues/get_homologues-x86_64-20210217/compare_clusters.pl -o $folder/all_intersection -m -d \
$foldOMCL,\
$foldBDBH,\
$foldCOG \
&> compare.log &

popd 	### back to main pwd

done

