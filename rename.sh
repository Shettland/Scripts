#!/bin/bash

for file in $1/*.fa;
do

mv -- "$file" "${file%.fa}.faa"

done
