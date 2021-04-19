#!/bin/bash

### Command to change the extension of every .fa file to .faa ###

for file in $1/*.fa;
do

mv -- "$file" "${file%.fa}.faa"

done
