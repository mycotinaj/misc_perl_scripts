#!/bin/bash
 
files=$(ls *.out)
 
for i in $files
do
echo "trimal -in ${i} -out ${i}.trim -gappyout"
done
 
