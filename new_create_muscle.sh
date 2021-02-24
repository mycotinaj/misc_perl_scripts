#!/bin/bash
 
files=$(ls *.fasta)
 
for i in $files
do
echo "muscle -in ${i} -out ${i}.out"
done
 
