#!/bin/bash
 
files=$(ls *.trim.renamed)
 
for i in $files
do
echo "raxml-ng-mpi --all -msa ${i} --prefix ${i}_raxml_v1_100bp --model WAG --seed 2 --bs-trees 100"

done
