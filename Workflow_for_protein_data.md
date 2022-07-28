## Phylogenomics Workflow For Using Protein Data

### This is an older version of the workflow which uses protein data as input. The Phyparts concordance analysis is also not included. I recommend using nucleotide data and the respective workflow instead. 

1. Change headers in files

For JGI files:
```
cat Marbr1_GeneCatalog_proteins_20141203.aa.fasta | sed 's/jgi|//g' | sed 's/|/_/g' | sed 's/Marbr1_/Marbr1|/g' > Marbr1.names.fas
```
For Funannotate proteins.fas:
```
cat Nagfri.proteins.fas | sed 's/>NAGFRI_/>NAGFRI|g/g' > Nagfri.proteins_rename.fas
```
For NCBI-Genbank files. POR is the 3 letter accession.
```
cat Tolpar.fas | sed 's/ \[Tolypocladium paradoxum\]//g' | sed 's/>POR/>Tolpar|POR/g' | sed 's/ /_/g' | sed 's/,//g' | sed 's/\///g' | sed 's/(//g' | sed 's/)//g' | sed 's/[][]//g' > Tolpar.names.fas
```
2. Make a new directory and move all renamed files
```
mkdir good_names
mv *.names.fas good_names
cd good_names
```
3. Get the names of all files in the directory for the next step
```
echo $(dir) 
```
4. Run [ProteinOrtho](http://www.bioinf.uni-leipzig.de/Software/proteinortho/)
```
proteinortho -clean -project=project_name -cpus=$SLURM_NTASKS input_protein_file_1.fas input_protein_file_2.fas input_protein_file_3.fas
```
5. Get single copy, shared orthologous clusters. This example has 3 input proteomes:
```
grep $'^3\t3' protortho_output.tsv > new_all_3.tsv
```
6. Grab proteins from files with [src_proteinortho_grab_proteins.pl](src_proteinortho_grab_proteins.pl). It's a good idea to take the header from the proteinortho output and include it at the top of your new tsv file. This can speed up the process. 
```
perl ./src_proteinortho_grab_proteins.pl -exact -tofiles new_all_18.tsv input_protein_file_1.fas input_protein_file_2.fas input_protein_file_3.fas
```
7. Make a new directory and move single copy files
```
mkdir single_copy
mv *.fasta single_copy
cd single_copy
```
8. Create an lb file from files in this directory for MUSCLE input using [new_create_muscle.sh](/new_create_muscle.sh)
```
./new_create_muscle.sh > lb_cmd_file
```
9. Run [MUSCLE](https://www.drive5.com/muscle/downloads.htm) with [Open MPI](https://www.open-mpi.org/software/ompi/v4.1/). 
```
mpirun lb lb_cmb_file
```
10. Move output files to new directory
```
mkdir muscle_out
mv *.fasta.out muscle_out
cd muscle_out
```
11. Create an lb file from files in this directory for trimAl input using [create_trimal.sh](/create_trimal.sh)
```
./create_trimal.sh > lb_cmd_trimal
```
12. Run [trimAl](http://trimal.cgenomics.org/)
```
mpirun lb lb_cmd_trimal
```
13. Make a new folder and move trimmed files to this directory
```
mkdir trimmed
mv *.out.trim trimmed
```
***

### Running [RAxML-NG](https://github.com/amkozlov/raxml-ng)
1. Create an expected file from your trimmed alignments. Check this output file for the presence of all of your taxa
```
cat trimmed/*.fasta.out.trim | grep "^>" | tr "|" "\t" | cut -f 1 | sort | uniq > expected.txt
```
2. Concatenate your alignment with [combine_fasaln.pl](/combine_fasaln.pl) and your expected file
```
perl ./combine_fasaln_v2019.pl -o allseqs.fas -of fasta -d ./trimmed/ --expected expected.txt
```
3. Make a new directory and move your alignment
```
mkdir raxmlng
mv allseqs.fas raxmlng
cd raxmlng
```
4. Run [RAxML-NG](https://github.com/amkozlov/raxml-ng) with WAG model (formerly PROTGAMMAWAG), 200 bs replicates and a seed value of 2 (this is arbitrary and can be changed)
```
raxml-ng-mpi --all -msa allseqs.fas --model WAG --opt-model on --seed 2 --bs-trees 200
```

### Running [ASTRAL](https://github.com/smirarab/ASTRAL)
***This process uses the files in trimmed/ to prepare individual gene trees to then estimate a species tree with ASTRAL***

1. Rename fasta headers across alignments in trimmed/
```
for filename in *.trim
do
    cat "${filename}" | awk -F "|" '{print $1}' > "${filename}.renamed"
done
```
2. Create an lb file using [create_astral.sh](/create_astral.sh)
```
./create_astral.sh > lb_cmd_astral
```
3. Run concatenated alignments of genes through RAxML to prepare individual gene trees
```
mpirun lb lb_cmd_astral
```
4. Concatenate all gene trees
```
cat raxml.bipartitions > all_gene_trees.tre
```
5. Optional: Filter gene trees with low support branches using Newick utilities. Example filters anything below 50%. 
```
all_gene_trees.tre 'i & b<=50' o > all_gene_tree_BS50.tre
```
6. Run [ASTRAL](https://github.com/smirarab/ASTRAL)
```
java -jar astral.5.7.3.jar -i all_gene_trees.tre -o output.tre
```
