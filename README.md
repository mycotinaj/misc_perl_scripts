## Sections:

> [Assembly](#assembly)
> 
> [Extracting target genome from metagenome](#extracting-target-genome-from-metagenome-with-an-ESOM)
> 
> [Gene Prediction](#gene-prediction-with-Funannotate)
> 
> [CAZyme Annotation](#cazyme-annotation)
> 
> [Phylogenomics workflow](#phylogenomics-workflow)


***

## Assembly

1. Trim reads with [Trimmomatic v0.39](http://www.usadellab.org/cms/?page=trimmomatic) with the following parameters
>  LEADING:3 TRAILING:15 SLIDINGWINDOW:4:10 MINLEN:36

2. Assembly with [SPAdes v3.11.1](http://cab.spbu.ru/files/release3.11.1/) using the -meta flag.

***

## Extracting target genome from metagenome with an ESOM

Modified and emended from the work out of Jill Banfield's lab. Greg Dick et al. 2009 Genome Biology. And Vincent Denef's version of the print_tetramer_freqs.pl perl script has been used.  

1. Remove small contigs from your input assembly (scaffolds.fasta from SPAdes output) with [remove_small_contigs.pl](/remove_small_contigs.pl). Contigs below 300 - 1000bp can be removed. We start with 300 and increase the minimum threshold size as needed.

```
perl remove_small_contigs.pl 300 scaffolds.fasta > species_300bp_spades.fasta
```

2. Using blastn create a best blast hit for each contig. Use tabular format (6) for blast output.

```
blastn -query species_300bp_spades.fasta -db /your_db_path -out species_300bp_to_ntdb_E5 -max_target_seqs 1 -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" -num_threads 12
```

Your best blastn hit could be spurious or as equally likely as the second or fourth best blastn hit, but in aggregate, we try to separate the data. For bacteria, this actually works quite well, for Eukaryotes the databases are sparser and therefore we give a best blastn hit less meaning. A top hit to ANY fungus is usually a definite for our targets, a top hit to animals (especially inverts) could still be a fungal target.

3. Multiple hits will be assigned to each contig, filter the results so only a "best" blast hit is assigned.

```
perl best_blast_hit.pl species_300bp_to_ntdb_E5_blast.out > best_species_300bp_to_ntdb_E5_blast.out 
```

4. Assign taxonomy based on the best blast hit with [new_ncbi_tax_report.py](/new_ncbi_tax_report.py).

```
python new_ncbi_tax_report_test.py best_species_300bp_to_ntdb_E5_blast.out -o best_species_300bp_to_ntdb_E5_taxreport
```

5. Using the taxonomy assignment, bin contigs manually. Reformat to the "annotation" format taken by ESOM and save it as a plain text file: species_300bp_metagenome.annotation

Here is an example of the top few lines of the annotation file, species_300bp_metagenome.annotation:
```
contig  annotation
NODE_1147_length_2667_cov_16.2489	10
NODE_707_length_4433_cov_14.0962	10
NODE_329_length_11939_cov_66.1998	1
NODE_326_length_12118_cov_62.7537	1
NODE_2533_length_1171_cov_47.2948	1
NODE_1903_length_1542_cov_12.076	10
NODE_1518_length_1955_cov_17.2579	10
NODE_989_length_3135_cov_14.3055	10
NODE_1005_length_3078_cov_17.2253	10
NODE_988_length_3137_cov_10.7748	10
NODE_1225_length_2482_cov_16.2073	10
NODE_867_length_3613_cov_15.5486	10
NODE_1675_length_1764_cov_10.7923	10
NODE_458_length_7907_cov_59.2615	3
```
Note, the "key" to these annotation (e.g.taxonomy for 1,10,3) is kept in the spreadsheet where this work was done - not in the esom input file.

Example of the key for annotation:
```
0 Unassigned
1 Ascomycota
2 Other fungi
3 Plant
4 Animal
5 Bacteria 
6 Virus
7 Amoeba
10 Other eukaryote
```
6. Generate input files for ESOM with [print_tetramer_freqs.pl](/print_tetramer_freqs.pl).

```
perl print_tetramer_freqs_esom.pl -s species_300bp_spades.fasta -m 300 -w 3000 -k 4 -a species_300bp_metagenome.annotation

#s = your assembly
#m = minimum contig size to be analyzed (make this higher for faster/ less CPU intensive visualization)
#w = window size to calculate kmer frequency (again, make this higher for faster viz)
#k = kmer size, default is 4, longer means slower!
#a = annotation file (tab delimited plain text) you generated (either from Blast taxonomy report as above, or other info e.g. coverage, RNASeq data)
```
7. Output, 4 files: .lrn, .names, .cls, .annotation files within your directory. 

***You can run the next step locally with ESOM, however, [Somoclu](https://somoclu.readthedocs.io/en/stable/download.html) is advantageous when you have a large dataset as this can be time and resource-consuming when running locally.***


8. Train the ESOM

```

somoclu -e 20 -l 0.5 -L 0.1 -m toroid -r 50 -x 450 -y 500 -v 2 species_300bp.fasta.lrn species_300bp_spades.fasta

#e = number of epochs
#l = starting learning rate 
#L = finishing learning rate
#m = map type
#r = start radius
#x = number of columns in map 
#y = number of rows in map 
#v = verbosity level (gives progress/status in output file)
```

***Somoclu will number the the .lrn file starting at 1 and the .bm file starting at 0. If you use them as is in ESOM, you will get a message about mismatches in observations when you run the "Getclassfasta" script on your generated .cls file. Renumber the .bm file starting with 1 before importing into ESOM.***

9. Renumbering is easy with awk
```
awk 'NR>2 {$1=NR-2}{print}' species_300bp.fasta.bm > renum.species_300bp.fasta.bm
```

10. Transfer files to your local computer. Load your .names, .umx, .wts, and fixed .bm file into your [Databionic ESOM Tools](http://databionic-esom.sourceforge.net/) GUI. 


11. Vizualize output locally and identify target genome

View -> UMatrix background, tiled display.
Use Zoom, Color, Bestmatch size to get your desired view.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A. Check "Draw Best Matches" in View pane.

  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B. File -> load .cls
 
   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;C. Under "Classes" tab, select classes to be displayed on the map.
 
 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;D. Once you have identified your target genome, select the "Data" tab. Then, left click around the target area of the map. 
 
  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;E. When you have enclosed all of the target area of the map, right click to select that target area. 

 ***This will only work correctly if you have selected the "Data" tab prior to beginning.***
 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;F. File, Selection, save as .cls. Name your selected .cls file to save to your folder. 
 
12. Obtain target contigs using [Getclassfasta.pl](/Getclassfasta.pl)
```
perl Getclassfasta.pl -fasta species_300bp_spades.fasta -names species_300bp_spades.fasta.names -num 1 -loyal 75 -cls species_300bp_spades.fasta.cls
```
This will generate a .fasta (your target genome assembly), and a .conf giving confidence values of specific contigs

For more in-depth steps and example files, see [Quandt Mycology's Running_ESOM](https://github.com/Quandt-Mycology-Lab/Lab_Codes_and_workflows/tree/13a29e99deba6a064015bf7791746de695150d0f/Running_ESOM) page.

***

## Gene Prediction with [Funannotate](https://funannotate.readthedocs.io/en/latest/)

1. Sort by contigs by size and rename headers
```
funannotate sort -i after_esom_assembly.fasta -o sorted_assembly.fasta -b profar
```
2. Mask repeats
```
funannotate mask -i sorted_assembly.fasta -o masked_sorted_assembly.fas -s fungi
```
3. Run gene prediction pipeline. ***Protein evidence used for this study***
```
funannotate predict -i masked_sorted_assembly.fas --name IDENTIFIER_ --species "Species name" --augustus_species botrytis_cinerea --protein_evidence /path/to/protein/evidence -o species_funannotate_out --cpus $SLURM_NTASKS 
```
***

### CAZyme Annotation

This will identify, filter, and count CAZyme families and subfamilies within a protein fasta with [dbCAN](http://bcb.unl.edu/dbCAN2/download/) and merge results of multiple genomes into a single CSV. This will find cazymes with coverage > 0.45 and e-value <1e-17 from proteins.fasta file using hmmer - these parameters are suggested for fungi.

1. Run dbCAN on .proteins.fa file generated from funannotate.
```
run_dbcan.py --out_pre yourspecies --hmm_cov 0.45 --hmm_eval 1e-17 --hmm_cpu 12 --db_dir /your/dbcan/database yourinputfile.proteins.fa protein
```
2. Remove header line. Use wc -l after this and other steps if you want to check if they worked
```
tail -n +2 hmmer.out > hmmer_step1.out
```

3. Remove all columns except cazy name and gene
```
awk '{ print $1, $3}' hmmer_step1.out > hmmer_step2.out
```

4. Remove any duplicate annotations on the same gene
```
awk '!x[$0]++' hmmer_step2.out > hmmer_step3.out
```

5. Generate counts of cazys
```
awk '{print $1}' hmmer_step3.out | sort | uniq -c > hmmer_step4.out
```

6. Swap columns
```
awk ' { t = $1; $1 = $2; $2 = t; print; } ' hmmer_step4.out > hmmer_step5.out
```

7. Add species name as column header (format as G.species no space).
```
awk 'BEGIN{print "Species\tG.species"}1' hmmer_step5.out > yourspecies_cazy_counts.out
```
Once you have all of your yourspecies_cazy_counts.out files, put them into a single folder. The next steps will merge the results into one file. [csvtk](https://bioinf.shenwei.me/csvtk/download/) needs to be installed.

1. Convert *.out files to tab-delimited
```
ls *.out | parallel 'csvtk space2tab {} > {.}.tsv'
```
2. Grab all possible values from the first column
```
cut -f 1 *.tsv | sort | uniq | tee all_cazys_step1.tsv
```

3. Join files and fill blanks with 0s
```
csvtk join -H -t -k -t *.tsv --na 0 > all_cazys_counts_step2.tsv
```
4. Clean .hmm from final result
```
sed 's/.hmm//' all_cazys_counts_step2.tsv > all_cazys_counts_step3.tsv
```

5. Flip so species names are headers
```
tac all_cazys_counts_step3.tsv > all_cazys_counts_final.tsv 
```

6. Convert to csv - possibly better for R?
```
csvtk tab2csv all_cazys_counts_final.tsv > all_cazys_counts_final.csv
```
```
Other commands that might help:

#Filter by evalue

awk '$5 < 1e-17 { print $0 }' hmmer.out

#Filter by coverage 
awk '$9 > 0.45 { print $0 }' hmmer.out 

#Sum of your cazy counts
awk '!x[$0]++' hmmer_step2.out
```
***

## Phylogenomics Workflow

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
6. Grab proteins from files with [src_proteinortho_grab_proteins.pl](src_proteinortho_grab_proteins.pl)
```
perl ./src_proteinortho_grab_proteins.pl -exact -tofiles new_all_18.tsv input_protein_file_1.fas input_protein_file_2.fas input_protein_file_3.fas
```
7. Make a new directory and move single copy files
```
mkdir single_copy
mv *.fasta single_copy
cd single copy
```
8. Create an lb file from files in this directory for MUSCLE input using [new_create_muscle.sh](/new_create_muscle.sh)
```
./new_create_muscle.sh > lb_cmd_file
```
9. Run [MUSCLE](https://www.drive5.com/muscle/downloads.htm) with [Open MPI](https://www.open-mpi.org/software/ompi/v4.1/). 
```
mpirun lb_cmb_file
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
mpirun lb lb_cmd_raxml
```
4. Concatenate all gene trees
```
cat raxml.bootstraps > all_gene_trees.tree
```
5. Run [ASTRAL](https://github.com/smirarab/ASTRAL)
```
java -jar astral.5.7.3.jar -i all_gene_trees.tree -o output.tree
```
