## IN PROGRESS 8/7/2022

## Sections:

> [Assembly](#assembly)
> 
> [Extracting target genome from metagenome](#extracting-target-genome-from-metagenome-with-an-ESOM)
> 
> [Gene Prediction](#gene-prediction-with-Funannotate)
> 
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

3. Multiple hits will be assigned to each contig, filter the results so only a "best" blast hit is assigned. Do this with [best_blast_hit.pl](/best_blast_hit.pl).

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

10. Transfer files to your local computer. Load your .names, .umx, .wts, and fixed .bm file into your [Databionic ESOM Tools](http://databionic-esom.sourceforge.net/) GUI. If you encounter issues with this program, see extra tips [here](/ESOM_installation_tips).


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
3. Run gene prediction pipeline.
```
funannotate predict -i masked_sorted_assembly.fas --name IDENTIFIER_ --species "Species name" --protein_evidence /path/to/protein/evidence -o species_funannotate_out --cpus $SLURM_NTASKS 
```
***
## This workflow uses nucleotide data (CDS files) as input.
## Phylogenomics Workflow


1. Concatenate your original cds files and convert to a one line fasta
```
cat *.fasta | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' > allcds_dna.fa
```
2. Make a new directory and move your concatenated file
```
mkdir cds_genes && mv allcds_dna.fa cds_genes
```
3. Translate sequences to aa with EMBOSS

```
for f in *.fasta
do
~/EMBOSS-6.6.0/emboss/transeq -trim -sequence ${f} -out ${f%.fasta}.preaa.fasta
done
```

4. Run protein ortho, 

Run [ProteinOrtho](http://www.bioinf.uni-leipzig.de/Software/proteinortho/)
```
proteinortho -clean -project=project_name -cpus=$SLURM_NTASKS *.aa.fasta
```
5. Get single copy, shared orthologous genes. This examples has 114 input genomes.
Optional: Check for a blank line up top. Remove if necessary with the sed command before grep. 
```
sed -n '1p' myproject.proteinortho.tsv > single_copy_114.tsv
grep $'^114\t114' myproject.proteinortho.tsv > single_copy_114.tsv
```
6. Grab proteins from files with [src_proteinortho_grab_proteins.pl](src_proteinortho_grab_proteins.pl). It's a good idea to take the header from the proteinortho output and include it at the top of your new tsv file. This can speed up the process. 
```
perl ./src_proteinortho_grab_proteins.pl -exact -tofiles new_all_114.tsv input_protein_file_1.fas input_protein_file_2.fas input_protein_file_3.fas
```
7.  Copy the new single copy files to cds_genes then cd to that directory
```
mv singlecopy114.tsv.OrthoGroup* cds_genes && cd $_
```

8. Align the single copy files with MAFFT
```
for f in *fasta
do
mafft --thread 12 --amino --maxiterate 1000 ${f} > ${f%.fasta}.aln
done

```
9. Might as well move the old stuff out of the way
```
mkdir protsinglecopy && mv *.fasta $_
```

10. Using your allcds_dna.fa, back translate to aa using RevTrans. This uses the match name to back translate. RevTrans is easy with a conda installation, but there's also a webserver.
```
for f in *aln
do 
revtrans.py -match name allcds_dna.fa ${f} > ${f%.aln}.revtrans.fasta
done
```


11. Run [trimAl](http://trimal.cgenomics.org/) with gappyout parameter
```
for f in *.revtrans.fasta
do 
trimal -in ${f} -out ${f%.revtrans}.trim
done

```

Note: If you encounter any problems from here, it’s because of the uneven sequence count making it “not an alignment”. All sequences must have the same number of bases to be considered an alignment. You can check a single file using this command:
```
bioawk -c fastx '{ print $name, length($seq) }' inputfilename
```
### Creating an ML tree

1. Concatenate your fasta files. Be sure to check your headers here before moving forward

```
catfasta2phyml.pl -f *trim > sequence.fasta
```
2. Run [RAxML-NG](https://github.com/amkozlov/raxml-ng)
```
raxml-ng --all sequence.fasta --model GTR+G --bs-trees 300

```
## Creating a gene trees and running ASTRAL

***This process uses the files in trimmed/ to prepare individual gene trees to then estimate a species tree with ASTRAL***

1. Run concatenated alignments of genes through RAxML to prepare individual gene trees
```
for f in *.trim
do
echo “raxml-ng --all ${f} --model GTR+G --bs-trees 100”
done

```
2. Concatenate the resulting gene trees and clean up our folder a bit
```
mkdir genetrees && cd $_
mv ../*.support .
cat *.support > Astral_in.tre
```

3. Optional: Filter gene trees with low support branches using Newick utilities. Example filters anything below 30%. See ASTRAL paper for more information.
```
nw_ed Astral_in.tre 'i & b<=30' o > Astral_in30.tre
```

4. Run [ASTRAL](https://github.com/smirarab/ASTRAL). Output can be visualized as a species coalescent tree. 
```
java -jar /projects/chjo1591/Astral/astral.5.7.8.jar -i Astral_in30.tre -o Astral_out30.tre
```
## Performing and visualizing condordance analysis
***Your ASTRAL output and .support files will be needed for this


1. Astral chooses an arbitrary taxon to root your tree to. It’s recommended that you manually reroot to your outgroup instead.

Example outgroup file:
```
Saccharomyces cerevisiae
```
Reroot your tree
```
reroot_trees.py Astral_30out.tre outgroup.list > rerooted_Astral30.tre

```
2. Run the Phyparts concordance analysis
```
java -jar phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -v -d ../support -m rerooted_Astral30.tre -o out

### -a = the type of analysis
### -v = verbose
### -m = target "species tree" from ASTRAL
### -d = directory with gene trees
```

3. Visualize the results with phypartspiecharts. The number is the orthologous clusters/gene trees. 

```
phypartspiecharts.py reroot_astral.tre out 84
```




