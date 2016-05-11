## Started 11 May 2016
## Christian Parobek
## To subset my P. falciparum and P. vivax gffs 
## Prior to splitting and use in PopGenome

pf_gff=/proj/julianog/refs/Pf3D7_v13.0/PlasmoDB-13.0_Pfalciparum3D7.gff
pv_gff=/proj/julianog/refs/PvSAL1_v13.0/PlasmoDB-13.0_PvivaxSal1.gff
nick=/proj/julianog/users/NickB/Projects/cambodiaWGS/data\&output


################################
######### For P. vivax #########
################################

grep "Pv_Sal1_chr\|##" $pv_gff | grep -v "AAKM\|NC\|FASTA\|^>" > PvSal1_features.gff
	## First, make a GFF with only gene entries and header
	## Excluding entries on AAKM contigs, mito, and the random "FASTA" line


bedtools intersect -v -a PvSal1_features.gff -b $nick/subtelo_60k_altered.bed |\
 bedtools intersect -v -a "stdin" -b $nick/neafseyExclude_altered.bed > PvSal1_filtered.gff
	## Then remove subtelomeres and Neafsey genes
	## Then make sure to paste back in the header

## Then need to split GFF by chromosome
for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
do

grep "^#\|Pv_Sal1_chr$i" PvSal1_filtered.gff > gff/chr$i.gff

done

#################################
####### For P. falciparum #######
#################################

grep "Pf3D7_\|##" $pf_gff | grep -v "IRAB\|M76611\|FASTA" > Pf3D7_features.gff
	## First, make a GFF with only gene entries and header
	## Excluding entries on AAKM contigs, mito, and the random "FASTA" line

bedtools intersect -v -a Pf3D7_features.gff -b $nick/subtelomeres_pf3d7.bed |\
 bedtools intersect -v -a "stdin" -b $nick/VAR_STEVOR_RIFIN_pf3d7.bed > Pf3D7_filtered.gff
	## Then remove subtelomeres and Neafsey genes
	## Then make sure to paste back in the header

## Then need to split GFF by chromosome
for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
do

grep "^#\|Pf3D7_$i\_v3" Pf3D7_filtered.gff > gff/chr$i.gff

done

## Then replace the original split GFF files in popgenome folder with these split gffs







