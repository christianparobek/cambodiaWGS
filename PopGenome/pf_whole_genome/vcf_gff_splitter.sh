## To split VCF and GFF files and put them where they belong
## for reading into PopGenome
## Formalized 08 October 2015

for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
do

grep "^#\|^Pf3D7_$i\_v3" /proj/julianog/users/ChristianP/cambodiaWGS/PopGenome/pf_whole_genome/our_goods_UG.pass.vcf > /proj/julianog/users/ChristianP/cambodiaWGS/PopGenome/pf_whole_genome/vcf/chr$i.vcf

grep "^#\|^Pf3D7_$i\_v3" /proj/julianog/refs/Pf3D7_v13.0/PlasmoDB-13.0_Pfalciparum3D7.gff > /proj/julianog/users/ChristianP/cambodiaWGS/PopGenome/pf_whole_genome/gff/chr$i.gff 

done


## Optional code bit if need to clean up the names in the VCF file
#sed -i -e '/#CHROM/ s/_20[0-9]*//g' pf_whole_genome/*
