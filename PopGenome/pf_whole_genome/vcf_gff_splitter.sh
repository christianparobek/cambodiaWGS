## To split VCF and GFF files and put them where they belong
## for reading into PopGenome
## Formalized 08 October 2015

##########
### SET USEFUL VARIABLES
##########

dadipath=/proj/julianog/users/ChristianP/cambodiaWGS/dadi/data/pf/cp_groups


## Optional code bit to split by CP group
module add vcftools
vcftools --vcf our_goods_UG.pass.vcf --keep $dadipath/cp1.txt --recode --out  our_goods_UG.pass.cp1
vcftools --vcf our_goods_UG.pass.vcf --keep $dadipath/cp2.txt --recode --out our_goods_UG.pass.cp2
vcftools --vcf our_goods_UG.pass.vcf --keep $dadipath/cp3.txt --recode --out our_goods_UG.pass.cp3
vcftools --vcf our_goods_UG.pass.vcf --keep $dadipath/cp4.txt --recode --out our_goods_UG.pass.cp4


## Split VCF and GFF by Chromosome

for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
do

grep "^#\|^Pf3D7_$i\_v3" /proj/julianog/users/ChristianP/cambodiaWGS/PopGenome/pf_whole_genome/our_goods_UG.pass.vcf > /proj/julianog/users/ChristianP/cambodiaWGS/PopGenome/pf_whole_genome/vcf/chr$i.vcf

done

grep "^#\|^Pf3D7_$i\_v3" /proj/julianog/refs/Pf3D7_v13.0/PlasmoDB-13.0_Pfalciparum3D7.gff > /proj/julianog/users/ChristianP/cambodiaWGS/PopGenome/pf_whole_genome/gff/chr$i.gff 

done


## Optional code bit if need to clean up the names in the VCF files
#sed -i -e '/#CHROM/ s/_20[0-9]*//g' pf_whole_genome/*


