## A script to do all the heavy liftOvering
## Want to convert Miotto's v2 coordinates to our v3 coordinates
## Started 19 November 2015

# 1. liftOver the entire Miotto bed
# 2. Split the VCF into CP group VCFs
# 3. Intersect liftedOver bed with by-group VCFs to produce by-group liftOvers
# 4. Intersect liftedOver bed with by-group VCFs to produce intersected liftover bed
# 5. Read into R, sum the allele frequencies, and do the comparisons

# 1. liftOver the entire Miotto bed
./liftOver beds/miotto.bed chainfile/2to3.liftOver beds/liftover.bed beds/liftover.trash.bed

# 2. Split the VCF into CP group VCFs
vcftools --vcf vcfs/our_goods_UG.pass.vcf --keep cp_groups/cp1.txt --recode --out vcfs/our_goods_cp1 # split out cp1
vcftools --vcf vcfs/our_goods_UG.pass.vcf --keep cp_groups/cp2+outliers.txt --recode --out vcfs/our_goods_cp2+outliers # split out cp2
vcftools --vcf vcfs/our_goods_UG.pass.vcf --keep cp_groups/cp3.txt --recode --out vcfs/our_goods_cp3 # split out cp3
vcftools --vcf vcfs/our_goods_UG.pass.vcf --keep cp_groups/cp4.txt --recode --out vcfs/our_goods_cp4 # split out cp4

# 3. Intersect liftedOver bed with by-group VCFs to produce intersected VCFs
bedtools intersect -a vcfs/our_goods_cp1.recode.vcf -b beds/liftover.bed > vcfs/cp1_intersected.vcf # itersect cp1 vcf
bedtools intersect -a vcfs/our_goods_cp2+outliers.recode.vcf -b beds/liftover.bed > vcfs/cp2_intersected.vcf # itersect cp2 vcf
bedtools intersect -a vcfs/our_goods_cp3.recode.vcf -b beds/liftover.bed > vcfs/cp3_intersected.vcf # itersect cp3 vcf
bedtools intersect -a vcfs/our_goods_cp4.recode.vcf -b beds/liftover.bed > vcfs/cp4_intersected.vcf # itersect cp4 vcf

# 4. Intersect liftedOver bed with by-group VCFs to produce intersected liftover bed
bedtools intersect -a beds/liftover_metadata.bed -b vcfs/our_goods_UG.pass.vcf > beds/liftover_intersected.bed

# 5. Read into R, sum the allele frequencies, and do the comparisons
