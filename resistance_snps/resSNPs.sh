## Script to identify and tally SNPs in resistance genes
## For both P. falciparum and P. vivax
## Started 2 May 2016

cp2vcf=/proj/julianog/users/ChristianP/cambodiaWGS/pf/variants/cp2.for_res_snp_calling_only.recode.vcf
monovcf=/proj/julianog/users/ChristianP/cambodiaWGS/pv/variants/mono.for_res_snp_calling_only.recode.vcf



#########################
## Pv RESISTANCE GENES ##
#########################

bedtools intersect -a $monovcf -b beds/pv_res_genes.bed > mono_res_genes.vcf

java -Xmx4g -jar ~/snpEff/snpEff.jar -v -o gatk PvSAL1v13 mono_res_genes.vcf > mono_res_genes.anno.vcf

Rscript vcf_digester.r mono_res_genes.anno.vcf mono_res_table.txt


#########################
### Pf RES GENES CP2 ####
#########################

bedtools intersect -a $cp2vcf -b beds/pf_res_genes.bed > cp2_res_genes.vcf

java -Xmx4g -jar ~/snpEff/snpEff.jar -v -o gatk Pf3D7v13 cp2_res_genes.vcf > cp2_res_genes.anno.vcf

Rscript vcf_digester.r cp2_res_genes.anno.vcf cp2_res_table.txt





#########################
#### Pv nSL TOP HITS ####
#########################

bedtools intersect -a $monovcf -b beds/pv_nsl_top_hits.bed > pv_nsl_genes.vcf

java -Xmx4g -jar ~/snpEff/snpEff.jar -v -o gatk PvSAL1v13 pv_nsl_genes.vcf > pv_nsl_genes.anno.vcf

Rscript vcf_digester.r pv_nsl_genes.anno.vcf pv_nsl_table.txt


