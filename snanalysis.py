## A snakefile for snakemake
## Started 28 January 2016
## Want to make sure all figures are based on most up-to-date datasets


workdir: '/proj/julianog/users/ChristianP/cambodiaWGS/'


####### Target #######
rule all:
	input: 'Fws/fws.svg'



rule calc_fws_boot :
	input: pv = 'pv/variants/our_goods_UG.pass.vcf', pf = 'pf/variants/our_goods_UG.pass.vcf'
	output: 'Fws/fws.svg'
	shell: 'Rscript Fws/fws.r {input.pv} {input.pf} {output}'

# bash Fws/fwsRunner.sh


#Rscript Fws/fws.r pv/variants/our_goods_UG.pass.vcf pf/variants/our_goods_UG.pass.vcf Fws/fws.svg


