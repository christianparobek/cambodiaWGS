## A snakefile for snakemake
## Started 28 January 2016
## Want to make sure all figures are based on most up-to-date datasets


workdir: '/proj/julianog/users/ChristianP/cambodiaWGS/'


####### Target #######
rule all:
	#input: 'Fws/fws.svg'
	#input: 'adegenet/grey_pca.svg'
	input: 'selscan/pv_all.svg'
	#input: 'selscan/pf_all.svg'
	#input: 'selscan/pv_mono.svg'
	#input: 'selscan/pf_cp2.svg'


rule pf_nsl_cp2 :
	input: vcf = 'pf/variants/cp2.pass.recode.vcf'
	output: 'selscan/pf_cp2.svg'
	shell: 'bash selscan/nslRunner.sh {input.vcf} pf_cp2 pf'

rule pf_nsl_all :
	input: vcf = 'pf/variants/our_goods_UG.pass.vcf'
	output: 'selscan/pf_all.svg'
	shell: 'bash selscan/nslRunner.sh {input.vcf} pf_all pf'

rule pv_nsl_mono :
	input: vcf = 'pv/variants/mono.pass.recode.vcf'
	output: 'selscan/pv_mono.svg'
	shell: 'bash selscan/nslRunner.sh {input.vcf} pv_mono pv'

rule pv_nsl_all :
	input: vcf = 'pv/variants/our_goods_UG.pass.vcf'
	output: 'selscan/pv_all.svg'
	shell: 'bash selscan/nslRunner.sh {input.vcf} pv_all pv'

rule adegenet_pca :
	input: pv = 'adegenet/our_goods_pv.pass.sansAAKM.vcf', pf = 'adegenet/our_goods_pf.pass.vcf'
	output: grey = 'adegenet/grey_pca.svg', color = 'adegenet/color_pca.svg'
	shell: 'Rscript adegenet/adegenet.r {input.pv} {input.pf} {output.grey} {output.color}'

rule calc_fws_boot :
	input: pv = 'pv/variants/our_goods_UG.pass.vcf', pf = 'pf/variants/our_goods_UG.pass.vcf'
	output: 'Fws/plot.svg'
	shell: 'Rscript Fws/fws.r {input.pv} {input.pf} {output}'

