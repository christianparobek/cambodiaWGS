## A snakefile for snakemake
## Started 28 January 2016
## Want to make sure all figures are based on most up-to-date datasets


workdir: '/proj/julianog/users/ChristianP/cambodiaWGS/'

#######################
####### TARGETS #######
#######################
rule all:
	#input: 'Fws/fws.svg'
	#input: 'adegenet/grey_pca.svg'
	#input: 'dadi/data/pv_all/dadi_output_afs/pv_all.afs'
	#input: 'dadi/data/pv_mono/dadi_output_afs/pv_mono.afs'
	#input: 'dadi/data/pf_all/dadi_output_afs/pf_all.afs'
	#input: 'dadi/data/pf_cp2/dadi_output_afs/pf_cp2.afs'
	#input: 'dadi/data/pfall_pvall.svg'
	#input: 'dadi/modelFitting/pv_mono/pv_mono.bg.fit'
	#input: 'dadi/modelFitting/pv_all/pv_all.bg.fit'
	#input: 'dadi/modelFitting/pf_cp2/pf_cp2.bg.fit'
	#input: 'dadi/modelFitting/pf_all/pf_all.bg.fit'
	#input: 'geneticMap/pf_cp2/logfile'
## SELSCAN IHS
	#input: 'selscan/ihs/pv_mono.svg'
	#input: 'selscan/ihs/pv_all.svg'
	#input: 'selscan/ihs/pf_cp2.svg'
	#input: 'selscan/ihs/pf_all.svg'
## SELSCAN nSL
	#input: 'selscan/nsl/pv_all.svg'
	#input: 'selscan/nsl/pf_all.svg'
	#input: 'selscan/nsl/pv_mono.svg'
	#input: 'selscan/nsl/pf_cp2.svg'
## SELSCAN EHH
	#input: 'selscan/ehh/pv_chr14_797870'
	input: 'selscan/ehh/pf_chr11_1294584'


#######################
##### SELSCAN EHH #####
#######################

rule pv_mono_chr14_797870:
	input: vcf = '/proj/julianog/users/ChristianP/cambodiaWGS/selscan/nsl/pv_mono/chr14.vcf', map = '/proj/julianog/users/ChristianP/cambodiaWGS/selscan/ihs/pv_mono/14.map'
	output: 'selscan/ehh/pv_chr14_797870'
	shell: 'bash selscan/ehh/ehhRunner.sh locus1511 {input.vcf} {input.map} {output}'

rule pf_cp2_chr11_1294584:
	input: vcf = '/proj/julianog/users/ChristianP/cambodiaWGS/selscan/nsl/pf_cp2/chr11.vcf', map = '/proj/julianog/users/ChristianP/cambodiaWGS/selscan/ihs/pf_cp2/11.map'
	output: 'selscan/ehh/pf_chr11_1294584'
	shell: 'bash selscan/ehh/ehhRunner.sh locus214 {input.vcf} {input.map} {output}'


#######################
##### SELSCAN iHS #####
#######################
## RUN ONE-AT-A-TIME ##
#######################

rule iHS_pf_all:
	input: vcf = '/proj/julianog/users/ChristianP/cambodiaWGS/pf/variants/our_goods_UG.pass.vcf'
	output: 'selscan/ihs/pf_all.svg'
	shell: 'bash selscan/ihs/iHScalc.sh {input.vcf} pf_all pf 0.35'

rule iHS_pf_cp2:
	input: vcf = '/proj/julianog/users/ChristianP/cambodiaWGS/pf/variants/cp2.pass.recode.vcf'
	output: 'selscan/ihs/pf_cp2.svg'
	shell: 'bash selscan/ihs/iHScalc.sh {input.vcf} pf_cp2 pf 0.20' # guessing, but don't know that we want to minimize this cutoff

rule iHS_pv_all:
	input: vcf = '/proj/julianog/users/ChristianP/cambodiaWGS/pv/variants/our_goods_UG.pass.vcf'
	output: 'selscan/ihs/pv_all.svg'
	shell: 'bash selscan/ihs/iHScalc.sh {input.vcf} pv_all pv 0.05'

rule iHS_pv_mono:
	input: vcf = '/proj/julianog/users/ChristianP/cambodiaWGS/pv/variants/mono.pass.recode.vcf'
	output: 'selscan/ihs/pv_mono.svg'
	shell: 'bash selscan/ihs/iHScalc.sh {input.vcf} pv_mono pv 0.05'



#######################
###### AFS PLOTS ######
#######################

rule plot_afs:
	input: pf_all = 'dadi/data/pf_all/dadi_output_afs/pf_all.afs', pf_cp2 = 'dadi/data/pf_cp2/dadi_output_afs/pf_cp2.afs', pv_all = 'dadi/data/pv_all/dadi_output_afs/pv_all.afs', pv_mono = 'dadi/data/pv_mono/dadi_output_afs/pv_mono.afs'
	output: full_datasets = 'dadi/data/pfall_pvall.svg', subset_datasets = 'dadi/data/pfcp2_pvmono.svg'
	shell: 'Rscript dadi/afsPlotter.r {input.pf_all} {input.pf_cp2} {input.pv_all} {input.pv_mono} {output.full_datasets} {output.subset_datasets}'

#######################
##### DADI PF ALL #####
#######################

rule dadi_fit_1pop_models_pf_all:
	input: afs = 'dadi/data/pf_all/dadi_input/pf_all.syn.dadi'
	output: 'dadi/modelFitting/pf_all/pf_all.bg.fit'
	shell: 'bash dadi/modelFitting-1pop.sh {input.afs} dadi/modelFitting/pf_all/ pf_all pf'

rule dadi_calc_afs_pf_all:
	input: annovcf = 'dadi/data/pf_all/vcfs/pf_all.anno.vcf', syn = 'dadi/data/pf_all/dadi_input/pf_all.syn.dadi', nonsyn = 'dadi/data/pf_all/dadi_input/pf_all.nonsyn.dadi', genic = 'dadi/data/pf_all/dadi_input/pf_all.genic.dadi', intergenic = 'dadi/data/pf_all/dadi_input/pf_all.intergenic.dadi'
	output: 'dadi/data/pf_all/dadi_output_afs/pf_all.afs'
	shell: 'bash dadi/afsPrinter.sh {input.annovcf} {input.syn} {input.nonsyn} {input.genic} {input.intergenic} {output}'

rule dadi_format_pf_all:
	input: vcf = 'pf/variants/our_goods_UG.pass.vcf'
	output: syn = 'dadi/data/pf_all/dadi_input/pf_all.syn.dadi', nonsyn = 'dadi/data/pf_all/dadi_input/pf_all.nonsyn.dadi', genic = 'dadi/data/pf_all/dadi_input/pf_all.genic.dadi', intergenic = 'dadi/data/pf_all/dadi_input/pf_all.intergenic.dadi'
	shell: 'bash dadi/vcf2dadi.sh {input.vcf} pf_all pf'

#######################
##### DADI PF CP2 #####
#######################

rule dadi_fit_1pop_models_pf_cp2:
	input: afs = 'dadi/data/pf_cp2/dadi_input/pf_cp2.syn.dadi'
	output: 'dadi/modelFitting/pf_cp2/pf_cp2.bg.fit'
	shell: 'bash dadi/modelFitting-1pop.sh {input.afs} dadi/modelFitting/pf_cp2/ pf_cp2 pf'

rule dadi_calc_afs_pf_cp2:
	input: annovcf = 'dadi/data/pf_cp2/vcfs/pf_cp2.anno.vcf', syn = 'dadi/data/pf_cp2/dadi_input/pf_cp2.syn.dadi', nonsyn = 'dadi/data/pf_cp2/dadi_input/pf_cp2.nonsyn.dadi', genic = 'dadi/data/pf_cp2/dadi_input/pf_cp2.genic.dadi', intergenic = 'dadi/data/pf_cp2/dadi_input/pf_cp2.intergenic.dadi'
	output: 'dadi/data/pf_cp2/dadi_output_afs/pf_cp2.afs'
	shell: 'bash dadi/afsPrinter.sh {input.annovcf} {input.syn} {input.nonsyn} {input.genic} {input.intergenic} {output}'

rule dadi_format_pf_cp2:
	input: vcf = 'pf/variants/cp2.pass.recode.vcf'
	output: syn = 'dadi/data/pf_cp2/dadi_input/pf_cp2.syn.dadi', nonsyn = 'dadi/data/pf_cp2/dadi_input/pf_cp2.nonsyn.dadi', genic = 'dadi/data/pf_cp2/dadi_input/pf_cp2.genic.dadi', intergenic = 'dadi/data/pf_cp2/dadi_input/pf_cp2.intergenic.dadi'
	shell: 'bash dadi/vcf2dadi.sh {input.vcf} pf_cp2 pf'

#######################
##### DADI PV ALL #####
#######################

rule dadi_fit_1pop_models_pv_all:
	input: afs = 'dadi/data/pv_all/dadi_input/pv_all.syn.dadi'
	output: 'dadi/modelFitting/pv_all/pv_all.bg.fit'
	shell: 'bash dadi/modelFitting-1pop.sh {input.afs} dadi/modelFitting/pv_all/ pv_all pv'

rule dadi_calc_afs_pv_all:
	input: annovcf = 'dadi/data/pv_all/vcfs/pv_all.anno.vcf', syn = 'dadi/data/pv_all/dadi_input/pv_all.syn.dadi', nonsyn = 'dadi/data/pv_all/dadi_input/pv_all.nonsyn.dadi', genic = 'dadi/data/pv_all/dadi_input/pv_all.genic.dadi', intergenic = 'dadi/data/pv_all/dadi_input/pv_all.intergenic.dadi'
	output: 'dadi/data/pv_all/dadi_output_afs/pv_all.afs'
	shell: 'bash dadi/afsPrinter.sh {input.annovcf} {input.syn} {input.nonsyn} {input.genic} {input.intergenic} {output}'

rule dadi_format_pv_all:
	input: vcf = 'pv/variants/our_goods_UG.pass.vcf'
	output: syn = 'dadi/data/pv_all/dadi_input/pv_all.syn.dadi', nonsyn = 'dadi/data/pv_all/dadi_input/pv_all.nonsyn.dadi', genic = 'dadi/data/pv_all/dadi_input/pv_all.genic.dadi', intergenic = 'dadi/data/pv_all/dadi_input/pv_all.intergenic.dadi'
	shell: 'bash dadi/vcf2dadi.sh {input.vcf} pv_all pv'

#######################
##### DADI PV MONO ####
#######################

rule dadi_fit_1pop_models_pv_mono:
	input: afs = 'dadi/data/pv_mono/dadi_input/pv_mono.syn.dadi'
	output: 'dadi/modelFitting/pv_mono/pv_mono.bg.fit'
	shell: 'bash dadi/modelFitting-1pop.sh {input.afs} dadi/modelFitting/pv_mono/ pv_mono pv'

rule dadi_calc_afs_pv_mono:
	input: annovcf = 'dadi/data/pv_mono/vcfs/pv_mono.anno.vcf', syn = 'dadi/data/pv_mono/dadi_input/pv_mono.syn.dadi', nonsyn = 'dadi/data/pv_mono/dadi_input/pv_mono.nonsyn.dadi', genic = 'dadi/data/pv_mono/dadi_input/pv_mono.genic.dadi', intergenic = 'dadi/data/pv_mono/dadi_input/pv_mono.intergenic.dadi'
	output: 'dadi/data/pv_mono/dadi_output_afs/pv_mono.afs'
	shell: 'bash dadi/afsPrinter.sh {input.annovcf} {input.syn} {input.nonsyn} {input.genic} {input.intergenic} {output}'

rule dadi_format_pv_mono:
	input: vcf = 'pv/variants/mono.pass.recode.vcf'
	output: syn = 'dadi/data/pv_mono/dadi_input/pv_mono.syn.dadi', nonsyn = 'dadi/data/pv_mono/dadi_input/pv_mono.nonsyn.dadi', genic = 'dadi/data/pv_mono/dadi_input/pv_mono.genic.dadi', intergenic = 'dadi/data/pv_mono/dadi_input/pv_mono.intergenic.dadi'
	shell: 'bash dadi/vcf2dadi.sh {input.vcf} pv_mono pv'

#######################
##### SELSCAN NSL #####
#######################
## RUN ONE-AT-A-TIME ##
#######################

rule nsl_pf_cp2:
	input: vcf = 'pf/variants/cp2.pass.recode.vcf'
	output: 'selscan/nsl/pf_cp2.svg'
	shell: 'bash selscan/nsl/nslRunner.sh {input.vcf} pf_cp2 pf'

rule nsl_pf_all:
	input: vcf = 'pf/variants/our_goods_UG.pass.vcf'
	output: 'selscan/nsl/pf_all.svg'
	shell: 'bash selscan/nsl/nslRunner.sh {input.vcf} pf_all pf'

rule nsl_pv_mono:
	input: vcf = 'pv/variants/mono.pass.recode.vcf'
	output: 'selscan/nsl/pv_mono.svg'
	shell: 'bash selscan/nsl/nslRunner.sh {input.vcf} pv_mono pv'

rule nsl_pv_all:
	input: vcf = 'pv/variants/our_goods_UG.pass.vcf'
	output: 'selscan/nsl/pv_all.svg'
	shell: 'bash selscan/nsl/nslRunner.sh {input.vcf} pv_all pv'

#######################
##### ADEGENET PCA ####
#######################

rule adegenet_pca:
	input: pv = 'adegenet/our_goods_pv.pass.sansAAKM.vcf', pf = 'adegenet/our_goods_pf.pass.vcf'
	output: grey = 'adegenet/grey_pca.svg', color = 'adegenet/color_pca.svg'
	shell: 'Rscript adegenet/adegenet.r {input.pv} {input.pf} {output.grey} {output.color}'

#######################
####### FWS ###########
#######################

rule calc_fws_boot:
	input: pv = 'pv/variants/our_goods_UG.pass.vcf', pf = 'pf/variants/our_goods_UG.pass.vcf'
	output: 'Fws/plot.svg'
	shell: 'Rscript Fws/fws.r {input.pv} {input.pf} {output}'

