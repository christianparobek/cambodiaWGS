## Started 25 September 2015
## To go from single, unannotated Pf and Pv VCFs
## And make annotated VCFs, split those VCFs by
## site, then make the split VCFs into dadi format.
## Christian Parobek

############################################################
################# DECLARE USEFUL PATHS #####################
############################################################

# directories
wd=/proj/julianog/users/ChristianP/cambodiaWGS/dadi/changFig1Maker
camWGS=/proj/julianog/users/ChristianP/cambodiaWGS

# data files
pfRef=/proj/julianog/refs/Pf3D7_v13.0/PlasmoDB-13.0_Pfalciparum3D7_Genome.fasta
pvRef=/proj/julianog/refs/PvSAL1_v13.0/PlasmoDB-13.0_PvivaxSal1_Genome.fasta

# scripts / programs
gatk=/nas02/apps/biojars-1.0/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar
snpEff=/nas02/home/p/r/prchrist/snpEff/snpEff.jar
SnpSift=/nas02/home/p/r/prchrist/snpEff/SnpSift.jar
vcf2dadi=/proj/julianog/users/ChristianP/cambodiaWGS/dadi/vcf2dadi.py

############################################################
################ ANNOTATE PV AND PF VCFS ###################
############################################################

## Run snpEff for Pf
java -Xmx4g -jar $snpEff -v -o gatk Pf3D7v13 $camWGS/pf/variants/our_goods_UG.pass.vcf > $wd/pf/our.pf.anno.vcf

## Run snpEff for Pv
java -Xmx4g -jar $snpEff -v -o gatk PvSAL1v13 $camWGS/pv/variants/our_goods_UG.pass.vcf > $wd/pv/our.pv.anno.vcf


############################################################
################ SPLIT VCFS BY SNP TYPE ####################
############################################################

for p in "pf" "pv"
do

	# Get the nonsyn snps
	java -jar $SnpSift filter "( EFF[*].EFFECT = 'NON_SYNONYMOUS_CODING' )" $wd/$p/our.$p.anno.vcf > $wd/$p/our.$p.nonsyn.vcf
	  
	# Get the syn snps
	java -jar $SnpSift filter "( EFF[*].EFFECT = 'SYNONYMOUS_CODING' )" $wd/$p/our.$p.anno.vcf > $wd/$p/our.$p.syn.vcf
	  
	# Get the genic snps
	java -jar $SnpSift filter "( EFF[*].EFFECT = 'SYNONYMOUS_CODING' ) | ( EFF[*].EFFECT = 'NON_SYNONYMOUS_CODING' )" $wd/$p/our.$p.anno.vcf > $wd/$p/our.$p.genic.vcf
	  
	# Get the intergenic snps
	java -jar $SnpSift filter "( EFF[*].EFFECT != 'SYNONYMOUS_CODING' ) & ( EFF[*].EFFECT != 'NON_SYNONYMOUS_CODING' )" $wd/$p/our.$p.anno.vcf > $wd/$p/our.$p.intergenic.vcf

done


############################################################
########### MAKE DADI FORMAT USING PY SCRIPT ###############
############################################################

module remove python
module add python/2.6.5

# For P falciparum
python $vcf2dadi $pfRef pf/our.pf.nonsyn.vcf pf/our.pf.nonsyn.dadi
python $vcf2dadi $pfRef pf/our.pf.syn.vcf pf/our.pf.syn.dadi
python $vcf2dadi $pfRef pf/our.pf.genic.vcf pf/our.pf.genic.dadi
python $vcf2dadi $pfRef pf/our.pf.intergenic.vcf pf/our.pf.intergenic.dadi

# for P vivax
python $vcf2dadi $pvRef pv/our.pv.nonsyn.vcf pv/our.pv.nonsyn.dadi
python $vcf2dadi $pvRef pv/our.pv.syn.vcf pv/our.pv.syn.dadi
python $vcf2dadi $pvRef pv/our.pv.genic.vcf pv/our.pv.genic.dadi
python $vcf2dadi $pvRef pv/our.pv.intergenic.vcf pv/our.pv.intergenic.dadi


