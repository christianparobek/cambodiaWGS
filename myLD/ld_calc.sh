## File to calculate LD using VCFtools
## It calculates point estimates for the CP groups
## As well as bootstraps and then calculates LD
## Started 5 November 2015
## Christian Parobek


########################
###### PARAMETERS ######
########################

#vcf2use=our_goods_pf.pass.vcf
vcf2use=our_goods_UG.pass80%.vcf

########################
## CP POINT ESTIMATES ##
########################

for i in 1 2 3 4
do

## split vcf
vcftools --vcf $vcf2use \
	--keep cp_groups/cp$i.txt \
	--recode --out pointestimates/vcfs/our.pf.cp$i

## calculate LD
vcftools --vcf pointestimates/vcfs/our.pf.cp$i.recode.vcf \
	--hap-r2 \
	--ld-window-bp 200000 \
	--out pointestimates/ld/cp$i.ld.1-200000

done


########################
####### BOOTSTRAP ######
########################

for strap in {1..100} # number of bootstraps to do
do

	## sample without replacement
	## sampling with replacement is impossible with vcftools
	## pull 18 isolates
	shuf -n18 cp_groups/all.txt >> bootstrap/boot/boot$strap.txt

	## subsample the VCF
	vcftools --vcf $vcf2use \
		--keep bootstrap/boot/boot$strap.txt \
		--recode \
		--out bootstrap/vcfs/boot$strap

	## calculate LD stats
	vcftools --vcf bootstrap/vcfs/boot$strap.recode.vcf \
		--hap-r2 \
		--ld-window-bp 200000 \
		--out bootstrap/ld/boot$strap.ld.1-200000

echo $strap
done

