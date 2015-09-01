## quick script to fix my broken symlinks

## first, identify broken pf links
find pf/aln/ -type l -xtype l -name "*merged.bam" | sed 's/pf\/aln\///' | sed 's/\.merged\.bam//' > pf_broken_syms


## then, identify broken pv links
find pv/aln/ -type l -xtype l -name "*merged.bam" | sed 's/pv\/aln\///' | sed 's/\.merged\.bam//' > pv_broken_syms


## fix the pf links
for name in `cat pf_broken_syms`
do

## remove the old link
rm /proj/julianog/users/ChristianP/cambodiaWGS/pf/aln/$name.merged.bam
## make the new link
ln -s /proj/julianog/users/ChristianP/cambodiaWGS/pf/aln/$name*.dedup.bam /proj/julianog/users/ChristianP/cambodiaWGS/pf/aln/$name.merged.bam

done


## fix the pv links
for name in `cat pv_broken_syms`
do

## remove the old link
rm /proj/julianog/users/ChristianP/cambodiaWGS/pv/aln/$name.merged.bam
## make the new link
ln -s /proj/julianog/users/ChristianP/cambodiaWGS/pv/aln/$name*.dedup.bam /proj/julianog/users/ChristianP/cambodiaWGS/pv/aln/$name.merged.bam

done


## remove my pf and pv name files
rm pf_broken_syms
rm pv_broken_syms
