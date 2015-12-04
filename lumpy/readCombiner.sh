## Need to combine reads so that I can skip the dedup step
## Because the headers get all messed up in dedup files
## So combine all reads that belong together, then realign
## Using bwa mem

readdir=/proj/julianog/users/ChristianP/cambodiaWGS/pv/symlinks
for name in `cat /proj/julianog/users/ChristianP/cambodiaWGS/pv/names/our_goods_5x@80%.txt`

do

echo $name

array=(`ls $readdir | grep $name`)
	## One lane
	if test ${#array[*]} = 2
	then
		echo "one run"
		ln -s $readdir/$name*_R1.fastq.gz catreads/$name\_R1.fq.gz
		ln -s $readdir/$name*_R2.fastq.gz catreads/$name\_R2.fq.gz
	else
		zcat $readdir/$name*_R1.fastq.gz > catreads/$name\_R1.fq.gz
		zcat $readdir/$name*_R2.fastq.gz > catreads/$name\_R2.fq.gz
	fi

done
