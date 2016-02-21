## Need to combine reads so that I can skip the dedup step
## Because the headers get all messed up in dedup files
## So combine all reads that belong together, then realign
## Using bwa mem

readdir=/proj/julianog/users/ChristianP/cambodiaWGS/pf/symlinks


for name in `cat /proj/julianog/users/ChristianP/cambodiaWGS/pf/names/our_goods_5x@60%.txt`
do

echo $name

array=(`ls $readdir | grep $name | grep R1`)
	## One lane
	if test ${#array[*]} = 1
	then
		echo "one run"
		ln -s $readdir/$name*_R1.fastq.gz catreads/$name\_R1.fq.gz
		ln -s $readdir/$name*_R2.fastq.gz catreads/$name\_R2.fq.gz
	elif test ${#array[*]} = 2
	then
		echo "multiple runs"
		zcat $readdir/$name*_R1.fastq.gz | gzip > catreads/$name\_R1.fq.gz
		zcat $readdir/$name*_R2.fastq.gz | gzip > catreads/$name\_R2.fq.gz
	else
		echo "something is wrong"
	fi

done
