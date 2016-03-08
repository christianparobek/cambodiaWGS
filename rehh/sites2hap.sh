## A shell script to convert 
## Started 07 March 2016
## Given a LDhat sites file, convert to the rehh HAP format
## The LDhat sites files were originally made for selscan iHS
## USAGE: bash sites2hap.sh <in.ldhat.sites> <out.hap>

in=$1
out=$2

#in=/proj/julianog/users/ChristianP/cambodiaWGS/selscan/ihs/pv_mono/01.ldhat.sites

tail -n+2 $in |\
	grep -v ">" |\
	sed 's/1/2/g' |\
	sed 's/0/1/g' |\
	sed 's/.\{1\}/& /g' |\
	awk '{printf("%01d %s\n", NR, $0)}' > $out

## This command does the following
##	(1) Skips the first line (header)
##	(2) Removes all the fasta info lines
##	(3) Converts 1s to 2s (derived), to conform to 1 or 2 format
##	(4) Converts 0s to 1s (ancestral), to conform to 1 or 2 format
##	(5) Adds spaces between each digit
##	(6) Adds line numbers to the beginning of each line

