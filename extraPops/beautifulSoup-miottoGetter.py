## BeautifulSoup work
## This package for python can parse HTML
## Need it to navigate the NCBI SRA website
## Want to download Miotto's 825 P. falciparum genomes
## Starting from the list of accessions provided by Derrick
## For Python 2.6+ and Python 3


## Import Libraries
from bs4 import BeautifulSoup
import urllib2
import re



## Get a list of all the ERS numbers, from csv Derrick provided
## Had to change \t to , and removed spaces that appear before commas

with open('thousandPf_list.csv') as f:
	lines = f.readlines()

ERSs = []
Cam_ERSs = []

## Miotto has multiple ERS numbers for some samples
## Maybe we just want to take the first ERS number for each sample?
## That would keep it simple
## Also, should I divide up by Country? RE match on Cambodia, Burkina, etc?

for line in lines:
	if re.search('ERS', line):
		line = line.split(',')
		ERSs.append(line[1])

# If I want to get all ERS numbers for any given population
# Replace the last line above with the following:
#		for i in line:
#			if re.search('ERS', i):
#				print i
#				Cam_ERSs.append(i)

## Now that I have the ERS numbers, get corresponding ERR numbers

ERRs = []
Cam_ERRs = []

for ers in ERSs:
	url = 'http://www.ncbi.nlm.nih.gov/sra/?term='+ers
	page = urllib2.urlopen(url)
	soup = BeautifulSoup(page.read())
	for result in soup.find_all('td'):
		if re.search('ERR', result.text):
			print result.text
			ERRs.append(result.text)
