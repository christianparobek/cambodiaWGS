## BeautifulSoup work
## This package for python can parse HTML
## Need it to navigate the NCBI SRA website
## Want to download Jane Carlton's ~175 P. vivax genomes
## Starting from the Pv HS BioProject page
## For Python 2.6+ and Python 3


## Import Libraries
from bs4 import BeautifulSoup
import urllib2
import re


## Get BioProject Accession Numbers (PRJNA)
url='http://www.ncbi.nlm.nih.gov/bioproject/240534'
page = urllib2.urlopen(url)
soup = BeautifulSoup(page.read())

prjnas = []
for prjna in soup.find_all('a'):
	if re.search('PRJNA', prjna.text):
		print(prjna.text)
		prjnas.append(prjna.text)


## Getting the SRX (isolate) numbers
srxs = []
for prjna in prjnas:
	url= 'http://www.ncbi.nlm.nih.gov/sra/?term='+prjna
	page = urllib2.urlopen(url)
	soup = BeautifulSoup(page.read())
	for result in soup.find_all('div', class_='rslt'):
		if re.search('HiSeq', result.text):
			for srx in result.find_all('dd'):
				print(srx.text)
				srxs.append(srx.text)


## Getting the SRR numbers
## It seems like there can be multiple SRR numbers for each isolate
srrs = []
for srx in srxs:
	url = 'http://www.ncbi.nlm.nih.gov/sra/?term='+srx
	page = urllib2.urlopen(url)
	soup = BeautifulSoup(page.read())
	for srr in soup.find_all('a'):
		if re.search('SRR', srr.text):
			print(srr.text)
			srrs.append(srr.text)


## Getting metadata
strains = []
dates = []
countries = []
for srr in srrs:
	print srr
	url = 'http://www.ncbi.nlm.nih.gov/sra/?term='+srr
	page = urllib2.urlopen(url)
	soup = BeautifulSoup(page.read())
	for meta in soup.find_all('span'):
		if re.search('vivax strain', meta.text):
			print(meta.text[24:])
			strains.append(meta.text[24:])
		if re.search('date', meta.text):
			print(meta.text[17:])
			dates.append(meta.text[17:])
		if re.search('country', meta.text):
			print(meta.text[50:])
			countries.append(meta.text[50:])






