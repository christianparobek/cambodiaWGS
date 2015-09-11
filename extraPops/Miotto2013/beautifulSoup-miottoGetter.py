## BeautifulSoup work
## This package for python can parse HTML
## Need it to navigate the NCBI SRA website
## Want to download Miotto's 825 P. falciparum genomes
## Starting from the list of accessions provided by Derrick
## For Python 2.6+ and Python 3

##########################
#### Import Libraries ####
##########################
from bs4 import BeautifulSoup
import urllib2
import re

###########################
## Define some functions ##
###########################

## errScraper should take a webpage and 
## scrape all ERR numbers off it and
## return a list of ERR numbers
def errScraper(url,listOfErrs):
	page = urllib2.urlopen(url)
	soup = BeautifulSoup(page.read())
	runs = []
	for result in soup.find_all('td'):
		if re.search('ERR', result.text):
			runs.append(result.text[0:9])
	if len(runs) == 1:
		listOfErrs.append(runs[0])
	elif len(runs) > 1:
		listOfErrs.append(runs)
	return


## erxScraper should take a webpage with
## ERX accession numbers and return a list
## of ERX numbers
def erxScraper(url,listOfErxs):
	page = urllib2.urlopen(url)
	soup = BeautifulSoup(page.read())
	for result in soup.find_all('dd'):
		if re.search('ERX', result.text):
			listOfErxs.append(result.text)
	return


## filePrinter should take a list of
## SRA accession numbers (eg. ERRs)
## and print to a pre-defined outfile
def filePrinter(array,output):
	for line in array:
		print>>output, line



###########################
##### Input & Output ######
###########################

## Infile
## Get a list of all the ERS numbers, from csv Derrick provided
## Had to change \t to , and removed spaces that appear before commas
with open('thousandPf_list.csv') as f:
	lines = f.readlines()

## Outfile
ERRs_googleDocs = open('ERRs_google.txt', 'w')
ERRs_forDownload = open('ERRs_dnload.txt', 'w')

###########################
#### Getting the ERRs #####
###########################

## Miotto has multiple ERS numbers for some samples
## Miotto has multiple libraries for some ERS numbers
## Miotto has multiple runs for some libraries
## Get all of them

## First, get all ERS numbers for each sample

ERSs = []

for line in lines:
	line = line.split(',')
	sampERSs = []
	for col in line:
		if re.search('ERS', col):
			sampERSs.append(col)
	ERSs.append(sampERSs)

## Now that I have all ERS nums for all samples, get corresponding ERR numbers
## Put the ERR numbers into list-of-lists format, useful for google docs

ERRs = []

for samp in ERSs:
	sampERRs = []
	for ers in samp:
		url = 'http://www.ncbi.nlm.nih.gov/sra/?term='+ers
		page = urllib2.urlopen(url)
		soup = BeautifulSoup(page.read())
		print url
		for result in soup.find_all('div'):
			if re.search('ERR', result.text):
				errScraper(url, sampERRs)
				break
			elif re.search('ERX', result.text): 
				sampERXs = []
				erxScraper(url, sampERXs)
				for erx in sampERXs:
					url = 'http://www.ncbi.nlm.nih.gov/sra/?term='+erx
					errScraper(url, sampERRs)
				break
	print sampERRs
	ERRs.append(sampERRs)

## Right now ERRs, is a list of lists, useful for the google doc
## However, for downloading purposes, we want to just get a list
ERRs_dnload = [item for sublist in ERRs for item in sublist]

filePrinter(ERRs,ERRs_googleDocs)
filePrinter(ERRs_dnload,ERRs_forDownload)


