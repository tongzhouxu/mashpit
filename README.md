[![Build Status](https://www.travis-ci.com/tongzhouxu/mashpit.svg?branch=master)](https://www.travis-ci.com/tongzhouxu/mashpit)
[![PyPI release](https://img.shields.io/pypi/v/mashpit)](https://pypi.python.org/pypi/mashpit/)
# Mashpit
Create a database of mash signatures and find the most similar genomes to a target sample 

## Installation
Install sratoolkit for Linux
```
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/sratoolkit.2.10.8-centos_linux64.tar.gz -O /tmp/sratoolkit.tar.gz
tar -xvf /tmp/sratoolkit.tar.gz
export PATH=$PATH:$PWD/sratoolkit.2.10.8-centos_linux64/bin
```

Install mashpit using pip:
  ```
  pip install mashpit
  ```
Or manually install after git clone:
  ```
  git clone https://github.com/tongzhouxu/mashpit.git
  python setup.py install 
  ```
Ngstool is needed to build mashpit database on Raspberry Pi:
1. Build and install ngs. Follow instructions at [https://github.com/ncbi/ngs/wiki/Building-and-Installing-from-Source](https://github.com/ncbi/ngs/wiki/Building-and-Installing-from-Source)
2. Install ncbi-vdb. Follow instructions at [https://github.com/ncbi/ncbi-vdb/wiki/Building-and-Installing-from-Source](https://github.com/ncbi/ncbi-vdb/wiki/Building-and-Installing-from-Source)
3. Build and install ngstools from source. Follow instructions at [https://github.com/ncbi/ngs-tools](https://github.com/ncbi/ngs-tools)

## Dependencies

- Python >= 3.8
- Sratoolkit 2.10.8

## Mashpit Database

A mashpit database is a directory containing:
- `$DB_NAME.db`
- `$DB_NAME.sig`

Two types of Mashpit database:

1. Standard Database (Pathogen Detection Database)

   A standard database is a collection of representative genomes from each cluster on [Pathogen Detection](https://www.ncbi.nlm.nih.gov/pathogens). By default mashpit will download the latest version of a specified species and find the centroid of each SNP cluter (SNP tree).
2. Custom Database  
   
   A custom database is a collection of genomes based on a proveded biosample accesion list or a keyword.

## Usage

### 0. Set up Entrez email and API key (For custom database only)
```
usage: mashpit config [-h] [-k KEY] email

Add Entrez email and key to environment variables

positional arguments:
  email              Entrez email address

optional arguments:
  -h, --help         show this help message and exit
  -k KEY, --key KEY  Entrez api key
```
- Example command
```
mashpit config example@example.com -k your_password
```
More information about Entrez API key can be found on [this page.](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/)
### 1. Build the database
```
usage: mashpit build [-h] [-s SPECIES] [-l LIST] [-t TERM]
                     {standard,biosample_list,keyword} name

Build mashpit database

positional arguments:
  {standard,biosample_list,keyword}
                        Database type: standard or custom
  name                  Mashpit database name

optional arguments:
  -h, --help            show this help message and exit
  -s SPECIES, --species SPECIES
                        Query keyword
  -l LIST, --list LIST  File containing a list of biosample accessions
  -t TERM, --term TERM  Query keyword
```
- Example command
```
mashpit build standard salmonella -s Salmonella
```
```
mashpit build biosample_list db_name -l file.list
```
```
mashpit build keyword db_name -t salmonella_enteritidis
```
### 2. Sketch the genomes
```
usage: mashpit sketch [-h] name

Build sketches for the records in the database

positional arguments:
  name        Mashpit database name

optional arguments:
  -h, --help  show this help message and exit
```
- Example command
```
mashpit sketch db_name
```
### 3. Run a query
```
usage: mashpit query [-h] [-n NUMBER] sample database

Find the most similar assemblies to the target sample

positional arguments:
  sample                file name of the query sample
  database              name of the database

optional arguments:
  -h, --help            show this help message and exit
```
- Example command
```
mashpit query sample.fasta db_name
```
