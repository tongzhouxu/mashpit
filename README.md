[![Build Status](https://www.travis-ci.com/tongzhouxu/mashpit.svg?branch=master)](https://www.travis-ci.com/tongzhouxu/mashpit)
[![PyPI release](https://img.shields.io/pypi/v/mashpit)](https://pypi.python.org/pypi/mashpit/)
# Mashpit
Create a database of mash signatures and find the most similar genomes to a target sample 

## Installation
- Sra-tools for Linux
  ```
  wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/sratoolkit.2.10.8-centos_linux64.tar.gz -O /tmp/sratoolkit.tar.gz
  ```
  ```
  tar -xvf /tmp/sratoolkit.tar.gz
  ```
  Add sratoolkit to the environment:
  ```
  export PATH=$PATH:$PWD/sratoolkit.2.10.8-centos_linux64/bin
  ```

- Mashpit can be downloaded using pip or conda:
  ```
  pip install mashpit
  ```

- Ngstool is needed to run mashpit on Raspberry Pi
  1. Build and install ngs. Follow instructions at [https://github.com/ncbi/ngs/wiki/Building-and-Installing-from-Source](https://github.com/ncbi/ngs/wiki/Building-and-Installing-from-Source)
  2. Install ncbi-vdb. Follow instructions at [https://github.com/ncbi/ncbi-vdb/wiki/Building-and-Installing-from-Source](https://github.com/ncbi/ncbi-vdb/wiki/Building-and-Installing-from-Source)
  3. Build and install ngstools from source. Follow instructions at [https://github.com/ncbi/ngs-tools](https://github.com/ncbi/ngs-tools)

## Dependencies

- Python >= 3.6 and <=3.8
- Sra-tools 2.10.8

## Usage

#### 1. Create the database
```
usage: mashpit create [-h] database

Create new mashpit database

positional arguments:
  database    Name for the database.

optional arguments:
  -h, --help  show this help message and exit
```
- Example command
```
mashpit create reading
```
#### 2. Set up Entrez email and API key
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
mashpit config email@email.com -k p@$$word
```
More information about Entrez API key can be found on [this page.](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/)
#### 3. Collect the metadata
```
usage: mashpit metadata [-h] [-l LIST] [-t TERM]
                        database {bioproject_list,biosample_list,keyword}

Collect metadata from NCBI based on bioproject/biosample accession or keywords

positional arguments:
  database              Name of the database
  {bioproject_list,biosample_list,keyword}
                        Metadata collecting method. Available options:
                        bioproject_list, biosample_list, keyword

optional arguments:
  -h, --help            show this help message and exit
  -l LIST, --list LIST  File name of a list of bioproject or biosample
  -t TERM, --term TERM  Query keyword
```
- Example command
  - Using BioProject list
  ```
  mashpit metadata reading bioproject_list -l list_file
  ```
  - Using BioSample list
  ```
  mashpit metadata reading biosample_list -l list_file
  ```
  - Using keyword
  ```
  mashpit metadata reading keyword -t salmonella_reading
  ```

#### 4. Get the assembly and the signature file for all the entries in the database
```
usage: mashpit sketch [-h] [-n NUMBER] database

Build sketches for the records in the database

positional arguments:
  database              Name of the database

optional arguments:
  -h, --help            show this help message and exit
  -n NUMBER, --number NUMBER
                        Number of genomes in a batch to be downloaded and
                        sketched. Default is 1000.
```
- Example command
```
mashpit sketch reading
```
#### 5.Split the large signature file into separate ones to speed up the query
```
usage: mashpit split [-h] [-n NUMBER] database

Split large signature file to speed up the query

positional arguments:
  database              Name of the database

optional arguments:
  -h, --help            show this help message and exit
  -n NUMBER, --number NUMBER
                        Number of files to be splited into. Default is 16.
```
- Example command
```
mashpit split reading -n 16
```
#### 6. Query in the database
```
usage: mashpit query [-h] [-n NUMBER] [-f] sample database

Find the most similar assemblies to the target sample

positional arguments:
  sample                target sample file path
  database              name of the database

optional arguments:
  -h, --help            show this help message and exit
  -n NUMBER, --number NUMBER
                        number of separated signature file
  -f, --force           overwrite query record if query table exists
```
- Example command
  ```
  mashpit query sample_file reading -n 16
  ```
