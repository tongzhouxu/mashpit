# Mashpit
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/mashpit/README.html)
![unittest](https://github.com/tongzhouxu/mashpit/actions/workflows/python-app.yml/badge.svg)
[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![PyPI release](https://img.shields.io/pypi/v/mashpit)](https://pypi.python.org/pypi/mashpit/)

Create a database of mash signatures and find the most similar genomes to a target sample. To contribute, please see the contributing guidelines [here](CONTRIBUTING.md).


## Installation
### Option 1. Install with Conda/Mamba (Recommended)
```
conda create -n mashpit -c conda-forge -c bioconda 'mashpit=0.9.6'
conda activate mashpit
```
### Option 2. Install with pip
#### 1. Dependency: Install NCBI datasets
```
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
chmod +x datasets
export PATH=$PATH:$PWD
```

#### 2. Install mashpit using pip:
  ```
  pip install mashpit
  ```
#### Or git clone from github:
  ```
  git clone https://github.com/tongzhouxu/mashpit.git
  cd mashpit
  pip install . 
  ```

## Mashpit Database

A mashpit database is a directory containing:
- `$DB_NAME.db`
- `$DB_NAME.sig`

Mashpit database can be built using:

1. A taxonomic name 
   A standard database is a collection of representative genomes from each cluster on [Pathogen Detection](https://www.ncbi.nlm.nih.gov/pathogens). By default mashpit will download the latest version of a specified species and find the centroid of each SNP cluter (SNP tree).
2. BioSample accessions  
   A custom database is a collection of genomes based on a proveded biosample accesion list.

## Usage

### 1. Build a mashpit database
```
usage: mashpit build [-h] [--quiet] [--number NUMBER] [--ksize KSIZE] [--species SPECIES] [--email EMAIL] [--key KEY] [--pd_version PD_VERSION] [--list LIST] {taxon,accession} name

positional arguments:
  {taxon,accession}     mashpit database type.
  name                  mashpit database name

optional arguments:
  -h, --help            show this help message and exit
  --quiet               disable logs
  --number NUMBER       maximum number of hashes for sourmash, default is 1000
  --ksize KSIZE         kmer size for sourmash, default is 31
  --species SPECIES     species name
  --email EMAIL         Entrez email
  --key KEY             Entrez api key
  --pd_version PD_VERSION
                        a specified Pathogen Detection version (PDG accession). Default is the latest.
  --list LIST           Path to a list of NCBI BioSample accessions
```
- Example command
```
mashpit build taxon salmonella --species Salmonella
```

Note: Supported species names can be found in this [list](https://ftp.ncbi.nlm.nih.gov/pathogen/Results/)

### 2. Query against a mashpit database
```
usage: mashpit query [-h] [--number NUMBER] [--threshold THRESHOLD] [--annotation ANNOTATION] sample database

positional arguments:
  sample                path to query sample
  database              path to the database folder

optional arguments:
  -h, --help            show this help message and exit
  --number NUMBER       number of isolates in the query output, default is 200
  --threshold THRESHOLD
                        minimum jaccard similarity for mashtree, default is 0.85
  --annotation ANNOTATION
                        mashtree tip annoatation, default is none
```
- Example command
```
mashpit query sample.fasta path/to/database
```
### Optional: Update the database
```
usage: mashpit update [-h] [--metadata METADATA] [--quiet] database name

positional arguments:
  database             path for the database folder
  name                 database name

optional arguments:
  -h, --help           show this help message and exit
  --metadata METADATA  metadata file in csv format
  --quiet              disable logs
  ```
- Example command
```
mashpit update path/to/database salmonella
```
### Webserver
A local host webserver can be started using:
```
mashpit webserver
```
