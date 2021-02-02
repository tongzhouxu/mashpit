[![Build Status](https://www.travis-ci.com/tongzhouxu/mashpit.svg?branch=master)](https://www.travis-ci.com/tongzhouxu/mashpit)
[![PyPI release](https://img.shields.io/pypi/v/mashpit)](https://pypi.python.org/pypi/mashpit/)
# Mashpit
Create a database of mash signatures and find the most similar genomes to a target sample 

## Usage

#### Create the database
Create an empty sqlite db file in the current working directory.
```
usage: create_db.py <database name>

positional arguments:
  database    <string>: name of the database

optional arguments:
  -h, --help  show this help message and exit
```
- Example command
```
create_db.py reading
```
#### Collect the metadata
```
usage: metadata_sra_db.py <database name> <method> <Entrez email address> [-key <Entrez API key>] [-list <accession list>] [-term <keyword>]

positional arguments:
  database              <string>: name of the database
  method                <string>: data collecting method. Available options:
                        bioproject_list, biosample_list, keyword
  email                 <string>: Entrez email address

optional arguments:
  -h, --help            show this help message and exit
  -k KEY, --key KEY     <string>: Entrez api key
  -l LIST, --list LIST  <string>: list file of bioproject or biosample
  -t TERM, --term TERM  <string>: query keyword
```
More information about Entrez API key can be found on [this page.](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/)
- Example command
  - Using BioProject list
  ```
  metadata_sra_db.py database_name bioproject_list your_ncbi_email --key your_key --list project_list.txt
  ```
  - Using BioSample list
  ```
  metadata_sra_db.py database_name biosample_list your_ncbi_email --key your_key --list biosample_list.txt
  ```
  - Using keyword
  ```
  metadata_sra_db.py database_name keyword your_ncbi_email --key your_key --term salmonella_reading
  ```

#### Get the assembly and the signature file for all the entries in the database
```
usage: sketch_db.py <database name>

positional arguments:
  database    <string>: name of the database

optional arguments:
  -h, --help  show this help message and exit
```

#### Query in the database
```
usage: query_against_db.py <sample name> <database name> [--force, -f]

positional arguments:
  sample_name  <string>: sample file name
  database     <string>: name of the database

optional arguments:
  -h, --help   show this help message and exit
  -f, --force  overwrite if query table exists
```
- Example command
  ```
  query_against_db.py sample_name reading
  ```

## Installation
- Sra-tools
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

- Mashpit can be downloaded using pip:
  ```
  pip install mashpit
  ```


## Dependencies

- Python >= 3.6 and <=3.8
- Sra-tools 2.10.8

## Step-by-step instructions on building a Salmonella Bareilly mashpit database

#### 1. Get mashpit and dependent softwares ready

- Sra-tools
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

- Mashpit can be downloaded using pip:
  ```
  pip install mashpit
  ```

#### 2. Start by creating an empty sqlite database

  Optionally, you can create a folder for the database
  ```
  mkdir mashpit_salmonella_bareilly
  cd mashpit_salmonella_bareilly
  ```
  And then run:
  ```
  create_db.py bareilly
  ```
  In this command, bareilly is the name for the database. This command should generate a file named **bareilly.db**

#### 3. Build up the metadata database by searching for all the information in NCBI using keyword salmonella bareilly

  ```
  metadata_sra_db.py bareilly keyword your_ncbi_account_email --key your_ncbi_api_key --term salmonella_bareilly
  ```
  
 #### 4. Get the SKESA assembly and sourmash signature file
 
  ```
  sketch_db.py bareilly
  ```

  This command should create a folder named "tmp" in the current path, and the assembly downloaded will be stored in this folder during the sketching process. A signatrue file named "bareilly.sig" should also be generated in which there are sourmash signatures for each assembly in the database folder.
 
  While running this command, you will see "Assembly downloaded for SRR_accession" if successful. If the skesa assembly does not exist or there is anything wrong, you will see "No assembly for SRR_accession" printing out.
  
 #### 5. Get the most similar genome in the database for a target sample.
 
   ```
   query_against_db.py test_sample_name bareilly
   ```
   
   In this command, test_sample_name should be the full name of the target sample. The target sample should be a genome assembly and exist in the current path.
   
   By running this command, you will see the top 50 results printing out sorted according to the similarity (jaccard_scores). It should also generate a csv file named test_sample_name_output.csv with full results record.
