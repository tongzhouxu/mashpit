# Mashpit
Create a database of mash signatures and find the most similar genomes to a target sample 

## Usage

#### Create the database
```
usage: create_db.py -database <database name>

optional arguments:
  -h, --help          show this help message and exit
  -database DATABASE  <string>: name of the database
```
- Example command
```
create_db.py -database reading
```
#### Collect the metadata
```
usage: metadata_sra_db.py -database <database name> -method <method> -list <accession list> (-term <keyword>) -email <Entrez email address> -key <Entrez API key>

optional arguments:
  -h, --help          show this help message and exit
  -database DATABASE  <string>: name of the database
  -method METHOD      <int>: data collecting method. 0: bioproject list, 1:
                      biosample list, 2: name
  -list LIST          <string>: full name of the list file
  -email EMAIL        <string>: entrez email address
  -key KEY            <string>: entrez api key
  -term TERM          <string>: specie/serovar name when using keyword to
                      collect data
```
More information about Entrez API key can be found on [this page.](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/)
- Example command
  - Using BioProject list
  ```
  metadata_sra_db.py -database database_name -method 0 -email your_email -key your_key -list project_list.txt
  ```
  - Using BioSample list
  ```
  metadata_sra_db.py -database database_name -method 1 -list biosample_list.txt -email your_email -key your_key
  ```
  - Using keyword
  ```
  metadata_sra_db.py -database database_name -method 2 -email your_email -key your_key -term salmonella_reading
  ```

#### Get the assembly and the signature file for all the entries in the database
```
usage: sketch_db.py -database <database name>

optional arguments:
  -h, --help          show this help message and exit
  -database DATABASE  <string>: name of the database
```

#### Query in the database
```
usage: query_against_db.py -n <sample name> -database <database name>

optional arguments:
  -h, --help          show this help message and exit
  -n N                <string>: sample file name
  -database DATABASE  <string>: name of the database
  -f, --force         overwrite if query table exists
```
- Example command
  ```
  query_against_db.py -n test_sample -database reading
  ```

## Dependencies

- Python 3.7
- Python packages:
  - Biopython
  - Pandas
  - Sourmash
- Sra-tools 2.10.7

## Step-by-step instructions on building a Salmonella Bareilly mashpit database

#### 1. Get mashpit and dependent softwares ready

- Python 3.7
Python can be downloaded from the [python official website.](https://www.python.org/downloads/)

- Python packages can be installed using commands:
  ```
  pip install biopython
  ```
  [More information about biopython.](https://biopython.org/wiki/Download)
  ```
  pip install pandas
  ```
  [More information about pandas.](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html)
  ```
  pip install sourmash
  ```
  [More information about sourmash.](https://pypi.org/project/sourmash/)

- Sra-tools
Installation information can be found [here.](https://github.com/ncbi/sra-tools)
Please make sure the **dump-ref-fasta** command in sratools is added to PATH before running mashpit.
  
- Mashpit can be downloaded using command:
  ```
  git clone https://github.com/tongzhouxu/mashpit.git
  ```
  python files in mashpit can be added to PATH using command:
  ```
  export PATH=/path-to-mashpit:$PATH
  ```

#### 2. Start by creating an empty sqlite database

  Optionally, you can create a folder for the database
  ```
  mkdir mashpit_salmonella_bareilly
  cd mashpit_salmonella_bareilly
  ```
  And then run:
  ```
  create_db.py -database bareilly
  ```
  In this command, bareilly is the name for the database. This command should generate a file named **bareilly.db**

#### 3. Build up the metadata database by searching for all the information in NCBI using keyword salmonella bareilly

  ```
  metadata_sra_db.py -database bareilly -method 2 -email your_ncbi_account_email -key your_ncbi_api_key -term salmonella_bareilly
  ```
  
  In this command, "-method 2" indicates we search for the information using a keyword. 

  While running this command, you will see all the metadata fetched printing out on the screen. If there is no SRA record for a biosample, the entry will be skipped and you will see "No SRA record, skipping it."
  
 #### 4. Get the SKESA assembly and sourmash signature file
 
  ```
  sketch_db.py -database bareilly
  ```

  This command should create a folder named "database" in the current path, and the assembly downloaded will be stored in this folder. A signatrue file named "database.sig" should also be generated in which there are sourmash signatures for each assembly in the database folder.
 
  While running this command, you will see "Downloaded assembly for SRR_accession" if successful. If the skesa assembly does not exist or there is anything wrong, you will see "Can't download SKESA assembly for SRR_accession" printing out.
  
 #### 5. Get the most similar genome in the database for a target sample.
 
   ```
   query_against_db.py -n test_sample_name -database bareilly
   ```
   
   In this command, test_sample_name should be the full name of the target sample. The target sample should be a genome assembly and exist in the current path.
   
   By running this command, you will see the top 50 results printing out sorted according to the similarity (jaccard_scores). It should also generate a csv file named test_sample_name_output.csv with full results record.
   
  
 
