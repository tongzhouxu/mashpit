# Mashpit
Create a database of mash sketches and find the most similar genome to the target

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
  metadata_sra_db.py -source 0 -email your_email -key your_key -list project_list.txt
  ```
  - Using BioSample list
  ```
  metadata_sra_db.py -source 1 -list biosample_list.txt -email your_email -key your_key
  ```
  - Using keyword
  ```
  metadata_sra_db.py -source 2 -email your_email -key your_key -term salmonella_reading
  ```

#### Get the assembly and sketch files for all the entries in the database
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
  python3 query_against_db.py -n test_sample -database reading
  ```

## Installation

#### Dependencies

- Python 3.7
- Python packages:
  - Biopython
  - sqlite3
  - pandas
- Sourmash 3.3.1
