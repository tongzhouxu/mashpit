# Mashpit
Create a database of mash sketches and find the most similar genome to the target

## Usage

#### Create the database
```
python3 create_db.py
```

#### Collect the metadata
```
usage: metadata_sra_db.py -source <source> -list <name> -email <Entrez email address> -key <Entrez API key>

optional arguments:
  -h, --help      show this help message and exit
  -source SOURCE  <int>: source of data. 0: bioproject list, 1: biosample list, 2: name
  -list LIST      <string>: full name of the list file
  -email EMAIL    <string>: entrez email address
  -key KEY        <string>: entrez api key
  -term TERM      <string>: specie/serovar name
```
More information about Entrez API key can be found on [this page.](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/)
- Example command
  - Using BioProject list
  ```
  python3 metadata_sra_db.py -source 0 -email your_email -key your_key -list project_list.txt
  ```
  - Using BioSample list
  ```
  python3 metadata_sra_db.py -source 1 -list biosample_list.txt -email your_email -key your_key
  ```
  - Using keyword
  ```
  python3 metadata_sra_db.py -source 2 -email your_email -key your_key -term salmonella_reading
  ```

#### Get the assembly and sketch files for all the entries in the database
```
python3 sketch_db.py
```

#### Query in the database
```
usage: query_against_db.py -n <sample name>

optional arguments:
  -h, --help  show this help message and exit
  -n N        <string>: sample SRR name
```
- Example command
  ```
  python3 query_against_db.py -n SRR11305724
  ```

## Installation

#### Dependencies

- Python 3.7
- Python packages:
  - Biopython
  - sqlite3
- Mash 2.2.2
