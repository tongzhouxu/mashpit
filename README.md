# Mashpit

Create a database of mash sketches

## Usage

### Create the database

    sqlite3 mashpit.sqlite < scripts/mashpit_init.sql

### Add to the database

    # Add one
    scripts/addFromNcbi.sh mashpit.sqlite SAMN02182865
    # Add the next ten
    tail -n +2 t/data/biosamples.txt | tail | xargs -P 1 -n 1 bash -c '
      scripts/addFromNcbi.sh mashpit.sqlite $0
    '

### View sketches

    ls mashpit.sqlite.sketches/

    sqlite3 mashpit.sqlite '
    SELECT DISTINCT * 
    FROM BIOSAMPLE 
    INNER JOIN SKETCH 
      ON BIOSAMPLE.biosample_acc=SKETCH.biosample_acc
    '

## Installation

### Dependencies

* sra toolkit >= v12.5
* edirect >= v2.9
* Mash >= v2.0
* sqlite v3

## TODO

* method for adding local assembly
* Wrapper script for comparing two genomes

