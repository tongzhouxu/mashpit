#!/bin/bash
# Library of functions for Mashpit shell scripting

# Run an sql statement with the right pragmas
# Arguments:
#   database.sqlite
#   SQL statement
function sql(){
  db=$1
  sql=$2
  sqlite3 $db "pragma foreign_keys='on';
    $sql
  "
}

# Safely add a new taxon
# Arguments:
#   database.sqlite
#   Genus
#   Species
#   Subspecies
# Prints:
#   taxid
# Returns:
#   0 if new taxid
#   1 if Genus/species/subspecies already exists
function addTaxon(){
  db=$1
  genus=$2
  species=$3
  subspecies=$4
  
  taxid=$(sql $db "
    SELECT taxid
    FROM TAXONOMY
    WHERE genus='$genus'
      AND species='$species'
      AND subspecies='$subspecies'
    LIMIT 1;"
  );

  if [ "$taxid" -gt 0 ]; then
    echo $taxid;
    return 0;
  fi
  
  taxid=$(sql $db "
    INSERT INTO TAXONOMY (genus species subspecies)
    VALUES ('$genus','$species','$subspecies');
    SELECT last_insert_rowid();
  ");

  echo $taxid;
  return 0;
}

# Safely add a new biosample 
# Arguments:
#   database.sqlite
#   ...
# Prints:
#   biosampleid 
# Returns:
#   0 if new biosampleid 
#   1 if biosample already existed
function addBiosample(){
}
