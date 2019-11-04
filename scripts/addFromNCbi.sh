#!/bin/bash

set -e
set -u

VERSION=1

function usage(){
  echo "Usage: $(basename $0) mashpit.sqlite3 biosample_acc"
  exit 0
}


while getopts "h" opt; do
  case ${opt} in
    h )
      usage;
      exit 0;
      ;;
    : )
      echo "ERROR with $OPTARG"
      exit 1
      ;;
    \? )
      echo "Invalid option specified!"
      exit 1
      ;;
  esac
done
      
# Take care of two positional variables
if [ $# -lt 2 ]; then
  usage
fi
DB=$1
BIOSAMPLE=$2

# First query: Get information about biosample
biosampleXML=$(esearch -db biosample -query "$BIOSAMPLE[accn]")

# Second query: Get taxonomy information
taxonXML=$(echo "$biosampleXML" | elink -target taxonomy | efetch -format xml)
taxonTSV=$(echo "$taxonXML" | xtract -pattern LineageEx -group Taxon -elg "\n" -element TaxId -element ScientificName -element Rank | sed s'/^\t//')

GENUS=""
SPECIES=""
SUBSPECIES=""
TAXID=0
GENUSTAXID=0
SPECIESTAXID=0
SUBSPECIESTAXID=0
while IFS=$'\t' read -r taxid name rank; do
  if [ "$rank" == "genus" ]; then
    GENUS=$name
    GENUSTAXID=$taxid
  elif [ "$rank" == "species" ]; then
    SPECIES=$name
    SPECIESTAXID=$taxid
  elif [ "$rank" == "subspecies" ]; then
    SUBSPECIES=$name
    SUBSPECIESTAXID=$taxid
  fi
done <<< "$taxonTSV"

# TODO logic on which taxid to take

