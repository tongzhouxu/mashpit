#!/bin/bash

set -e
set -u

VERSION=2

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
BIOSAMPLE_ACC=$2

# Ensure that the sketches directory exists
mkdir -pv $DB.sketches

LOCKFILE=$DB.sketches/.lock
# Remove the lockfile if any crash.
trap ' { rm -rf $LOCKFILE; } ' SIGHUP SIGINT SIGTERM
lockfile $LOCKFILE

# Create a temporary directory
TMPDIR=$(mktemp -d MASHPIT.XXXXXX)
trap ' { rm -rf $TMPDIR; } ' EXIT

# First query: Get information about biosample
biosampleXML=$(esearch -db biosample -query "$BIOSAMPLE_ACC[accn]")

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

# Which taxid are we using?
if [ $SUBSPECIESTAXID -gt 0 ]; then
  TAXID=$SUBSPECIESTAXID
elif [ $SPECIESTAXID -gt 0 ]; then
  TAXID=$SPECIESTAXID
elif [ $GENUSTAXID -gt 0 ]; then
  TAXID=$GENUSTAXID
fi

# Add to taxonomy database if the taxon doesn't
# already exist. 
# First check if the taxon exists.
TAXID_DB=$(sqlite3 $DB "
  SELECT taxid 
  FROM TAXONOMY 
  WHERE taxid='$TAXID'
    AND genus='$GENUS'
    AND species='$SPECIES'
    AND subspecies='$SUBSPECIES'
  "
);
# If the taxon doesn't exist, make an entry. The key for
# this entry is $TAXONID
if [ ! "$TAXID_DB" ]; then
  sqlite3 $DB "
    INSERT INTO TAXONOMY (taxid, genus, species, subspecies)
    VALUES ($TAXID, '$GENUS', '$SPECIES', '$SUBSPECIES');
  ";
fi

# Figure out the biosample
BIOSAMPLE_ACC_DB=$(sqlite3 $DB "
  SELECT biosample_acc
  FROM BIOSAMPLE
  WHERE biosample_acc='$BIOSAMPLE_ACC'
  ";
);
# If the biosample isn't in the database already, then
# make an entry in the database.
if [ ! "$BIOSAMPLE_ACC_DB" ]; then
  biosampleXML_FULL=$(echo "$biosampleXML" | efetch -format xml);
  #echo "$biosampleXML_FULL" | xtract -pattern BioSample -group Attributes -elg "...\n..." -element Attribute@attribute_name Attribute -element Attribute@harmonized_name -element Attribute@display_name -element Attribute| sed 's/^\t//'
  #STRAIN=$(echo "$biosampleXML_FULL" | xtract -pattern BioSample -group Attributes -block Attribute -if Attribute@attribute_name -equals "sample name" -element Attribute)
  STRAIN=$(echo "$biosampleXML_FULL" | xtract -pattern BioSample -group Attributes -block Attribute -if Attribute@attribute_name -equals "strain" -element Attribute)
  ISOLATE=0;
  COLLECTED_BY=$(echo "$biosampleXML_FULL" | xtract -pattern BioSample -group Attributes -block Attribute -if Attribute@attribute_name -equals "collected-by" -element Attribute)
  LATITUDE=$(echo "$biosampleXML_FULL" | xtract -pattern BioSample -group Attributes -block Attribute -if Attribute@attribute_name -equals "lat_lon" -element Attribute)
  LONGITUDE=$(echo "$biosampleXML_FULL" | xtract -pattern BioSample -group Attributes -block Attribute -if Attribute@attribute_name -equals "lat_lon" -element Attribute)
  HOST=$(echo "$biosampleXML_FULL" | xtract -pattern BioSample -group Attributes -block Attribute -if Attribute@attribute_name -equals "specific host" -element Attribute)
  HOST_TAXID=0; # Need to translate $HOST back to taxid
  HOST_DISEASE=$(echo "$biosampleXML_FULL" | xtract -pattern BioSample -group Attributes -block Attribute -if Attribute@attribute_name -equals "host disease" -element Attribute)
  ISOLATION_SOURCE=0; # TODO
  SEROVAR=$(echo "$biosampleXML_FULL" | xtract -pattern BioSample -group Attributes -block Attribute -if Attribute@attribute_name -equals "serogroup" -element Attribute)

  sqlite3 $DB "
    INSERT INTO BIOSAMPLE (biosample_acc, strain, isolate, taxid, collected_by, collection_date, latitude, longitude, host_taxid, host_disease, isolation_source, serovar)
    VALUES ('$BIOSAMPLE_ACC', '$STRAIN', '$ISOLATE', '$TAXID', '$COLLECTED_BY', '$LATITUDE', '$LONGITUDE', '$HOST', '$HOST_TAXID', '$HOST_DISEASE', '$ISOLATION_SOURCE', '$SEROVAR')
    ;
  "
fi

# Sketch the sample

# To sketch the sample, I need the SRR identifier
# TODO capture the SRR somehow in the database
srrXML=$(echo "$biosampleXML" | elink -target SRA | efetch -format xml);
SRR=$(echo "$srrXML" | xtract -pattern EXPERIMENT_PACKAGE -group RUN_SET -element RUN@accession | tail -n 1)

sqlite3 $DB "
  INSERT INTO SRA (srr, biosample_acc)
  VALUES ('$SRR', '$BIOSAMPLE_ACC');
"

# grab the assembly if it exists
# https://github.com/ncbi/SKESA/issues/12#issuecomment-431503915
rp=$(srapath -f names -r $SRR.realign | awk '-F|' 'NF>8 && $(NF-1)==200 { print $8;}'  ) ; 
dump-ref-fasta "$rp" > $TMPDIR/$BIOSAMPLE_ACC.$SRR.fasta 

SKETCH_PATH="$DB.sketches/$BIOSAMPLE_ACC.$SRR.fasta.msh"
SKETCH_ID=""
if [ -s "$SKETCH_PATH" ]; then
  SKETCH_ID=$(sqlite3 $DB "
    SELECT sketchid
    FROM SKETCH
    WHERE path='$SKETCH_PATH';"
  );
else
  # Make a subshell to temporarily change directory
  (
    cd $TMPDIR
    mash sketch $BIOSAMPLE_ACC.$SRR.fasta > sketch.log 2>&1 || \
      cat sketch.log
  )
  mv $TMPDIR/$BIOSAMPLE_ACC.$SRR.fasta.msh $SKETCH_PATH
  SKETCH_BASENAME=$(basename $SKETCH_PATH .msh)
  # INSERT INTO DB
  SKETCH_ID=$(sqlite3 $DB "
    INSERT INTO SKETCH (biosample_acc, srr, path, source, software, seed)
    VALUES ('$BIOSAMPLE_ACC', '$SRR', '$SKETCH_BASENAME', 'NCBI download', 'Mash', '42')
    ;
    SELECT last_insert_rowid();
    "
  )
fi

#echo "SKETCH_ID => $SKETCH_ID"

rm -f $LOCKFILE

