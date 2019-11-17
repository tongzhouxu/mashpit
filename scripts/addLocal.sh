#!/bin/bash

set -e
set -u

VERSION=1

function usage(){
  echo "Usage: $(basename $0) [options] mashpit.sqlite3 infile"
  echo "  infile can be a genome assembly in fasta format"
  echo "  or a raw reads set in fastq.gz format."
  echo "  -t taxid"
  echo "  -b biosampleId"
  echo "  -s SRR"
  exit 0
}


while getopts "ht:b:s:" opt; do
  case ${opt} in
    s )
      SRR=$OPTARG
      ;;
    b )
      BIOSAMPLE_ACC=$OPTARG
      ;;
    t )
      TAXID=$OPTARG
      ;;
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
shift $((OPTIND -1))
      
# Take care of two positional variables
if [ $# -lt 2 ]; then
  usage
fi
DB=$1
INPATH=$2

# Ensure that the sketches directory exists
mkdir -pv $DB.sketches

# Create a temporary directory
TMPDIR=$(mktemp -d MASHPIT.XXXXXX)
trap ' { rm -rf $TMPDIR; } ' EXIT

# Sketch the assembly
SKETCH_BASENAME=$(basename $INPATH)
SKETCH_PATH="$DB.sketches/$SKETCH_BASENAME.msh"
# Make a subshell to temporarily change directory
(
  ln -s $(realpath $INPATH) $TMPDIR/
  cd $TMPDIR
  ls -lh
  mash sketch $SKETCH_BASENAME > sketch.log 2>&1 || \
    (echo "ERROR:"; cat sketch.log; exit 1;)
)
mv $TMPDIR/$SKETCH_BASENAME.msh $SKETCH_PATH
# INSERT INTO DB
SKETCH_ID=$(sqlite3 $DB "
  INSERT INTO SKETCH (biosample_acc, srr, path, source, software, seed)
  VALUES ('$BIOSAMPLE_ACC', '$SRR', '$SKETCH_BASENAME', 'local', 'Mash', '42')
  ;
  SELECT last_insert_rowid();
  "
)

#echo "SKETCH_ID => $SKETCH_ID"


