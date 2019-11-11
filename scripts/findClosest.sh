#!/bin/bash

set -e
set -u

VERSION=1

function usage(){
  echo "Usage: $(basename $0) [options] mashpit.sqlite3 biosample_acc"
  echo "  -h help"
  echo "  -m maximum mash distance [Default: 0.1]"
  exit 0
}


MAXDIST=0.1
while getopts "hm:" opt; do
  case ${opt} in
    h )
      usage;
      exit 0;
      ;;
    m )
      MAXDIST=$OPTARG
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
BIOSAMPLE_ACC=$2

# Ensure that the sketches directory exists
mkdir -pv $DB.sketches

# Create a temporary directory
TMPDIR=$(mktemp -d MASHPIT.XXXXXX)
trap ' { rm -rf $TMPDIR; } ' EXIT

# Get ref mash sketch location
BASENAME=$(sqlite3 $DB "
  SELECT SKETCH.path 
  FROM SKETCH
  WHERE biosample_acc='$BIOSAMPLE_ACC'
  LIMIT 1
  ";
)
REFFILE="$DB.sketches/$BASENAME.msh"
ls $REFFILE


# Get a list of files that the database captures
FOFN="$TMPDIR/fofn.txt"
sqlite3 $DB "SELECT path FROM SKETCH" |\
  while read BASENAME; do
    echo "$DB.sketches/$BASENAME.msh"
  done > $FOFN

DISTFILE="$TMPDIR/dist.tsv"
mash dist -t $REFFILE -l $FOFN | grep -v '^#' | sort -k2,2n > $DISTFILE

BIOSAMPLE_MATCHES=""
cat $DISTFILE | while read -r path dist; do
  if (( $(echo "$dist > $MAXDIST" | bc -l) )); then
    break;
  fi
  MATCH=$(sqlite3 $DB "
    SELECT biosample_acc
    FROM SKETCH
    WHERE path='$path'
    LIMIT 1"
  )
  echo $path $MATCH
  BIOSAMPLE_MATCHES="$BIOSAMPLE_MATCHES $MATCH"
done

echo "$BIOSAMPLE_MATCHES"

