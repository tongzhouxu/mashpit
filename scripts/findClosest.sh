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
# TODO might want to get taxid to help filter for later query
BASENAME=$(sqlite3 $DB "
  SELECT SKETCH.path 
  FROM SKETCH
  WHERE biosample_acc='$BIOSAMPLE_ACC'
  LIMIT 1
  ";
)
REFFILE="$DB.sketches/$BASENAME.msh"

# Get a list of files that the database captures
# TODO might want to filter by taxid
FOFN="$TMPDIR/fofn.txt"
sqlite3 $DB "SELECT path FROM SKETCH" |\
  while read BASENAME; do
    echo "$DB.sketches/$BASENAME.msh"
  done > $FOFN

DISTFILE="$TMPDIR/dist.tsv"
mash dist -t $REFFILE -l $FOFN | grep -v '^#' | sort -k2,2n > $DISTFILE

FILE_MATCHES=""
SRR=""
while read -r path dist; do
  if (( $(echo "$dist > $MAXDIST" | bc -l) )); then
    break;
  fi
  FILE_MATCHES="$FILE_MATCHES $path"
  IFS=. read -r biosample srr ext <<< "$path"
  SRR="$SRR $srr"
done < $DISTFILE

function join_by { local ifs="$1"; shift; str=""; for i in $@; do str="$str$i$ifs"; done; str=${str/%$ifs/}; echo $str;}

IN=$(join_by "','" $SRR)
sqlite3 -nullvalue . -separator $'\t' -header $DB "
  /*SELECT S.sketchid, S.biosample_acc, S.srr, S.path, B.strain, B.isolate, B.serovar, B.host_taxid*/
  SELECT B.*, S.*
  FROM SKETCH as S
  LEFT JOIN BIOSAMPLE as B ON B.biosample_acc=S.biosample_acc
  WHERE SRR IN ('$IN')
"

