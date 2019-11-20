#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw/GetOptions/;
use Data::Dumper;
use Bio::DB::EUtilities;
use Sys::Hostname qw/hostname/;

use FindBin qw/$RealBin/;
use lib "$RealBin/../lib";
use Mashpit;

my $email = $ENV{USER} . '@' . hostname();
exit(main());

sub main{
  my $settings={};
  GetOptions($settings, qw(help)) or die $!;
  usage() if($$settings{help} || @ARGV < 2);
  
  my($db, $biosample_acc) = @ARGV;

  die "ERROR: database does not exist: $db" if(!-e $db);

  my $mp = Mashpit->new($db);
  my $biosampleInfo = getBiosampleInfo($biosample_acc, $settings);

  return 0;
}

sub getBiosampleInfo{
  my($biosample_acc, $settings) = @_;
  my $esearch = Bio::DB::EUtilities->new(-eutil  => "esearch",
                                         -db     => "biosample",
                                         -term   => "$biosample_acc",
                                         -email  => $email,
                                         -usehistory => "y",
                                        );
  my @ids = $esearch->get_ids;

  while(my $docsum = $esearch->next_DocSum){
    die Dumper $docsum, 1;
  }
  die;
  for my $id(@ids){
    my $efetch = Bio::DB::EUtilities->new(-eutil   => "efetch",
                                          -db      => "biosample",
                                          -id      => $id,
                                         );
    die Dumper $efetch, $efetch->get_ids;
  }

  die "ERROR: no biosample returned using $biosample_acc";
}

sub usage{
  print "Usage: $0 mashpitDb biosample_acc
  \n";
  exit 0;
}

__END__

# Ensure that the sketches directory exists
mkdir -pv $DB.sketches

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


