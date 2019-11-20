#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw/GetOptions/;
use Data::Dumper;
use Bio::DB::EUtilities;
use Sys::Hostname qw/hostname/;
use List::MoreUtils qw/uniq/;

use FindBin qw/$RealBin/;
use lib "$RealBin/../lib";
use Mashpit qw/logmsg/;
use XML::Hash;

my $email = $ENV{USER} . '@' . hostname();
my $MP; # Global mashpit instance, set in main()
exit(main());

sub main{
  my $settings={};
  GetOptions($settings, qw(help)) or die $!;
  usage() if($$settings{help} || @ARGV < 2);
  
  my($db, $biosample_acc) = @ARGV;

  die "ERROR: database does not exist: $db" if(!-e $db);
  $MP = Mashpit->new($db);

  my $taxid         = getTaxonInfo($biosample_acc, $settings);
  my $biosampleInfo = getBiosampleInfo($biosample_acc, $settings);

  return 0;
}

# Get the taxonomy id and if it's in the database, great.
# If not in the database, add it.
sub getTaxonInfo{
  my($biosample_acc, $settings) = @_;
  my $taxid = 0;

  my $esearch = Bio::DB::EUtilities->new(-eutil  => "esearch",
                                         -db     => "biosample",
                                         -term   => "$biosample_acc",
                                         -email  => $email,
                                         -usehistory => "y",
                                        );
  my @id    = $esearch->get_ids;

  my $elink = Bio::DB::EUtilities->new(-eutil   => "elink",
                                       -db      => "taxonomy",
                                       -target  => "taxonomy",
                                       -dbFrom  => "biosample",
                                       -usehistory=>"y",
                                       -email   => $email,
                                       -id      => \@id,
                                      );
  my @taxid = uniq $elink->get_ids;

  if(@taxid > 1){
    logmsg "WARNING: found more than one taxid for $biosample_acc. Will only use $taxid[0].";
    logmsg "All taxids found: ". join(" ", @taxid);
  }

  my $efetch = Bio::DB::EUtilities->new(-eutil   => "efetch",
                                        -db      => "taxonomy",
                                        -email   => $email,
                                        -id      => \@taxid,
                                        -rettype => "xml",
                                       );
  my $count = $esearch->get_count;
  my($retry, $retmax, $retstart) = (0, 500, 0);
  RETRIEVE_TAXONOMY:
  while($retstart < $count){
    $efetch->set_parameters(-retmax   => $retmax,
                            -retstart => $retstart);
    my $xml = "";
    eval{
      $efetch->get_Response(-cb=>sub{
          my($data) = @_;
          $xml .= $data;
      });
    };
    if($@ || !$xml){
      die "Server error: $@. Try again later" if($retry==5);
      logmsg "Server error, redo #$retry";
      $retry++ && next RETRIEVE_TAXONOMY;
    }

    my $xmlConverter = XML::Hash->new();
    my $taxonomyHash = $xmlConverter->fromXMLStringtoHash($xml);

    # Now let's get a hash of %rank=> {genus=>{text=>Salmonella, taxid=>590}, species=>{...}, ...}
    my %rank;
    for my $rankInfo(@{ $$taxonomyHash{TaxaSet}{Taxon}{LineageEx}{Taxon} }){
      my($taxid, $rankKey, $rankValue) = ($$rankInfo{TaxId}{text}, $$rankInfo{Rank}{text}, $$rankInfo{ScientificName}{text});
      next if($rankKey eq 'no rank' || !$rankKey);
      $rank{$rankKey} = {taxid=>$taxid, text=>$rankValue};
    }
    
    # For each genus, species, subspecies, add the taxonomy
    # to the database if it isn't already in there.
    my $mpGenusTaxid = $MP->getTaxonomy($rank{genus}{taxid});
    if(!$mpGenusTaxid){
      $mpGenusTaxid = $MP->addTaxonomy($rank{genus}{taxid}, $rank{genus}{text});
    }
    my $mpSpeciesTaxid = $MP->getTaxonomy($rank{species}{taxid});
    if(!$mpSpeciesTaxid){
      $mpSpeciesTaxid = $MP->addTaxonomy($rank{species}{taxid}, $rank{genus}{text}, $rank{species}{text});
    }
    my $mpSubspeciesTaxid = $MP->getTaxonomy($rank{subspecies}{taxid});
    if(!$mpSubspeciesTaxid){
      $mpSubspeciesTaxid = $MP->addTaxonomy($rank{subspecies}{taxid}, $rank{genus}{text}, $rank{species}{text}, $rank{subspecies}{text});
    }
    my $taxid = $rank{subspecies}{taxid} || $rank{species}{taxid} || $rank{genus}{taxid}
                || die "ERROR: no taxonomy id was found for $biosample_acc";

    return $taxid;

    $retstart++;
  }
  
  die "Internal error getting taxonomy ID";
}

sub getBiosampleInfo{
  my($biosample_acc, $settings) = @_;

  my @biosampleInfo;
  
  my $esearch = Bio::DB::EUtilities->new(-eutil  => "esearch",
                                         -db     => "biosample",
                                         -term   => "$biosample_acc",
                                         -email  => $email,
                                         -usehistory => "y",
                                        );
  my $count = $esearch->get_count;
  my $hist = $esearch->next_History || die "No history data returned";
  $esearch->set_parameters(-eutil   => "efetch",
                           -rettype => "xml",
                           -history => $hist,
                          );
  
  my($retry, $retmax, $retstart) = (0, 500, 0);
  RETRIEVE_BIOSAMPLE:
  while($retstart < $count){
    $esearch->set_parameters(-retmax   => $retmax,
                             -retstart => $retstart);
    my $xml = "";
    eval{
      $esearch->get_Response(-cb=>sub{
          my($data) = @_;
          $xml .= $data;
      });
    };
    if($@ || !$xml){
      die "Server error: $@. Try again later" if($retry==5);
      logmsg "Server error, redo #$retry";
      $retry++ && next RETRIEVE_BIOSAMPLE;
    }

    my $mashpitBiosample = parseBiosampleXml($xml);

    #print Dumper $biosampleInfo,$retry, $retstart, $count,'===';
    push(@biosampleInfo, $mashpitBiosample);
    $retstart++;
  }

  return \@biosampleInfo;
}

# Parse biosample information for our sql schema
# Needed: biosample_acc, strain, isolate, taxid, collected_by,
# collection_date, latitude, longitude, host_taxid, host_disease,
# isolation_source, serovar
sub parseBiosampleXml{
  my($xml) = @_;
  my $xmlConverter = XML::Hash->new();
  my $biosampleHash = $xmlConverter->fromXMLStringtoHash($xml);

  my %h = ();
  
  if(defined($$biosampleHash{BioSampleSet}{BioSample}{Attributes}{Attribute})){
    my $attribute = $$biosampleHash{BioSampleSet}{BioSample}{Attributes}{Attribute};
    for my $a(@$attribute){
      $h{$$a{attribute_name}} = $$a{text};
    }
  }
  
  my $organism = $$biosampleHash{BioSampleSet}{BioSample}{Description}{Organism};
  $h{taxid} = $$organism{taxonomy_id};
  $h{latitude}  = $h{lat_lon}; # TODO parse correctly
  $h{longitude} = $h{lat_lon}; # TODO parse correctly
  $h{host_taxid}= 'missing';   # TODO
  $h{host_disease}='missing';  # TODO

  return \%h;
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


