#!/usr/bin/env perl
use strict;
use warnings;
use File::Temp qw/tempdir/;
use File::Copy qw/mv/;
use File::Basename qw/basename/;
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
local $0 = basename $0;
exit(main());

sub main{
  my $settings={};
  GetOptions($settings, qw(help tempdir=s)) or die $!;
  usage() if($$settings{help} || @ARGV < 2);
  $$settings{tempdir} ||= tempdir("MASHPIT.XXXXXX", TMPDIR=>1, CLEANUP=>1);
  
  my($dbPath, $biosample_acc) = @ARGV;

  die "ERROR: database does not exist: $dbPath" if(!-e $dbPath);
  $MP = Mashpit->new($dbPath);

  # Taxonomy info
  my $taxid         = getTaxonInfo($biosample_acc, $settings);
  # Biosample info
  my $biosampleInfo = getBiosampleInfo($biosample_acc, $settings);
  # SRA info
  my $srrInfo       = getSraInfo($biosample_acc, $settings);

  my $sketchPath    = sketchGenome($dbPath, $taxid, $srrInfo, $biosampleInfo, $settings);
  die "TODO add to db";

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

  # Does it already exist?  If so, just return.
  my $dbBiosample = $MP->getBiosample($biosample_acc);
  if(defined($dbBiosample)){
    die Dumper $dbBiosample, "Line ".__LINE__;
    return $dbBiosample;
  }

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
    my $tmp = $MP->addBiosample($mashpitBiosample);
    die Dumper $tmp;

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
  $h{biosample_acc}=$$biosampleHash{BioSampleSet}{BioSample}{accession};

  # Convert keys to lowercase
  my %lcH;
  while(my($key,$value) = each(%h)){
    $lcH{lc($key)} = $value;
  }

  # Get rid of some fields
  for(qw(genus species subspecies geo_loc_name lat_lon)){
    delete($lcH{$_});
  }

  return \%lcH;
}

sub getSraInfo{
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
                                       -db      => "sra",
                                       -target  => "sra",
                                       -dbFrom  => "biosample",
                                       -usehistory=>"y",
                                       -email   => $email,
                                       -id      => \@id,
                                      );
  my @srr = uniq $elink->get_ids;

  if(@srr > 1){
    logmsg "WARNING: found more than one run id for $biosample_acc. Will only use $srr[0].";
    logmsg "All run ids found: ". join(" ", @srr);
  }

  my $efetch = Bio::DB::EUtilities->new(-eutil   => "efetch",
                                        -db      => "sra",
                                        -email   => $email,
                                        -id      => \@srr,
                                        -rettype => "xml",
                                       );
  my $count = $esearch->get_count;
  my($retry, $retmax, $retstart) = (0, 500, 0);
  RETRIEVE_SRA:
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
      $retry++ && next RETRIEVE_SRA;
    }

    my $xmlConverter = XML::Hash->new();
    my $sraHash = $xmlConverter->fromXMLStringtoHash($xml);
    my $runHash = $$sraHash{EXPERIMENT_PACKAGE_SET}{EXPERIMENT_PACKAGE}{RUN_SET}{RUN};

    #die Dumper $runHash;
    my $accession = $$runHash{accession};
    my($asmUrl, $rawUrl, $asmPath, $rawPath);
    for my $remoteFile(@{ $$runHash{SRAFiles}{SRAFile} }){
      if($$remoteFile{filename} eq "$accession.realign"){
        $asmUrl = $$remoteFile{url};
        $asmPath = "$$settings{tempdir}/$biosample_acc.$accession.fasta";
        system("dump-ref-fasta $asmUrl > $asmPath");
        die "ERROR dumping $$remoteFile{url}" if $?;
      } elsif ($$remoteFile{filename} eq $accession){
        $rawUrl = $$remoteFile{url};
      }
    }
    return {accession=>$accession, asmUrl=>$asmUrl, rawUrl=>$rawUrl,
            asmPath=>$asmPath, rawPath=>$rawPath};
    
    $retstart++;
  }

  die "ERROR: no SRA results for $biosample_acc";
}

# srrInfo has keys:
#  asmPath, asmUrl,
#  accession,
#  rawPath, rawUrl
sub sketchGenome{
  my($dbPath, $taxid, $srrInfo, $biosampleInfo, $settings) = @_;

  my $sketchPath="";
  if($$srrInfo{asmPath}){
    my $basename = basename($$srrInfo{asmPath});
    my $newPath = "$dbPath/sketches/$basename";
    $sketchPath = $newPath . ".msh";
    mv($$srrInfo{asmPath}, $newPath);
    system("cd $dbPath/sketches; mash sketch $basename");
    die "ERROR sketching $$srrInfo{asmPath}" if $?;
    unlink($newPath);
  }
  
  die "MADE IT TO ".__LINE__;
}

sub usage{
  print "Usage: $0 mashpitDb biosample_acc [biosample_acc2...]
  \n";
  exit 0;
}

