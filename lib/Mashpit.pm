#!/usr/bin/env perl
package Mashpit;
use strict;
use warnings;
use Exporter qw(import);
use File::Basename qw/fileparse basename dirname/;
use Data::Dumper;
use Bio::Kmer;
use DBI;
use Carp qw/croak carp/;

my $here = dirname($INC{"Mashpit.pm"});
use lib dirname($INC{"Mashpit.pm"});

our @EXPORT_OK = qw(
           @fastqExt @fastaExt @richseqExt @mshExt
           $MASHPIT_VERSION
         );

local $0=basename $0;

######
# CONSTANTS

our $VERSION = "0.3";
our $MASHPIT_VERSION=$VERSION;
our @fastqExt=qw(.fastq.gz .fastq .fq .fq.gz);
our @fastaExt=qw(.fasta .fna .faa .mfa .fas .fsa .fa);
our @mshExt=qw(.msh);
# Richseq extensions were obtained mostly from bioperl under
# the genbank, embl, and swissprot entries, under
# the source for Bio::SeqIO
our @richseqExt=qw(.gb .gbank .genbank .gbk .gbs .gbf .embl .ebl .emb .dat .swiss .sp);

sub logmsg{print STDERR "$0: @_\n";}
# Must supply the path as the first argument
sub new{
  my($class,$dir,$settings)=@_;

  my $dbFile = "$dir/db.sqlite";
  if(! -d $dir){
    mkdir ($dir) 
      or croak("ERROR: could not make directory $dir: $!");
  }

  my $self={
    dir      => $dir,   # path of mashpit directory
    dbFile   => $dbFile,# database path
    dbh      => undef,  # DBI object for $dbFile
  };
  bless($self,$class);

  $self->selectDb($dbFile);
  return $self;
}

# Create an SQLite database for genome distances.
sub selectDb{
  my($self, $dbFile)=@_;

  $self->{dbFile}=$dbFile;

  $self->connect();

  if(-e $dbFile && -s $dbFile > 0){
    return 0;
  }

  # grab the CREATE TABLE statements
  my $CREATE = "";
  open(my $fh, "<", "$here/mashpit_init.sql") or croak "ERROR: could not open $here/mashpit_init.sql: $!";
  while(<$fh>){
    $CREATE.=$_;
  }
  close $fh;

  my $dbh=$self->{dbh};
  $dbh->begin_work;
  # Separate every command by semicolon
  while($CREATE =~ /(.+?;)/gs){
    my $sth = $dbh->prepare($1)
      or die "ERROR: $DBI::errstr";
    $sth->execute()
      or die "ERROR: $DBI::errstr";
  }
  $dbh->commit;

  return 1;
}

sub connect{
  my($self)=@_;

  my $dbFile=$self->{dbFile};
  my $dbh=DBI->connect("dbi:SQLite:dbname=$dbFile","","",{
      RaiseError => 1,
  });
  $dbh->do("PRAGMA foreign_keys = ON");
  
  $self->{dbh}=$dbh;
  
  return $dbh;
}

