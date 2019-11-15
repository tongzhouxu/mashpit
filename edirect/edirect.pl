#!/usr/bin/env perl

# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#            National Center for Biotechnology Information (NCBI)
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government do not place any restriction on its use or reproduction.
#  We would, however, appreciate having the NCBI and the author cited in
#  any work or product based on this material.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
# ===========================================================================
#
# File Name:  edirect.pl
#
# Author:  Jonathan Kans
#
# Version Creation Date:   8/20/12
#
# ==========================================================================

# Entrez Direct - EDirect

# use strict;
use warnings;

my ($LibDir, $ScriptName);

use File::Spec;

# EDirect version number

$version = "12.5";

BEGIN
{
  my $Volume;
  ($Volume, $LibDir, $ScriptName) = File::Spec->splitpath($0);
  $LibDir = File::Spec->catpath($Volume, $LibDir, '');
  if (my $RealPathname = eval {readlink $0}) {
    do {
      $RealPathname = File::Spec->rel2abs($RealPathname, $LibDir);
      ($Volume, $LibDir, undef) = File::Spec->splitpath($RealPathname);
      $LibDir = File::Spec->catpath($Volume, $LibDir, '')
    } while ($RealPathname = eval {readlink $RealPathname});
  } else {
    $LibDir = File::Spec->rel2abs($LibDir)
  }
  $LibDir .= '/aux/lib/perl5';
}
use lib $LibDir;

# usage - edirect.pl -function arguments

use Data::Dumper;
use Encode;
use Getopt::Long;
use HTML::Entities;
use JSON::PP;
use LWP::Simple;
use LWP::UserAgent;
use MIME::Base64;
use Net::FTP;
use Net::hostent;
use POSIX;
use Time::HiRes;
use URI::Escape;
use XML::Simple;

# required first argument is name of function to run

$fnc = shift or die "Must supply function name on command line\n";

if ( $fnc eq "-internal" ) {
  $fnc = shift or die "Must supply function name on command line\n";
}

# get starting time

$begin_time = Time::HiRes::time();

# definitions

use constant false => 0;
use constant true  => 1;

# URL address components

$base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/";

$ecitmat  = "ecitmatch.cgi";
$efetch   = "efetch.fcgi";
$einfo    = "einfo.fcgi";
$elink    = "elink.fcgi";
$epost    = "epost.fcgi";
$esearch  = "esearch.fcgi";
$espell   = "espell.fcgi";
$esummary = "esummary.fcgi";

# utility subroutines

sub clearflags {
  @rest = ();
  %labels = ();
  %macros = ();
  $alias = "";
  $author = "";
  $basx = "";
  $batch = false;
  $chr_start = -1;
  $chr_stop = -1;
  $cited = false;
  $cites = false;
  $class = "";
  $clean = false;
  $cmd = "";
  $compact = false;
  $complexity = 0;
  $country = "";
  $db = "";
  $dbase = "";
  $dbs = "";
  $dbto = "";
  $dcsm = false;
  $debug = false;
  $drop = false;
  $dttype = "";
  $emaddr = "";
  $email = "";
  $err = "";
  $extend = -1;
  $extrafeat = -1;
  $field = "";
  $fields = false;
  $feature = "";
  $gtype = "";
  $help = false;
  $holding = "";
  $http = "";
  $id = "";
  $input = "";
  $internal = false;
  $journal = "";
  $json = false;
  $just_num = false;
  $key = "";
  $kind = "";
  $lbl = "";
  $links = false;
  $location = "";
  $log = false;
  $max = 0;
  $meadow = "";
  $min = 0;
  $mndate = "";
  $mode = "";
  $molecule = "";
  $mxdate = "";
  $name = "";
  $neighbor = false;
  $num = "";
  $organism = "";
  $output = "";
  $page = "";
  $pair = "";
  $pipe = false;
  $pathway = "";
  $pub = "";
  $query = "";
  $raw = false;
  $related = false;
  $result = 0;
  $revcomp = false;
  $rldate = 0;
  $seq_start = 0;
  $seq_stop = 0;
  $showgi = false;
  $silent = false;
  $sort = "";
  $source = "";
  $spell = false;
  $split = "";
  $status = "";
  $stp = "";
  $stpminusone = 0;
  $strand = "";
  $style = "";
  $tool = "";
  $trim = false;
  $trunc = false;
  $tuul = "";
  $type = "";
  $verbose = false;
  $volume = "";
  $web = "";
  $released = "";
  $word = false;
  $year = "";

  $stop_words="#a#about#again#all#almost#also#although#always#among#an#and#" .
  "another#any#are#as#at#be#because#been#before#being#between#both#but#by#can#" .
  "could#did#do#does#done#due#during#each#either#enough#especially#etc#for#" .
  "found#from#further#had#has#have#having#here#how#however#i#if#in#into#is#it#" .
  "its#itself#just#kg#km#made#mainly#make#may#mg#might#ml#mm#most#mostly#" .
  "must#nearly#neither#no#nor#obtained#of#often#on#our#overall#perhaps#pmid#" .
  "quite#rather#really#regarding#seem#seen#several#should#show#showed#shown#" .
  "shows#significantly#since#so#some#such#than#that#the#their#theirs#them#" .
  "then#there#therefore#these#they#this#those#through#thus#to#upon#use#used#" .
  "using#various#very#was#we#were#what#when#which#while#with#within#without#would#";

  $os = "$^O";

  $api_key = "";
  $api_key = $ENV{NCBI_API_KEY} if defined $ENV{NCBI_API_KEY};

  $abbrv_flag = false;
  if (defined $ENV{EDIRECT_DO_AUTO_ABBREV} && $ENV{EDIRECT_DO_AUTO_ABBREV} eq "true" ) {
    $abbrv_flag = true;
  }
}

sub do_sleep {

  if ( $internal ) {
    Time::HiRes::usleep(1000);
    return;
  }

  if ( $api_key ne "" ) {
    if ( $log ) {
      print STDERR "sleeping 1/10 second\n";
    }
    Time::HiRes::usleep(110000);
    return;
  }

  if ( $log ) {
    print STDERR "sleeping 1/3 second\n";
  }
  Time::HiRes::usleep(350000);
}

# gets a live UID for any database

sub get_zero_uid {

  my $db = shift (@_);

  my $val = "";

  %zeroUidHash = (
    'annotinfo'        =>  '122134',
    'assembly'         =>  '443538',
    'biocollections'   =>  '7370',
    'bioproject'       =>  '146229',
    'biosample'        =>  '3737421',
    'biosystems'       =>  '1223165',
    'blastdbinfo'      =>  '998664',
    'books'            =>  '1371014',
    'cdd'              =>  '274590',
    'clinvar'          =>  '10510',
    'clone'            =>  '18646800',
    'dbvar'            =>  '6173073',
    'gap'              =>  '872875',
    'gapplus'          =>  '136686',
    'gds'              =>  '200022309',
    'gencoll'          =>  '398148',
    'gene'             =>  '3667',
    'genome'           =>  '52',
    'geoprofiles'      =>  '16029743',
    'grasp'            =>  '2852486',
    'gtr'              =>  '559277',
    'homologene'       =>  '510',
    'ipg'              =>  '422234',
    'medgen'           =>  '162753',
    'mesh'             =>  '68007328',
    'ncbisearch'       =>  '3158',
    'nlmcatalog'       =>  '0404511',
    'nuccore'          =>  '1322283',
    'nucleotide'       =>  '1322283',
    'omim'             =>  '176730',
    'orgtrack'         =>  '319950',
    'pcassay'          =>  '1901',
    'pccompound'       =>  '16132302',
    'pcsubstance'      =>  '126522451',
    'pmc'              =>  '209839',
    'popset'           =>  '27228303',
    'probe'            =>  '9997691',
    'protein'          =>  '4557671',
    'proteinclusters'  =>  '2945638',
    'pubmed'           =>  '2539356',
    'seqannot'         =>  '9561',
    'snp'              =>  '137853337',
    'sparcle'          =>  '10022454',
    'sra'              =>  '190091',
    'structure'        =>  '61024',
    'taxonomy'         =>  '562',
    'unigene'          =>  '1132160',
  );

  if ( defined $zeroUidHash{$db} ) {
    $val = $zeroUidHash{$db};
  }

  return $val;
}

# support for substitution of (#keyword) to full query phrase or URL component

sub map_labels {

  my $qury = shift (@_);

  if ( $query !~ /\(#/ ) {
    return $qury;
  }

  if ( scalar (keys %labels) > 0 ) {
    for ( keys %labels ) {
      $ky = $_;
      $vl = $labels{$_};
      $qury =~ s/\((#$ky)\)/\($vl#\)/g;
    }
    $qury =~ s/\((\w+)#\)/\(#$1\)/g;
  }

  return $qury;
}

sub map_macros {

  my $qury = shift (@_);

  if ( $qury !~ /\(#/ ) {
    return $qury;
  }

  if ( scalar (keys %macros) > 0 ) {
    for ( keys %macros ) {
      $ky = $_;
      $vl = $macros{$_};
      $qury =~ s/\((\#$ky)\)/$vl/g;
    }
  }

  return $qury;
}

sub read_aliases {

  if ( $alias ne "" ) {
    if (open (my $PROXY_IN, $alias)) {
      while ( $thisline = <$PROXY_IN> ) {
        $thisline =~ s/\r//;
        $thisline =~ s/\n//;
        $thisline =~ s/ +/ /g;
        $thisline =~ s/> </></g;

        if ( $thisline =~ /(.+)\t(.+)/ ) {
          $ky = $1;
          $vl = $2;
          $vl =~ s/\"//g;
          $macros{"$ky"} = "$vl";
        }
      }
      close ($PROXY_IN);
    } else {
      print STDERR "Unable to open alias file '$alias'\n";
    }
  }
}

# base EUtils URL can be overridden for access to test versions of server

sub adjust_base {

  if ( $basx ne "" ) {
    $internal = false;
  }

  if ( $basx eq "" ) {

    if ( $internal ) {
      $base = "https://eutils-internal.ncbi.nlm.nih.gov/entrez/eutils/";
      return;
    }

    # if base not overridden, check URL of previous query, stick with main or colo site,
    # since history server data for EUtils does not copy between locations, by design

    if ( $web ne "" and $web =~ /NCID_\d+_\d+_(\d+)\.\d+\.\d+\.\d+_\d+_\d+_\d+/ ) {
      if ( $1 == 130 ) {
        $base = "https://eutils.be-md.ncbi.nlm.nih.gov/entrez/eutils/";
      } elsif ( $1 == 165 ) {
        $base = "https://eutils.st-va.ncbi.nlm.nih.gov/entrez/eutils/";
      }
    }
    return;
  }

  # shortcut for eutilstest base
  if ( $basx eq "test" ) {
    $basx = "https://eutilstest.ncbi.nlm.nih.gov/entrez/eutils";
  }

  if ( $basx =~ /\(#/ ) {
    $basx = map_macros ($basx);
  }

  if ( $basx !~ /^https:\/\// ) {
    $basx = "https://" . $basx;
  }

  if ( $basx !~ /\/$/ ) {
    $basx .= "/";
  }

  if ( $basx !~ /\/entrez\/eutils\/$/ ) {
    $basx .= "entrez/eutils/";
  }

  $base = $basx;
}

# ensure that ENTREZ_DIRECT data structure contains required fields

sub test_edirect {

  my $dbsx = shift (@_);
  my $webx = shift (@_);
  my $keyx = shift (@_);
  my $numx = shift (@_);
  my $label = shift (@_);

  if ( $dbsx eq "" ) {
    close (STDOUT);
    die "Db value not found in $label input\n";
  }
  if ( $webx eq "" ) {
    close (STDOUT);
    die "WebEnv value not found in $label input\n";
  }
  if ( $keyx eq "" ) {
    close (STDOUT);
    die "QueryKey value not found in $label input\n";
  }
  if ( $numx eq "" ) {
    close (STDOUT);
    die "Count value not found in $label input\n";
  }
}

sub get_email {

  # adapted from code provided by Aaron Ucko

  my $addr = "";
  if (defined $ENV{EMAIL}) {
    $addr = $ENV{EMAIL};
  } else {
    # Failing that, try to combine the username from USER or whoami
    # with the contents of /etc/mailname if available or the system's
    # qualified host name.  (Its containing domain may be a better
    # choice in many cases, but anyone contacting abusers can extract
    # it if necessary.)
    my $lhs = $ENV{USER} || `whoami`;
    my $rhs = "";
    if (-r '/etc/mailname') {
      $rhs = `cat /etc/mailname`;
    } else {
      my @uname = POSIX::uname();
      $rhs = $uname[1];
      if ($rhs !~ /\./) {
        # clearly unqualified, try to resolve back and forth
        my $h = gethostbyname($rhs);
        if (defined $h  &&  $h->name =~ /\./) {
          $rhs = $h->name;
        }
      }
    }
    $addr = $lhs . '@' . $rhs;
  }

  return $addr;
}

# elink and epost currently need a separate ESearch to get the correct result count

sub get_count {

  my $dbsx = shift (@_);
  my $webx = shift (@_);
  my $keyx = shift (@_);
  my $tulx = shift (@_);
  my $emlx = shift (@_);

  test_edirect ( $dbsx, $webx, $keyx, "1", "count" );

  $url = $base . $esearch;
  $url .= "?db=$dbsx&query_key=$keyx&WebEnv=$webx";
  $url .= "&retmax=0&usehistory=y";

  $url .= "&edirect=$version";

  if ( $os ne "" ) {
    $url .= "&edirect_os=$os";
  }

  if ( $api_key ne "" ) {
    $url .= "&api_key=$api_key";
  }

  if ( $tulx eq "" ) {
    $tulx = "entrez-direct";
  }
  if ( $tulx ne "" ) {
    $tulx =~ s/\s+/\+/g;
    $url .= "&tool=$tulx-count";
  }

  if ( $emlx eq "" ) {
    $emlx = get_email ();
  }
  if ( $emlx ne "" ) {
    $emlx =~ s/\s+/\+/g;
    $url .= "&email=$emlx";
  }

  $keyx = "";
  $numx = "";
  $errx = "";

  do_sleep ();

  $output = get ($url);

  if ( ! defined $output ) {
    print STDERR "Failure of get_count '$url'\n";
    $result = 1;
    return "", "";
  }

  if ( $output eq "" ) {
    print STDERR "No get_count output returned from '$url'\n";
    $result = 1;
    return "", ""
  }

  if ( $debug ) {
    print STDERR "$output\n";
  }

  $keyx = $1 if ($output =~ /<QueryKey>(\S+)<\/QueryKey>/);
  $numx = $1 if ($output =~ /<Count>(\S+)<\/Count>/);
  $errx = $1 if ($output =~ /<Error>(.+?)<\/Error>/i);

  if ( $errx ne "" ) {
    close (STDOUT);
    $result = 1;
    die "ERROR in count output: $errx\nURL: $url\n\n";
  }

  if ( $numx eq "" ) {
    $result = 1;
    die "Count value not found in count output - WebEnv $webx\n";
  }

  return $numx, $keyx;
}

sub get_uids {

  my @working = ();

  my $dbsx = shift (@_);
  my $webx = shift (@_);
  my $keyx = shift (@_);
  my $sttx = shift (@_);
  my $chkx = shift (@_);
  my $numx = shift (@_);
  my $tulx = shift (@_);
  my $emlx = shift (@_);

  test_edirect ( $dbsx, $webx, $keyx, $numx, "uids" );

  # adjust retmax if -stop has been overridden
  # e.g., efetch -format uid -start 2 -stop 9

  if ( $sttx + $chkx > $numx ) {
    $chkx = $numx - $sttx;
  }

  $url = $base . $esearch . "?db=$dbsx&query_key=$keyx&WebEnv=$webx";
  $url .= "&rettype=uilist&retmode=text";
  $url .= "&retstart=$sttx&retmax=$chkx";

  $url .= "&edirect=$version";

  if ( $os ne "" ) {
    $url .= "&edirect_os=$os";
  }

  if ( $api_key ne "" ) {
    $url .= "&api_key=$api_key";
  }

  if ( $tulx eq "" ) {
    $tulx = "edirect";
  }
  if ( $tulx ne "" ) {
    $tulx =~ s/\s+/\+/g;
    $url .= "&tool=$tulx";
  }

  if ( $emlx eq "" ) {
    $emlx = get_email ();
  }
  if ( $emlx ne "" ) {
    $emlx =~ s/\s+/\+/g;
    $url .= "&email=$emlx";
  }

  if ( $debug ) {
    print STDERR "$url\n";
  }

  $keep_trying = true;
  for ( $try = 0; $try < 3 && $keep_trying; $try++) {

    do_sleep ();

    $data = get ($url);

    if ( defined $data ) {
      $keep_trying = false;
    } else {
      print STDERR "Failure of get_uids '$url'\n";
      $result = 1;
    }
  }
  if ( $keep_trying ) {
    return @working;
  }

  if ( $data eq "" ) {
    print STDERR "No get_uids output returned from '$url'\n";
    return @working;
  }

  if ( $debug ) {
    print STDERR "$data\n";
  }

  my @ids = ($data =~ /<Id>(\d+)<\/Id>/g);
  foreach $uid (@ids) {
    push (@working, $uid);
  }

  return @working;
}

# send actual query

sub do_post_yielding_ref {

  my $urlx = shift (@_);
  my $argx = shift (@_);
  my $tulx = shift (@_);
  my $emlx = shift (@_);
  my $intr = shift (@_);

  if ( $os ne "" ) {
    $argx .= "&edirect_os=$os";
  }

  if ( $api_key ne "" ) {
    $argx .= "&api_key=$api_key";
  }

  $argx .= "&edirect=$version";

  if ( $intr ) {
    if ( $tulx eq "" ) {
      $tulx = "edirect";
    }
    if ( $tulx ne "" ) {
      $tulx =~ s/\s+/\+/g;
      $argx .= "&tool=$tulx";
    }

    if ( $emlx eq "" ) {
      $emlx = get_email ();
    }
    if ( $emlx ne "" ) {
      $emlx =~ s/\s+/\+/g;
      $argx .= "&email=$emlx";
    }
  }

  my $empty = '';
  $rslt = \$empty;

  if ( $debug or $log ) {
    print STDERR "$urlx?$argx\n";
  }

  if ( $http eq "get" or $http eq "GET" ) {
    if ( $argx ne "" ) {
      $urlx .= "?";
      $urlx .= "$argx";
    }

    do_sleep ();

    $rslt = get ($urlx);

    if ( ! defined $rslt ) {
      print STDERR "Failure of do_get '$urlx'\n";
      $result = 1;
      return "";
    }

    if ( $rslt eq "" ) {
      print STDERR "No do_get output returned from '$urlx'\n";
      $result = 1;
      return "";
    }

    if ( $debug ) {
      print STDERR "$rslt\n";
    }

    return \$rslt;
  }

  $usragnt = new LWP::UserAgent (timeout => 300);
  $usragnt->env_proxy;

  $req = new HTTP::Request POST => "$urlx";
  $req->content_type('application/x-www-form-urlencoded');
  $req->content("$argx");

  do_sleep ();

  $res = $usragnt->request ( $req );

  if ( $res->is_success) {
    $rslt = $res->content_ref;
  } else {
    $stts = $res->status_line;
    print STDERR $stts . "\n";
    if ( $stts eq "429 Too Many Requests" ) {
      if ( $api_key eq "" ) {
        print STDERR "PLEASE REQUEST AN API_KEY FROM NCBI\n";
      } else {
        print STDERR "TOO MANY REQUESTS EVEN WITH API_KEY\n";
      }
    }
  }

  if ( $$rslt eq "" ) {
    if ( $argx ne "" ) {
      $urlx .= "?";
      $urlx .= "$argx";
    }

    do_sleep ();

    print STDERR "No do_post output returned from '$urlx'\n";
    print STDERR "Result of do_post http request is\n";
    print STDERR Dumper($res);
    print STDERR "\n";
  }

  if ( $debug ) {
    print STDERR $$rslt, "\n";
  }

  return $rslt;
}

sub do_post {
  my $rslt = do_post_yielding_ref(@_);
  return $$rslt;
}

# read ENTREZ_DIRECT data structure

sub read_edirect {

  my $dbsx = "";
  my $webx = "";
  my $keyx = "";
  my $numx = "";
  my $stpx = "";
  my $errx = "";
  my $tulx = "";
  my $emlx = "";
  my $inlabel = false;
  my $inmacro = false;
  my $ky = "";
  my $vl = "";

  @other = ();
  $has_num = false;
  $all_num = true;

  while ( defined($thisline = <STDIN>) ) {
    $thisline =~ s/\r//;
    $thisline =~ s/\n//;
    $thisline =~ s/^\s+//;
    $thisline =~ s/\s+$//;
    if ( $thisline =~ /<Labels>/ ) {
      $inlabel = true;
    } elsif ( $thisline =~ /<\/Labels>/ ) {
      $inlabel = false;
    } elsif ( $thisline =~ /<Macros>/ ) {
      $inmacro = true;
    } elsif ( $thisline =~ /<\/Macros>/ ) {
      $inmacro = false;
    } elsif ( $inlabel ) {
      if ( $thisline =~ /<Label>/ ) {
        $ky = "";
        $vl = "";
      } elsif ( $thisline =~ /<\/Label>/ ) {
        if ( $vl ne "" and $ky ne "" ) {
          $labels{"$ky"} = "$vl";
        }
      } elsif ( $thisline =~ /<Key>(.+)<\/Key>/ ) {
        $ky = $1;
      } elsif ( $thisline =~ /<Val>(.+)<\/Val>/ ) {
        $vl = $1;
      }
    } elsif ( $inmacro ) {
      if ( $thisline =~ /<Macro>/ ) {
        $ky = "";
        $vl = "";
      } elsif ( $thisline =~ /<\/Macro>/ ) {
        if ( $vl ne "" and $ky ne "" ) {
          $macros{"$ky"} = "$vl";
        }
      } elsif ( $thisline =~ /<Key>(.+)<\/Key>/ ) {
        $ky = $1;
      } elsif ( $thisline =~ /<Val>(.+)<\/Val>/ ) {
        $vl = $1;
      }
    }
    if ( $thisline =~ /<Db>(\S+)<\/Db>/ ) {
      $dbsx = $1;
    }
    if ( $thisline =~ /<WebEnv>(\S+)<\/WebEnv>/ ) {
      $webx = $1;
    }
    if ( $thisline =~ /<QueryKey>(\S+)<\/QueryKey>/ ) {
      $keyx = $1;
    }
    if ( $thisline =~ /<Count>(\S+)<\/Count>/ ) {
      $numx = $1;
    }
    if ( $thisline =~ /<Step>(\S+)<\/Step>/ ) {
      $stpx = $1 + 1;
    }
    if ( $thisline =~ /<Error>(.+?)<\/Error>/i ) {
      $errx = $1;
    }
    if ( $thisline =~ /<Tool>(.+?)<\/Tool>/i ) {
      $tulx = $1;
    }
    if ( $thisline =~ /<Email>(.+?)<\/Email>/i ) {
      $emlx = $1;
    }
    if ( $thisline =~ /<Silent>Y<\/Silent>/i ) {
      $silent = true;
    }
    if ( $thisline =~ /<Silent>N<\/Silent>/i ) {
      $silent = false;
    } elsif ( $thisline =~ /<Verbose>Y<\/Verbose>/i ) {
      $verbose = true;
    }
    if ( $thisline =~ /<Verbose>N<\/Verbose>/i ) {
      $verbose = false;
    }
    if ( $thisline =~ /<Debug>Y<\/Debug>/i ) {
      $debug = true;
      $silent = false;
    }
    if ( $thisline =~ /<Debug>N<\/Debug>/i ) {
      $debug = false;
    }
    if ( $thisline =~ /<Log>Y<\/Log>/i ) {
      $log = true;
    }
    if ( $thisline =~ /<Log>N<\/Log>/i ) {
      $log = false;
    }
    if ( $thisline =~ /<ENTREZ_DIRECT>/i ) {
    } elsif ( $thisline =~ /<\/ENTREZ_DIRECT>/i ) {
    } elsif ( $thisline =~ /<Labels>/i ) {
    } elsif ( $thisline =~ /<\/Labels>/i ) {
    } elsif ( $thisline =~ /<Macros>/i ) {
    } elsif ( $thisline =~ /<\/Macros>/i ) {
    } elsif ( $thisline =~ /^(\d+)$/ ) {
      push (@other, $1);
      $has_num = true;
    } elsif ( $thisline =~ /^(\d+).\d+$/ ) {
      push (@other, $1);
      $has_num = true;
    } elsif ( $thisline =~ /^(.+)$/ ) {
      push (@other, $1);
      $all_num = false;
    }
  }

  return ( $dbsx, $webx, $keyx, $numx, $stpx, $errx,
           $tulx, $emlx, $has_num && $all_num, @other );
}

# write ENTREZ_DIRECT data structure

sub write_edirect {

  my $dbsx = shift (@_);
  my $webx = shift (@_);
  my $keyx = shift (@_);
  my $numx = shift (@_);
  my $stpx = shift (@_);
  my $errx = shift (@_);
  my $tulx = shift (@_);
  my $emlx = shift (@_);

  my $seconds = "";
  my $end_time = Time::HiRes::time();
  my $elapsed = $end_time - $begin_time;
  if ( $elapsed > 0.0005 ) {
    if ( $elapsed =~ /(\d+)\.(\d\d\d)\d+/ ) {
      $seconds = "$1.$2";
    }
  }

  if ( $stpx eq "" and $webx ne "" ) {
    $stpx = "1";
  }

  if ( $tulx ne "" ) {
    $tulx =~ s/\s+/\+/g;
  }
  if ( $emlx ne "" ) {
    $emlx =~ s/\s+/\+/g;
  }

  if ( $verbose ) {
    print STDERR "\n";
    print STDERR "edirutil";

    if ( $dbsx ne "" ) {
      print STDERR " -db $dbsx";
    }
    if ( $webx ne "" ) {
      print STDERR " -web $webx";
    }
    if ( $keyx ne "" ) {
      print STDERR " -key $keyx";
    }
    if ( $numx ne "" ) {
      print STDERR " -count $numx";
    }
    if ( $stpx ne "" ) {
      print STDERR " -step $stpx";
    }
    if ( $seconds ne "" ) {
      print STDERR " -seconds $seconds";
    }

    print STDERR "\n\n";
  }

  if ( false ) {
    print STDERR "\n";
    print STDERR "<ENTREZ_DIRECT>\n";

    if ( $dbsx ne "" ) {
      print STDERR "  <Db>$dbsx</Db>\n";
    }
    if ( $webx ne "" ) {
      print STDERR "  <WebEnv>$webx</WebEnv>\n";
    }
    if ( $keyx ne "" ) {
      print STDERR "  <QueryKey>$keyx</QueryKey>\n";
    }
    if ( $numx ne "" ) {
      print STDERR "  <Count>$numx</Count>\n";
    }
    if ( $stpx ne "" ) {
      print STDERR "  <Step>$stpx</Step>\n";
    }
    if ( $errx ne "" ) {
      print STDERR "  <Error>$errx</Error>\n";
    }
    if ( $tulx ne "" ) {
      print STDERR "  <Tool>$tulx</Tool>\n";
    }
    if ( $emlx ne "" ) {
      print STDERR "  <Email>$emlx</Email>\n";
    }

    print STDERR "</ENTREZ_DIRECT>\n";
    print STDERR "\n";
  }

  if ( $compact ) {
    print STDERR "<ENTREZ_DIRECT>";

    if ( $dbsx ne "" ) {
      print STDERR " <Db>$dbsx</Db>";
    }
    if ( $webx ne "" ) {
      print STDERR " <WebEnv>$webx</WebEnv>";
    }
    if ( $keyx ne "" ) {
      print STDERR " <QueryKey>$keyx</QueryKey>";
    }
    if ( $numx ne "" ) {
      print STDERR " <Count>$numx</Count>";
    }

    print STDERR " </ENTREZ_DIRECT>\n";
  }

  print "<ENTREZ_DIRECT>\n";

  if ( $dbsx ne "" ) {
    print "  <Db>$dbsx</Db>\n";
  }
  if ( $webx ne "" ) {
    print "  <WebEnv>$webx</WebEnv>\n";
  }
  if ( $keyx ne "" ) {
    print "  <QueryKey>$keyx</QueryKey>\n";
  }
  if ( $numx ne "" ) {
    print "  <Count>$numx</Count>\n";
  }
  if ( $stpx ne "" ) {
    print "  <Step>$stpx</Step>\n";
  }
  if ( $errx ne "" ) {
    print "  <Error>$errx</Error>\n";
  }
  if ( $tulx ne "" ) {
    print "  <Tool>$tulx</Tool>\n";
  }
  if ( $emlx ne "" ) {
    print "  <Email>$emlx</Email>\n";
  }
  if ( $silent ) {
    print "  <Silent>Y</Silent>\n";
  }
  if ( $verbose ) {
    print "  <Verbose>Y</Verbose>\n";
  }
  if ( $debug ) {
    print "  <Debug>Y</Debug>\n";
  }
  if ( $log ) {
    print "  <Log>Y</Log>\n";
  }
  if ( scalar (keys %labels) > 0 ) {
    print "  <Labels>\n";
    for ( keys %labels ) {
      print "    <Label>\n";
      print "      <Key>$_</Key>\n";
      print "      <Val>$labels{$_}</Val>\n";
      print "    </Label>\n";
    }
    print "  </Labels>\n";
  }
  if ( scalar (keys %macros) > 0 ) {
    print "  <Macros>\n";
    for ( keys %macros ) {
      print "    <Macro>\n";
      print "      <Key>$_</Key>\n";
      print "      <Val>$macros{$_}</Val>\n";
      print "    </Macro>\n";
    }
    print "  </Macros>\n";
  }

  print "</ENTREZ_DIRECT>\n";
}

# wrapper to detect command line errors

my $abbrev_help = qq{
  To enable argument auto abbreviation resolution, run:

    export EDIRECT_DO_AUTO_ABBREV="true"

  in the terminal, or add that line to your .bash_profile configuration file.

};

sub MyGetOptions {

  my $help_msg = shift @_;

  if ( $abbrv_flag ) {
    Getopt::Long::Configure("auto_abbrev");
  } else {
    Getopt::Long::Configure("no_auto_abbrev");
  }

  if ( !GetOptions(@_) ) {
    if ( $abbrv_flag ) {
      die $help_msg;
    } else {
      print $help_msg;
      die $abbrev_help;
    }
  } elsif (@ARGV) {
    die ("Entrez Direct does not support positional arguments.\n"
         . "Please remember to quote parameter values containing\n"
         . "whitespace or shell metacharacters.\n");
  }
}

# subroutines for each -function

# ecntc prepares the requested tool and email arguments for an EUtils pipe

my $cntc_help = qq{
  -email    Contact person's address
  -tool     Name of script or program

};

sub ecntc {

  # ... | edirect.pl -contact -email darwin@beagle.edu -tool edirect_test | ...

  clearflags ();

  MyGetOptions(
    $cntc_help,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "internal" => \$internal,
    "log" => \$log,
    "compact" => \$compact,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "econtact $version\n";
    print $cntc_help;
    return;
  }

  if ( -t STDIN and not @ARGV ) {
  } else {
    ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();
  }

  read_aliases ();
  adjust_base ();

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }

  write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );
}

# efilt performs ESearch query refinement on the EUtils history server

my $filt_help = qq{
Query Specification

  -query       Query string

Document Order

  -sort        Result presentation order

Date Constraint

  -days        Number of days in the past
  -datetype    Date field abbreviation
  -mindate     Start of date range
  -maxdate     End of date range

Limit by Field

  -field       Query words individually in field
  -pairs       Query overlapping word pairs

Spell Check

  -spell       Correct misspellings in query

Publication Filters

  -pub         abstract, clinical, english, free, historical,
               journal, medline, preprint, published, review,
               structured
  -journal     pnas, "j bacteriol", ...
  -released    last_week, last_month, last_year, prev_years

Sequence Filters

  -country     usa:minnesota, united_kingdom, "pacific ocean", ...
  -feature     gene, mrna, cds, mat_peptide, ...
  -location    mitochondrion, chloroplast, plasmid, plastid
  -molecule    genomic, mrna, trna, rrna, ncrna
  -organism    animals, archaea, bacteria, eukaryotes, fungi,
               human, insects, mammals, plants, prokaryotes,
               protists, rodents, viruses
  -source      genbank, insd, pdb, pir, refseq, swissprot, tpa

Gene Filters

  -status      alive
  -type        coding, pseudo

SNP Filters

  -class       acceptor, donor, frameshift, indel, intron,
               missense, nonsense, synonymous

Biosystems Filters

  -kind        pathway
  -pathway     reactome, wikipathways

Miscellaneous Arguments

  -label       Alias for query step

};

sub process_extras {

  my $frst = shift (@_);
  my $publ = shift (@_);
  my $rlsd = shift (@_);
  my $jrnl = shift (@_);
  my $ctry = shift (@_);
  my $fkey = shift (@_);
  my $locn = shift (@_);
  my $bmol = shift (@_);
  my $orgn = shift (@_);
  my $sorc = shift (@_);
  my $stat = shift (@_);
  my $gtyp = shift (@_);
  my $clss = shift (@_);
  my $kind = shift (@_);
  my $ptwy = shift (@_);

  $publ = lc($publ);
  $rlsd = lc($rlsd);
  $jrnl = lc($jrnl);
  $ctry = lc($ctry);
  $fkey = lc($fkey);
  $bmol = lc($bmol);
  $locn = lc($locn);
  $orgn = lc($orgn);
  $sorc = lc($sorc);
  $stat = lc($stat);
  $gtyp = lc($gtyp);
  $clss = lc($clss);
  $kind = lc($kind);
  $ptwy = lc($ptwy);

  %pubHash = (
    'abstract'     =>  'has abstract [FILT]',
    'clinical'     =>  'clinical trial [FILT]',
    'english'      =>  'english [FILT]',
    'free'         =>  'freetext [FILT]',
    'historical'   =>  'historical article [FILT]',
    'journal'      =>  'journal article [FILT]',
    'last_month'   =>  'published last month [FILT]',
    'last month'   =>  'published last month [FILT]',
    'last_week'    =>  'published last week [FILT]',
    'last week'    =>  'published last week [FILT]',
    'last_year'    =>  'published last year [FILT]',
    'last year'    =>  'published last year [FILT]',
    'medline'      =>  'medline [FILT]',
    'preprint'     =>  'ahead of print [FILT]',
    'review'       =>  'review [FILT]',
    'structured'   =>  'hasstructuredabstract [WORD]',
    'trial'        =>  'clinical trial [FILT]',
  );

  %releasedHash = (
    'last_month'   =>  'published last month [FILT]',
    'last month'   =>  'published last month [FILT]',
    'last_week'    =>  'published last week [FILT]',
    'last week'    =>  'published last week [FILT]',
    'last_year'    =>  'published last year [FILT]',
    'last year'    =>  'published last year [FILT]',
  );

  @featureArray = (
    "-10_signal",
    "-35_signal",
    "3'clip",
    "3'utr",
    "5'clip",
    "5'utr",
    "allele",
    "assembly_gap",
    "attenuator",
    "c_region",
    "caat_signal",
    "cds",
    "centromere",
    "conflict",
    "d_segment",
    "d-loop",
    "enhancer",
    "exon",
    "gap",
    "gc_signal",
    "gene",
    "idna",
    "intron",
    "j_segment",
    "ltr",
    "mat_peptide",
    "misc_binding",
    "misc_difference",
    "misc_feature",
    "misc_recomb",
    "misc_rna",
    "misc_signal",
    "misc_structure",
    "mobile_element",
    "modified_base",
    "mrna",
    "mutation",
    "n_region",
    "ncrna",
    "old_sequence",
    "operon",
    "orit",
    "polya_signal",
    "polya_site",
    "precursor_rna",
    "prim_transcript",
    "primer_bind",
    "promoter",
    "propeptide",
    "protein_bind",
    "rbs",
    "regulatory",
    "rep_origin",
    "repeat_region",
    "repeat_unit",
    "rrna",
    "s_region",
    "satellite",
    "scrna",
    "sig_peptide",
    "snorna",
    "snrna",
    "source",
    "stem_loop",
    "sts",
    "tata_signal",
    "telomere",
    "terminator",
    "tmrna",
    "transit_peptide",
    "trna",
    "unsure",
    "v_region",
    "v_segment",
    "variation"
  );

  %locationHash = (
    'mitochondria'   =>  'mitochondrion [FILT]',
    'mitochondrial'  =>  'mitochondrion [FILT]',
    'mitochondrion'  =>  'mitochondrion [FILT]',
    'chloroplast'    =>  'chloroplast [FILT]',
    'plasmid'        =>  'plasmid [FILT]',
    'plastid'        =>  'plastid [FILT]',
  );

  %moleculeHash = (
    'genomic'  =>  'biomol genomic [PROP]',
    'mrna'     =>  'biomol mrna [PROP]',
    'trna'     =>  'biomol trna [PROP]',
    'rrna'     =>  'biomol rrna [PROP]',
    'ncrna'    =>  'biomol ncrna [PROP]',
  );

  %organismHash = (
    'animal'           =>  'animals [FILT]',
    'animals'          =>  'animals [FILT]',
    'archaea'          =>  'archaea [FILT]',
    'archaeal'         =>  'archaea [FILT]',
    'archaean'         =>  'archaea [FILT]',
    'archaebacteria'   =>  'archaea [FILT]',
    'archaebacterial'  =>  'archaea [FILT]',
    'bacteria'         =>  'bacteria [FILT]',
    'bacterial'        =>  'bacteria [FILT]',
    'bacterium'        =>  'bacteria [FILT]',
    'eubacteria'       =>  'bacteria [FILT]',
    'eubacterial'      =>  'bacteria [FILT]',
    'eukaryota'        =>  'eukaryota [ORGN]',
    'eukaryote'        =>  'eukaryota [ORGN]',
    'eukaryotes'       =>  'eukaryota [ORGN]',
    'fungal'           =>  'fungi [FILT]',
    'fungi'            =>  'fungi [FILT]',
    'fungus'           =>  'fungi [FILT]',
    'human'            =>  'human [ORGN]',
    'humans'           =>  'human [ORGN]',
    'insect'           =>  'insecta [ORGN]',
    'insecta'          =>  'insecta [ORGN]',
    'insects'          =>  'insecta [ORGN]',
    'mammal'           =>  'mammals [FILT]',
    'mammalia'         =>  'mammals [FILT]',
    'mammalian'        =>  'mammals [FILT]',
    'mammals'          =>  'mammals [FILT]',
    'man'              =>  'human [ORGN]',
    'metaphyta'        =>  'plants [FILT]',
    'metazoa'          =>  'animals [FILT]',
    'monera'           =>  'prokaryota [ORGN]',
    'plant'            =>  'plants [FILT]',
    'plants'           =>  'plants [FILT]',
    'prokaryota'       =>  'prokaryota [ORGN]',
    'prokaryote'       =>  'prokaryota [ORGN]',
    'prokaryotes'      =>  'prokaryota [ORGN]',
    'protist'          =>  'protists [FILT]',
    'protista'         =>  'protists [FILT]',
    'protists'         =>  'protists [FILT]',
    'rodent'           =>  'rodents [ORGN]',
    'rodentia'         =>  'rodents [ORGN]',
    'rodents'          =>  'rodents [ORGN]',
    'viral'            =>  'viruses [FILT]',
    'virus'            =>  'viruses [FILT]',
    'viruses'          =>  'viruses [FILT]',
  );

  %snpHash = (
    'acceptor'    =>  'splice acceptor variant [FXN]',
    'donor'       =>  'splice donor variant [FXN]',
    'frameshift'  =>  'frameshift [FXN]',
    'indel'       =>  'cds indel [FXN]',
    'intron'      =>  'intron [FXN]',
    'missense'    =>  'missense [FXN]',
    'nonsense'    =>  'nonsense [FXN]',
    'synonymous'  =>  'synonymous codon [FXN]',
  );

  %sourceHash = (
    'ddbj'       =>  'srcdb ddbj [PROP]',
    'embl'       =>  'srcdb embl [PROP]',
    'genbank'    =>  'srcdb genbank [PROP]',
    'insd'       =>  'srcdb ddbj/embl/genbank [PROP]',
    'pdb'        =>  'srcdb pdb [PROP]',
    'pir'        =>  'srcdb pir [PROP]',
    'refseq'     =>  'srcdb refseq [PROP]',
    'swissprot'  =>  'srcdb swiss prot [PROP]',
    'tpa'        =>  'srcdb tpa ddbj/embl/genbank [PROP]',
  );

  %statusHash = (
    'alive'   =>  'alive [PROP]',
    'live'    =>  'alive [PROP]',
    'living'  =>  'alive [PROP]',
  );

  %typeHash = (
    'coding'  =>  'genetype protein coding [PROP]',
    'pseudo'  =>  'genetype pseudo [PROP]',
  );

  %kindHash = (
    'pathway'  =>  'pathway [TYPE]',
  );

  %pathwayHash = (
    'reactome'      =>  'src reactome [FILT]',
    'wikipathways'  =>  'src wikipathways [FILT]',
  );

  my @working = ();

  my $suffix1 = "";
  my $suffix2 = "";

  my $is_published = false;
  my $is_prev_year = false;

  if ( $frst ne "" ) {
    push (@working, $frst);
  }

  if ( $publ ne "" ) {
    # -pub can use comma-separated list
    my @pbs = split (',', $publ);
    foreach $pb (@pbs) {
      if ( defined $pubHash{$pb} ) {
        $val = $pubHash{$pb};
        push (@working, $val);
      } elsif ( $pb eq "published" ) {
        $is_published = true;
      } else {
        die "\nUnrecognized -pub argument '$pb', use efilter -help to see available choices\n\n";
      }
    }
  }

  if ( $rlsd ne "" ) {
    if ( defined $releasedHash{$rlsd} ) {
      $val = $releasedHash{$rlsd};
      push (@working, $val);
    } elsif ( $rlsd eq "prev_years" ) {
      $is_prev_year = true;
    } elsif ( $rlsd =~ /^\d\d\d\d$/ ) {
      $val = $rlsd . " [PDAT]";
      push (@working, $val);
    } else {
      die "\nUnrecognized -released argument '$rlsd', use efilter -help to see available choices\n\n";
    }
  }

  if ( $jrnl ne "" ) {
    $val = $jrnl . " [JOUR]";
    push (@working, $val);
  }

  if ( $ctry ne "" ) {
    $val = "country " . $ctry . " [TEXT]";
    push (@working, $val);
  }

  if ( $fkey ne "" ) {
    # -feature can use comma-separated list
    my @fts = split (',', $fkey);
    foreach $ft (@fts) {
      if ( grep( /^$ft$/, @featureArray ) ) {
        $val = $ft . " [FKEY]";
        push (@working, $val);
      } else {
        die "\nUnrecognized -feature argument '$ft', use efilter -help to see available choices\n\n";
      }
    }
  }

  if ( $locn ne "" ) {
    if ( defined $locationHash{$locn} ) {
      $val = $locationHash{$locn};
      push (@working, $val);
    } else {
      die "\nUnrecognized -location argument '$locn', use efilter -help to see available choices\n\n";
    }
  }

  if ( $bmol ne "" ) {
    if ( defined $moleculeHash{$bmol} ) {
      $val = $moleculeHash{$bmol};
      push (@working, $val);
    } else {
      die "\nUnrecognized -molecule argument '$bmol', use efilter -help to see available choices\n\n";
    }
  }

  if ( $orgn ne "" ) {
    if ( defined $organismHash{$orgn} ) {
      $val = $organismHash{$orgn};
      push (@working, $val);
    } else {
      # allow any organism
      $val = $orgn . " [ORGN]";
      push (@working, $val);
      # die "\nUnrecognized -organism argument '$orgn', use efilter -help to see available choices\n\n";
    }
  }

  if ( $sorc ne "" ) {
    if ( defined $sourceHash{$sorc} ) {
      $val = $sourceHash{$sorc};
      push (@working, $val);
    } else {
      die "\nUnrecognized -source argument '$sorc', use efilter -help to see available choices\n\n";
    }
  }

  if ( $stat ne "" ) {
    if ( defined $statusHash{$stat} ) {
      $val = $statusHash{$stat};
      push (@working, $val);
    } else {
      die "\nUnrecognized -status argument '$stat', use efilter -help to see available choices\n\n";
    }
  }

  if ( $gtyp ne "" ) {
    if ( defined $typeHash{$gtyp} ) {
      $val = $typeHash{$gtyp};
      push (@working, $val);
    } else {
      die "\nUnrecognized -type argument '$gtyp', use efilter -help to see available choices\n\n";
    }
  }

  if ( $clss ne "" ) {
    if ( defined $snpHash{$clss} ) {
      $val = $snpHash{$clss};
      push (@working, $val);
    } else {
      die "\nUnrecognized -class argument '$clss', use efilter -help to see available choices\n\n";
    }
  }

  if ( $kind ne "" ) {
    if ( defined $kindHash{$kind} ) {
      $val = $kindHash{$kind};
      push (@working, $val);
    } else {
      die "\nUnrecognized -kind argument '$kind', use efilter -help to see available choices\n\n";
    }
  }

  if ( $ptwy ne "" ) {
    if ( defined $pathwayHash{$ptwy} ) {
      $val = $pathwayHash{$ptwy};
      push (@working, $val);
    } else {
      die "\nUnrecognized -pathway argument '$ptwy', use efilter -help to see available choices\n\n";
    }
  }

  my $xtras = join (" AND ", @working);

  if ( $is_published ) {
    $xtras = $xtras . " NOT ahead of print [FILT]";
  }

  if ( $is_prev_year ) {
    $xtras = $xtras . " NOT published last year [FILT]";
  }

  return $xtras;
}

# correct misspellings in query

sub spell_check_query {

  my $db = shift (@_);
  my $qury = shift (@_);

  my $url = $base . $espell;

  my $enc = uri_escape($query);
  $arg = "db=$db&term=$enc";

  my $data = do_post ($url, $arg, $tool, $email, true);

  Encode::_utf8_on($data);

  $qury = $1 if ( $data =~ /<CorrectedQuery>(.+)<\/CorrectedQuery>/ );

  return $qury;
}

sub efilt {

  # ... | edirect.pl -filter -query "bacteria [ORGN]" -days 365 | ...

  clearflags ();

  MyGetOptions(
    $filt_help,
    "query=s" => \$query,
    "q=s" => \$query,
    "sort=s" => \$sort,
    "days=i" => \$rldate,
    "mindate=s" => \$mndate,
    "maxdate=s" => \$mxdate,
    "datetype=s" => \$dttype,
    "label=s" => \$lbl,
    "field=s" => \$field,
    "spell" => \$spell,
    "pairs=s" => \$pair,
    "journal=s" => \$journal,
    "pub=s" => \$pub,
    "released=s" => \$released,
    "country=s" => \$country,
    "feature=s" => \$feature,
    "location=s" => \$location,
    "molecule=s" => \$molecule,
    "organism=s" => \$organism,
    "source=s" => \$source,
    "status=s" => \$status,
    "type=s" => \$gtype,
    "class=s" => \$class,
    "kind=s" => \$kind,
    "pathway=s" => \$pathway,
    "api_key=s" => \$api_key,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "internal" => \$internal,
    "log" => \$log,
    "compact" => \$compact,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "efilter $version\n";
    print $filt_help;
    return;
  }

  # process special filter flags, add to query string
  $query = process_extras ( $query, $pub, $released, $journal, $country, $feature, $location, $molecule, $organism, $source, $status, $gtype, $class, $kind, $pathway );

  if ( -t STDIN ) {
    if ( $query eq "" ) {
      die "Must supply -query or -days or -mindate and -maxdate arguments on command line\n";
    }
    print "efilter -query \"$query\"\n";
    return;
  }

  ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();

  read_aliases ();
  adjust_base ();

  if ( $err ne "" ) {
    write_edirect ( "", "", "", "", "", $err, "", "" );
    close (STDOUT);
    if ( ! $silent ) {
      die "ERROR in filt input: $err\n\n";
    }
    return;
  }

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }

  if ( $query eq "" && $sort eq "" && $rldate < 1 and $mndate eq "" and $mxdate eq "" ) {
    die "Must supply -query or -sort or -days or -mindate and -maxdate arguments on command line\n";
  }

  binmode STDOUT, ":utf8";

  if ( $dbase ne "" and $web ne "" and $key eq "" and $num eq "0" ) {
    write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );
    close (STDOUT);
    die "QueryKey value not found in filter input\n";
    return;
  }

  # warn on mismatch between filter argument and database
  if ( $dbase ne "pubmed" ) {
    if ( $pub ne "" ) {
      print STDERR "\nUnexpected use of pubmed filter argument\n\n";
    }
  }
  if ( $dbase ne "nucleotide" and
       $dbase ne "nuccore" and
       $dbase ne "est" and
       $dbase ne "gss" and
       $dbase ne "protein" ) {
    if ( $feature ne "" or
         $location ne "" or
         $molecule ne "" or
         $organism ne "" or
         $source ne "" ) {
      print STDERR "\nUnexpected use of sequence filter argument\n\n";
    }
  }

  test_edirect ( $dbase, $web, $key, $num, "filter" );

  # -field combines -drop and -split (-field TITL produces same behavior as Web PubMed)
  if ( $field ne "" ) {
    $query = remove_stop_words ($query);
    $query = field_each_word ($field, $query);
  }

  # -pairs separately fields query word pairs, breaking chain at stop words
  if ( $pair ne "" ) {
    $query = remove_punctuation ($query);
    if ( $query =~ /^ +(.+)$/ ) {
      $query = $1;
    }
    if ( $query =~ /^(.+) +$/ ) {
      $query = $1;
    }
    $query =~ s/ +/ /g;
    $query = field_each_pair ($pair, $query);
  }

  # spell check each query word
  if ( $spell ) {
    $query = spell_check_query ($dbase, $query);
  }

  $url = $base . $esearch;

  $arg = "db=$dbase&query_key=$key&WebEnv=$web";
  if ( $sort ne "" ) {
    if ( $sort eq "Relevance" ) {
      $sort = "relevance";
    }
    $arg .= "&sort=$sort";
  }
  $arg .= "&retmax=0&usehistory=y";
  if ( $query ne "" ) {
    $query = map_labels ($query);
    $query = map_macros ($query);
    $enc = uri_escape($query);
    $arg .= "&term=$enc";
  }
  if ( $rldate > 0 ) {
    $arg .= "&reldate=$rldate";
    if ( $dttype eq "" ) {
      $dttype = "PDAT";
    }
  }
  if ( $mndate ne "" and $mxdate ne "" ) {
    $arg .= "&mindate=$mndate&maxdate=$mxdate";
    if ( $dttype eq "" ) {
      $dttype = "PDAT";
    }
  }
  if ( $dttype ne "" ) {
    $arg .= "&datetype=$dttype";
  }

  $wb = $web;

  $web = "";
  $key = "";
  $num = "";
  $err = "";
  my $trn = "";

  $output = "";

  for ( $tries = 0; $tries < 3 && $web eq ""; $tries++) {
    $output = do_post ($url, $arg, $tool, $email, true);

    $err = $1 if ($output =~ /<Error>(.+?)<\/Error>/i);
    if ( $err eq "" ) {
      $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
      $key = $1 if ($output =~ /<QueryKey>(\S+)<\/QueryKey>/);
      $num = $1 if ($output =~ /<Count>(\S+)<\/Count>/);
      $trn = $1 if ($output =~ /<QueryTranslation>(.+?)<\/QueryTranslation>/i);
    } else {
      if ( ! $silent ) {
        print STDERR "Retrying efilter, step $stp: $err\n";
      }
      sleep 3;
    }
  }

  if ( $err ne "" ) {
    write_edirect ( "", "", "", "", "", $err, "", "" );
    close (STDOUT);
    if ( ! $silent) {
      die "ERROR in filt output: $err\nURL: $arg\nResult: $output\n\n";
    }
    return;
  }

  if ( $web eq "" ) {
    die "WebEnv value not found in filt output - WebEnv1 $wb\n";
  }
  if ( $key eq "" ) {
    die "QueryKey value not found in filt output - WebEnv1 $wb\n";
  }

  if ( $web ne $wb ) {
    $err = "WebEnv mismatch in filt output - WebEnv1 $wb, WebEnv2 $web";
    write_edirect ( "", $wb, "", "", "", $err, "", "" );
    close (STDOUT);
    die "WebEnv value changed in filt output - WebEnv1 $wb\nWebEnv2 $web\n";
  }

  if ( $num eq "" ) {
    die "Count value not found in filt output - WebEnv1 $wb\n";
  }

  if ( $lbl ne "" and $key ne "" ) {
    $labels{"$lbl"} = "$key";
  }

  write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );

  if ( $log ) {
    if ( $trn ne "" ) {
      print STDERR "$trn\n";
    }
  }
}

# efetch -format docsum calls esmry to retrieve document summaries

sub fix_mixed_content {

  my $x = shift (@_);

  while ( $x =~ /\&amp\;/ || $x =~ /\&lt\;/ || $x =~ /\&gt\;/ ) {
    HTML::Entities::decode_entities($x);
  }
  # removed mixed content tags
  $x =~ s|<b>||g;
  $x =~ s|<i>||g;
  $x =~ s|<u>||g;
  $x =~ s|<sup>||g;
  $x =~ s|<sub>||g;
  $x =~ s|</b>||g;
  $x =~ s|</i>||g;
  $x =~ s|</u>||g;
  $x =~ s|</sup>||g;
  $x =~ s|</sub>||g;
  $x =~ s|<b/>||g;
  $x =~ s|<i/>||g;
  $x =~ s|<u/>||g;
  $x =~ s|<sup/>||g;
  $x =~ s|<sub/>||g;
  # Reencode any resulting less-than or greater-than entities to avoid breaking the XML.
  $x =~ s/</&lt;/g;
  $x =~ s/>/&gt;/g;

  return $x;
}

my %fields_to_fix = (
  'biosample' => ['SampleData'],
  'medgen'    => ['ConceptMeta'],
  'sra'       => ['ExpXml', 'Runs']
);

sub fix_one_encoding {

  my $dbase = shift (@_);
  my $data = shift (@_);

  if ( $dbase eq "pubmed" ) {
    if ( $data =~ /<Title>(.+?)<\/Title>/ ) {
      my $x = $1;
      if ( $x =~ /\&amp\;/ || $x =~ /\&lt\;/ || $x =~ /\&gt\;/ || $x =~ /\</ || $x =~ /\>/ ) {
        $x = fix_mixed_content($x);
        $data =~ s/<Title>(.+?)<\/Title>/<Title>$x<\/Title>/;
      }
    }
  } elsif ( $dbase eq "gene" ) {
    if ( $data =~ /<Summary>(.+?)<\/Summary>/ ) {
      my $x = $1;
      if ( $x =~ /\&amp\;/ ) {
        HTML::Entities::decode_entities($x);
        # Reencode any resulting less-than or greater-than entities to avoid breaking the XML.
        $x =~ s/</&lt;/g;
        $x =~ s/>/&gt;/g;
        $data =~ s/<Summary>(.+?)<\/Summary>/<Summary>$x<\/Summary>/;
      }
    }
    # $data =~ s/(\s+?)<ChrStart>(\d+)<\/ChrStart>/$1<ChrStart>$2<\/ChrStart>$1<SeqStart>${\($2 + 1)}<\/SeqStart>/g;
    # $data =~ s/(\s+?)<ChrStop>(\d+)<\/ChrStop>/$1<ChrStop>$2<\/ChrStop>$1<SeqStop>${\($2 + 1)}<\/SeqStop>/g;
  } elsif ( $dbase eq "assembly" ) {
    if ( $data =~ /<Meta>(.+?)<\/Meta>/ ) {
      my $x = $1;
      if ( $x =~ /<!\[CDATA\[\s*(.+?)\s*\]\]>/ ) {
        $x = $1;
        if ( $x =~ /</ and $x =~ />/ ) {
            # If CDATA contains embedded XML, simply remove CDATA wrapper
          $data =~ s/<Meta>(.+?)<\/Meta>/<Meta>$x<\/Meta>/;
        }
      }
    }
  } elsif ( $dbase eq "gtr" ) {
    if ( $data =~ /<Extra>(.+?)<\/Extra>/ ) {
      my $x = $1;
      if ( $x =~ /<!\[CDATA\[\]\]>/ ) {
        # Remove CDATA
        $data =~ s/<Extra>(.+?)<\/Extra>/<Extra><\/Extra>/;
      }
    }
  } elsif (defined $fields_to_fix{$dbase}) {
    foreach $f (@{$fields_to_fix{$dbase}}) {
      if ( $data =~ /<$f>(.+?)<\/$f>/ ) {
        my $x = $1;
        if ( $x =~ /\&lt\;/ and $x =~ /\&gt\;/ ) {
          HTML::Entities::decode_entities($x);
          $data =~ s/<$f>(.+?)<\/$f>/<$f>$x<\/$f>/;
        }
      }
    }
  }

  return $data;
}

sub fix_bad_encoding {

  my $dbase = shift (@_);
  my $data = shift (@_);

  if ( $dbase eq "pubmed" || $dbase eq "gene" || $dbase eq "assembly" || $dbase eq "gtr" || defined $fields_to_fix{$dbase} ) {

    my @accum = ();
    my @working = ();
    my $prefix = "";
    my $suffix = "";
    my $docsumset_attrs = '';

    if ( $data =~ /(.+?)<DocumentSummarySet(\s+.+?)?>(.+)<\/DocumentSummarySet>(.+)/s ) {
      $prefix = $1;
      $docsumset_attrs = $2;
      my $docset = $3;
      $suffix = $4;

      my @vals = ($docset =~ /<DocumentSummary>(.+?)<\/DocumentSummary>/sg);
      foreach $val (@vals) {
        push (@working, "<DocumentSummary>");
        push (@working, fix_one_encoding ( $dbase, $val) );
        push (@working, "</DocumentSummary>");
      }
    }

    if ( scalar @working > 0 ) {
      push (@accum, $prefix);
      push (@accum, "<DocumentSummarySet$docsumset_attrs>");
      push (@accum, @working);
      push (@accum, "</DocumentSummarySet>");
      push (@accum, $suffix);
      $data = join ("\n", @accum);
      $data =~ s/\n\n/\n/g;
    }
  }

  return $data;
}

sub esmry {

  my $dbase = shift (@_);
  my $web = shift (@_);
  my $key = shift (@_);
  my $num = shift (@_);
  my $id = shift (@_);
  my $mode = shift (@_);
  my $min = shift (@_);
  my $max = shift (@_);
  my $tool = shift (@_);
  my $email = shift (@_);
  my $silent = shift (@_);
  my $verbose = shift (@_);
  my $debug = shift (@_);
  my $internal = shift (@_);
  my $log = shift (@_);
  my $http = shift (@_);
  my $alias = shift (@_);
  my $basx = shift (@_);

  $dbase = lc($dbase);

  if ( $dbase ne "" and $id ne "" ) {

    if ( $id eq "0" ) {

      # id "0" returns a live UID for any database

      $id = get_zero_uid ($dbase);

      if ( $id eq "0" ) {

        # id "0" is an unrecognized accession

        return;
      }
    }

    $url = $base . $esummary;

    if ( $dbase eq "pubmed" ) {
      $id =~ s/(\d+)\.\d+/$1/g;
    }
    $arg = "db=$dbase&id=$id";
    $arg .= "&version=2.0";
    if ( $mode ne "" and $mode ne "text" ) {
      $arg .= "&retmode=$mode";
    }

    $data = "";
    $retry = true;

    for ( $tries = 0; $tries < 3 && $retry; $tries++) {
      $data = do_post ($url, $arg, $tool, $email, true);

      if ($data =~ /<eSummaryResult>/i and $data =~ /<ERROR>(.+?)<\/ERROR>/i) {
        if ( ! $silent ) {
          print STDERR "Retrying esummary, step $stp: $err\n";
        }
        sleep 3;
      } else {
        $retry = false;
      }
    }

    if ($data =~ /<eSummaryResult>/i and $data =~ /<ERROR>(.+?)<\/ERROR>/i) {
      $err = $1;
      if ( $err ne "" ) {
        if ( ! $silent ) {
          print STDERR "ERROR in esummary: $err\n";
        }
      }
      if ( $err =~ "Invalid uid" ) {
        # remove Invalid uid error message from XML
        $data =~ s/<ERROR>.+?<\/ERROR>//g;
      } else {
        return;
      }
    }

    if (! $raw) {
      if ($data !~ /<Id>\d+<\/Id>/) {
        $data =~ s/<DocumentSummary uid=\"(\d+)\">/<DocumentSummary><Id>$1<\/Id>/g;
      }
      if ( $dbase eq "gtr" ) {
        $data =~ s/<DocumentSummary uid=\"\d+\">/<DocumentSummary>/g;
      }
    }

    Encode::_utf8_on($data);

    if (! $raw) {
      $data = fix_bad_encoding($dbase, $data);
    }

    # remove eSummaryResult wrapper
    $data =~ s/<!DOCTYPE eSummaryResult PUBLIC/<!DOCTYPE DocumentSummarySet PUBLIC/g;
    $data =~ s/<eSummaryResult>//g;
    $data =~ s/<\/eSummaryResult>//g;

    if ( $json ) {
      $data = xml_to_json($data);
    }

    print "$data";

    return;
  }

  if ( $dbase ne "" and $web ne "" and $key eq "" and $num eq "0" ) {
    die "QueryKey value not found in summary input\n";
    return;
  }

  test_edirect ( $dbase, $web, $key, $num, "summary" );

  if ( $stp ne "" ) {
    $stpminusone = $stp - 1;
  }

  if ( $max == 0 ) {
    if ( $silent ) {
      return;
    }
  }

  # use larger chunk for document summaries
  $chunk = 1000;
  if ( $mode eq "json" ) {
    # json has a 500 limit
    $chunk = 500;
    if ( $dbase eq "gtr" ) {
      # use smaller chunk for GTR JSON
      $chunk = 50;
    }
  }
  for ( $start = $min; $start < $max; $start += $chunk ) {
    $url = $base . $esummary;

    $chkx = $chunk;
    if ( $start + $chkx > $max ) {
      $chkx = $max - $start;
    }

    $arg = "db=$dbase&query_key=$key&WebEnv=$web";
    $arg .= "&retstart=$start&retmax=$chkx&version=2.0";
    if ( $mode ne "" and $mode ne "text" ) {
      $arg .= "&retmode=$mode";
    }

    $data = "";
    $retry = true;

    for ( $tries = 0; $tries < 3 && $retry; $tries++) {
      $data = do_post ($url, $arg, $tool, $email, true);

      if ($data =~ /<eSummaryResult>/i and $data =~ /<ERROR>(.+?)<\/ERROR>/i ) {
        if ( ! $silent ) {
          print STDERR "Retrying esummary, step $stp: $err\n";
        }
        sleep 3;
      } else {
        $retry = false;
      }
    }

    if ($data =~ /<eSummaryResult>/i and $data =~ /<ERROR>(.+?)<\/ERROR>/i) {
      $err = $1;
      if ( $err ne "" ) {
        if ( ! $silent ) {
          my $from = $start + 1;
          my $to = $start + $chunk;
          if ( $to > $max ) {
            $to = $max;
          }
          print STDERR "ERROR in esummary ($from-$to / $max): $err\n";
          print STDERR "Replicate for debugging with:\n";
          print STDERR "  edirutil -db $dbase -web $web -key $key -count $num";
          if ( $stpminusone > 0 ) {
            print STDERR " -step $stpminusone";
          }
          my $seconds = "";
          my $end_time = Time::HiRes::time();
          my $elapsed = $end_time - $begin_time;
          if ( $elapsed > 0.0005 ) {
            if ( $elapsed =~ /(\d+)\.(\d\d\d)\d+/ ) {
              $seconds = "$1.$2";
            }
          }
          if ( $seconds ne "" ) {
            print STDERR " -seconds $seconds";
          }
          print STDERR " | efetch -format docsum -start $from -stop $to\n";
        }
      }
    } else {
      if ( $verbose ) {
        my $from = $start + 1;
        my $to = $start + $chunk;
        if ( $to > $max ) {
          $to = $max;
        }
        print STDERR "( edirutil -db $dbase -web $web -key $key -count $num";
        if ( $stpminusone > 0 ) {
          print STDERR " -step $stpminusone";
        }
        my $seconds = "";
        my $end_time = Time::HiRes::time();
        my $elapsed = $end_time - $begin_time;
        if ( $elapsed > 0.0005 ) {
          if ( $elapsed =~ /(\d+)\.(\d\d\d)\d+/ ) {
            $seconds = "$1.$2";
          }
        }
        if ( $seconds ne "" ) {
          print STDERR " -seconds $seconds";
        }
        print STDERR " | efetch -format docsum -start $from -stop $to )\n";
        $begin_time = $end_time;
      }

      if (! $raw) {
        if ($data !~ /<Id>\d+<\/Id>/) {
          $data =~ s/<DocumentSummary uid=\"(\d+)\">/<DocumentSummary><Id>$1<\/Id>/g;
        }
        if ( $dbase eq "gtr" ) {
          $data =~ s/<DocumentSummary uid=\"\d+\">/<DocumentSummary>/g;
        }
      }

      Encode::_utf8_on($data);

      if (! $raw) {
        $data = fix_bad_encoding($dbase, $data);
      }

      # remove eSummaryResult wrapper
      $data =~ s/<!DOCTYPE eSummaryResult PUBLIC/<!DOCTYPE DocumentSummarySet PUBLIC/g;
      $data =~ s/<eSummaryResult>//g;
      $data =~ s/<\/eSummaryResult>//g;

      if ( $json ) {
        $data = xml_to_json($data);
      }

      print "$data";
    }
  }
}

# eftch can read all arguments from the command line or participate in an EUtils pipe

my $ftch_help = qq{
Format Selection

  -format        Format of record or report
  -mode          text, xml, asn.1, json
  -style         withparts, conwithfeat

Direct Record Selection

  -db            Database name
  -id            Unique identifier or accession number

Sequence Range

  -seq_start     First sequence position to retrieve
  -seq_stop      Last sequence position to retrieve
  -strand        1 = forward DNA strand, 2 = reverse complement
  -revcomp       Shortcut for strand 2

Gene Range

  -chr_start     Sequence range from 0-based coordinates
  -chr_stop        in gene docsum GenomicInfoType object

Sequence Flags

  -complexity    0 = default, 1 = bioseq, 3 = nuc-prot set
  -extend        Extend sequence retrieval in both directions
  -extrafeat     Bit flag specifying extra features

Subset Retrieval

  -start         First record to fetch
  -stop          Last record to fetch

Miscellaneous

  -raw           Skip database-specific XML modifications
  -json          Convert adjusted XML output to JSON

Format Examples

  -db            -format            -mode    Report Type
  ___            _______            _____    ___________

  (all)
                 docsum                      DocumentSummarySet XML
                 docsum             json     DocumentSummarySet JSON
                 full                        Same as native except for mesh
                 uid                         Unique Identifier List
                 url                         Entrez URL
                 xml                         Same as -format full -mode xml

  bioproject
                 native                      BioProject Report
                 native             xml      RecordSet XML

  biosample
                 native                      BioSample Report
                 native             xml      BioSampleSet XML

  biosystems
                 native             xml      Sys-set XML

  gds
                 native             xml      RecordSet XML
                 summary                     Summary

  gene
                 full_report                 Detailed Report
                 gene_table                  Gene Table
                 native                      Gene Report
                 native             asn.1    Entrezgene ASN.1
                 native             xml      Entrezgene-Set XML
                 tabular                     Tabular Report

  homologene
                 alignmentscores             Alignment Scores
                 fasta                       FASTA
                 homologene                  Homologene Report
                 native                      Homologene List
                 native             asn.1    HG-Entry ASN.1
                 native             xml      Entrez-Homologene-Set XML

  mesh
                 full                        Full Record
                 native                      MeSH Report
                 native             xml      RecordSet XML

  nlmcatalog
                 native                      Full Record
                 native             xml      NLMCatalogRecordSet XML

  pmc
                 bioc                        PubTator Central BioC XML
                 medline                     MEDLINE
                 native             xml      pmc-articleset XML

  pubmed
                 abstract                    Abstract
                 bioc                        PubTator Central BioC XML
                 medline                     MEDLINE
                 native             asn.1    Pubmed-entry ASN.1
                 native             xml      PubmedArticleSet XML

  (sequences)
                 acc                         Accession Number
                 est                         EST Report
                 fasta                       FASTA
                 fasta              xml      TinySeq XML
                 fasta_cds_aa                FASTA of CDS Products
                 fasta_cds_na                FASTA of Coding Regions
                 ft                          Feature Table
                 gb                          GenBank Flatfile
                 gb                 xml      GBSet XML
                 gbc                xml      INSDSet XML
                 gene_fasta                  FASTA of Gene
                 gp                          GenPept Flatfile
                 gp                 xml      GBSet XML
                 gpc                xml      INSDSet XML
                 gss                         GSS Report
                 ipg                         Identical Protein Report
                 ipg                xml      IPGReportSet XML
                 native             text     Seq-entry ASN.1
                 native             xml      Bioseq-set XML
                 seqid                       Seq-id ASN.1

  snp
                 json                        Reference SNP Report

  sra
                 native             xml      EXPERIMENT_PACKAGE_SET XML
                 runinfo            xml      SraRunInfo XML

  structure
                 mmdb                        Ncbi-mime-asn1 strucseq ASN.1
                 native                      MMDB Report
                 native             xml      RecordSet XML

  taxonomy
                 native                      Taxonomy List
                 native             xml      TaxaSet XML

};

sub fix_sra_xml_encoding {

  my $data = shift (@_);

  $data =~ s/<!--[^<]+</</g;
  $data =~ s/>\s*-->/>/g;

  return $data;
}

sub fix_pubmed_xml_encoding {

  my $x = shift (@_);

  # my $markup = '(?:[biu]|su[bp])';
  # my $attrs = ' ?';
  my $markup = '(?:[\w.:_-]*:)?[[:lower:]-]+';
  my $attrs = '(?:\s[\w.:_-]+=[^>]*)?';

  # check for possible newline artifact
  $x =~ s|</$markup>\n||g;
  $x =~ s|\n<$markup$attrs>||g;

  # removed mixed content tags
  $x =~ s|</$markup>||g;
  $x =~ s|<$markup$attrs/?>||g;
  $x =~ s|</?DispFormula$attrs>| |g;

  # check for encoded tags
  if ( $x =~ /\&amp\;/ || $x =~ /\&lt\;/ || $x =~ /\&gt\;/ ) {
    # remove runs of amp
    $x =~ s|&amp;(?:amp;)+|&amp;|g;
    # fix secondary encoding
    $x =~ s|&amp;lt;|&lt;|g;
    $x =~ s|&amp;gt;|&gt;|g;
    $x =~ s|&amp;#(\d+);|&#$1;|g;
    # temporarily protect encoded scientific symbols, e.g., PMID 9698410 and 21892341
    $x =~ s|(?<= )(&lt;)(=*$markup&gt;)(?= )|$1=$2|g;
    # remove encoded markup
    $x =~ s|&lt;/$markup&gt;||g;
    $x =~ s|&lt;$markup$attrs/?&gt;||g;
    # undo temporary protection of scientific symbols adjacent to space
    $x =~ s|(?<= )(&lt;)=(=*$markup&gt;)(?= )|$1$2|g;
  }

  # compress runs of horizontal whitespace
  $x =~ s/\h+/ /g;

  # remove lines with just space
  $x =~ s/\n \n/\n/g;

  # remove spaces just outside of angle brackets
  $x =~ s|> |>|g;
  $x =~ s| <|<|g;

  # remove spaces just inside of parentheses
  $x =~ s|\( |\(|g;
  $x =~ s| \)|\)|g;

  # remove newlines flanking spaces
  $x =~ s|\n ||g;
  $x =~ s| \n| |g;

  return $x;
}

sub convert_bools {
    my %unrecognized;

    local *_convert_bools = sub {
        my $ref_type = ref($_[0]);
        if (!$ref_type) {
            # Nothing.
        }
        elsif ($ref_type eq 'HASH') {
            _convert_bools($_) for values(%{ $_[0] });
        }
        elsif ($ref_type eq 'ARRAY') {
            _convert_bools($_) for @{ $_[0] };
        }
        elsif (
               $ref_type eq 'JSON::PP::Boolean'           # JSON::PP
            || $ref_type eq 'Types::Serialiser::Boolean'  # JSON::XS
        ) {
            $_[0] = $_[0] ? 1 : 0;
        }
        else {
            ++$unrecognized{$ref_type};
        }
    };

    &_convert_bools;
}

sub xml_to_json {

  my $data = shift (@_);

  my $xc = new XML::Simple(KeepRoot => 1);
  my $conv = $xc->XMLin($data);
  convert_bools($conv);
  my $jc = JSON::PP->new->ascii->pretty->allow_nonref;
  $data = $jc->encode($conv);

  return $data;
}

sub bioc_do_get {

  my $urlx = shift (@_);
  my $argx = shift (@_);

  if ( $argx ne "" ) {
    $urlx .= "?";
    $urlx .= "$argx";
  }

  $rslt = get ($urlx);

    if ( ! defined $rslt ) {
      print STDERR "Failure of bioc_do_get '$urlx'\n";
      $result = 1;
      return "";
    }

    if ( $rslt eq "" ) {
      print STDERR "No bioc_do_get output returned from '$urlx'\n";
      $result = 1;
      return "";
    }


  if ( $rslt eq "" && $debug ) {
    print STDERR "No bioc_do_get output returned from '$urlx'\n";
  }

  if ( $debug ) {
    print STDERR "$rslt\n";
  }

  return $rslt;
}

sub do_icite_link {

  my $urlx = shift (@_);
  my $argx = shift (@_);

  if ( $argx ne "" ) {
    $urlx .= "?";
    $urlx .= "$argx";
  }

  $item = get ($urlx);

  if ( ! defined $item ) {
    die "Failure of get https://icite.od.nih.gov/api/pubs/\n";
  }

  if ( $item eq "" ) {
    die "No get output returned from get https://icite.od.nih.gov/api/pubs/\n";
  }

  # convert json to xml
  my $jc = JSON::PP->new->ascii->pretty->allow_nonref;
  my $conv = $jc->decode($item);
  convert_json($conv);
  my $result = XMLout($conv, SuppressEmpty => undef);

  my @unsorted = ();
  my @ids = ();
  if ( $cited ) {
    @ids = ($result =~ /<cited_by>(\d+)<\/cited_by>/g);
  } elsif ( $cites ) {
    @ids = ($result =~ /<references>(\d+)<\/references>/g);
  }

  foreach $uid (@ids) {
    push (@unsorted, $uid);
  }

  @sorted = sort { $a <=> $b } @unsorted;

  return @sorted;
}

# for id in 9698410 16271163 17282049 20968289 21892341 22785267 25435818 27672066 28635620 28976125 29547395
# do
#   efetch -db pubmed -format xml -id "$id" |
#   xtract -pattern PubmedArticle -plg "\n\n" -sep "\n\n" -tab "\n\n" \
#     -element MedlineCitation/PMID ArticleTitle AbstractText
# done

sub eftch {

  # ... | edirect.pl -fetch -format gp | ...

  clearflags ();

  MyGetOptions(
    $ftch_help,
    "db=s" => \$db,
    "id=s" => \$id,
    "format=s" => \$type,
    "docsum" => \$dcsm,
    "style=s" => \$style,
    "mode=s" => \$mode,
    "json" => \$json,
    "seq_start=i" => \$seq_start,
    "seq_stop=i" => \$seq_stop,
    "strand=s" => \$strand,
    "revcomp" => \$revcomp,
    "complexity=i" => \$complexity,
    "chr_start=i" => \$chr_start,
    "chr_stop=i" => \$chr_stop,
    "showgi" => \$showgi,
    "extend=i" => \$extend,
    "extrafeat=i" => \$extrafeat,
    "start=i" => \$min,
    "stop=i" => \$max,
    "api_key=s" => \$api_key,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "pipe" => \$pipe,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "internal" => \$internal,
    "log" => \$log,
    "compact" => \$compact,
    "raw" => \$raw,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "efetch $version\n";
    print $ftch_help;
    return;
  }

  # convert spaces between UIDs to commas

  $id =~ s/ /,/g;
  $id =~ s/,,/,/g;

  # efetch -docsum typo now repaired as efetch -format docsum

  if ( $type eq "" && $dcsm ) {
    $type = "docsum";
  }

  # "-format xml" is a shortcut for "-format full -mode xml"

  if ( $type eq "xml" and $mode eq "" ) {
    $type = "full";
    $mode = "xml";
  }

  if ( $style eq "normal" or $style eq "none" ) {
    $style = "";
  }

  if ( $style eq "conwithfeats" or $style eq "gbconwithfeat" or $style eq "gbconwithfeats" ) {
    $style = "conwithfeat";
  } elsif ( $style eq "withpart" or $style eq "gbwithpart" or $style eq "gbwithparts" ) {
    $style = "withparts";
  }

  if ( $type eq "gbconwithfeat" or $type eq "gbconwithfeats" ) {
    $type = "gb";
    $style = "conwithfeat";
  } elsif ( $type eq "gbwithparts" or $type eq "gbwithpart" ) {
    $type = "gb";
    $style = "withparts";
  }

  if ( $type eq "gbc" and $mode eq "" ) {
    $mode = "xml";
  }

  if ( -t STDIN and not @ARGV ) {
  } elsif ( $db ne "" and $id ne "" ) {
  } else {
    ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();
  }

  read_aliases ();
  adjust_base ();

  if ( $err ne "" ) {
    if ( ! $silent ) {
      die "ERROR in fetch input: $err\n\n";
    }
    return;
  }

  if ( $dbase eq "" ) {
    $dbase = $db;
  }

  $dbase = lc($dbase);

  if ( $type eq "" and $db ne "" ) {
    if ( get_zero_uid ($db) eq "" ) {
      die "Must supply -format report type on command line\n";
    }
    $type = "native";
  }

  if ( $mode eq "" ) {
    $mode = "text";
  }

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }

  if ( $revcomp ) {
    $strand = "2";
  }
  if ( $strand eq "plus" or $strand eq "+" ) {
    $strand = "1";
  }
  if ( $strand eq "minus" or $strand eq "-" or $strand eq "revcomp" ) {
    $strand = "2";
  }

  binmode STDOUT, ":utf8";

  # arguments can override loop start and stop

  if ( $min > 0 ) {
    $min--;
  }
  if ( $max == 0 ) {
    $max = $num
  }

  if ( $type eq "docsum" or $fnc eq "-summary" ) {

    esmry ( $dbase, $web, $key, $num, $id, $mode, $min, $max, $tool, $email,
            $silent, $verbose, $debug, $internal, $log, $http, $alias, $basx );

    return;
  }

  if ( $dbase eq "structure" and $type eq "mmdb" ) {

    emmdb ( $dbase, $web, $key, $num, $id, $tool, $email );

    return;
  }

  if ( $dbase ne "" and ( $type eq "UID" or $type eq "uid" ) ) {

    if ( $id ne "" ) {

      my @ids = split (',', $id);
      foreach $uid (@ids) {
        $uid =~ s/(\d+)\.\d+/$1/g;
        print "$uid\n";
      }

      return;
    }

    if ( $web eq "" ) {
      die "WebEnv value not found in fetch input\n";
    }

    if ( $pipe ) {
      write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );
    }

    if ( $max == 0 ) {
      if ( $silent ) {
        return;
      }
    }

    # use larger chunk for UID format
    $chunk = 25000;
    for ( $start = $min; $start < $max; $start += $chunk ) {

      my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $max, $tool, $email );
      foreach $uid (@ids) {
        print "$uid\n";
      }
    }

    return;
  }

  if ( $dbase ne "" and ( $type eq "URL" or $type eq "url" ) ) {

    if ( $id ne "" ) {

      my @ids = split (',', $id);
      $url = "https://www.ncbi.nlm.nih.gov/";
      $url .= "$dbase/";
      my $pfx = "";
      foreach $uid (@ids) {
        $uid =~ s/(\d+)\.\d+/$1/g;
        $url .= "$pfx$uid";
        $pfx = ",";
      }
      print "$url\n";

      return;
    }

    if ( $web eq "" ) {
      die "WebEnv value not found in fetch input\n";
    }

    if ( $pipe ) {
      write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );
    }

    if ( $max == 0 ) {
      if ( $silent ) {
        return;
      }
    }

    # use larger chunk for URL format
    $chunk = 25000;
    for ( $start = $min; $start < $max; $start += $chunk ) {

      my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $max, $tool, $email );
      $url = "https://www.ncbi.nlm.nih.gov/";
      $url .= "$dbase/";
      my $pfx = "";
      foreach $uid (@ids) {
        $url .= "$pfx$uid";
        $pfx = ",";
      }
      print "$url\n";
    }

    return;
  }

  if ( $dbase ne "" and ( $type eq "URLS" or $type eq "urls" ) ) {

    if ( $id ne "" ) {

      my @ids = split (',', $id);
      foreach $uid (@ids) {
        $uid =~ s/(\d+)\.\d+/$1/g;
        print "https://www.ncbi.nlm.nih.gov/$dbase/$uid\n";
      }

      return;
    }

    if ( $web eq "" ) {
      die "WebEnv value not found in fetch input\n";
    }

    if ( $pipe ) {
      write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );
    }

    if ( $max == 0 ) {
      if ( $silent ) {
        return;
      }
    }

    # use larger chunk for URL format
    $chunk = 25000;
    for ( $start = $min; $start < $max; $start += $chunk ) {

      my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $max, $tool, $email );
      foreach $uid (@ids) {
        print "https://www.ncbi.nlm.nih.gov/$dbase/$uid\n";
      }
    }

    return;
  }

  if ( $dbase eq "pubmed" and $type eq "bioc" and $id ne "" and $id ne "0" ) {

    # efetch -db pubmed -id 28483577 -format bioc

    $url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocxml";
    $arg = "pmids=$id";

    $output = bioc_do_get ($url, $arg);

    if ( $output !~ /\n$/ ) {
      $output .= "\n";
    }

    print $output;

    return;
  }

  if ( $dbase eq "pmc" and $type eq "bioc" and $id ne "" and $id ne "0" ) {

    # efetch -db pmc -id 6207735,3681088 -format bioc

    my @accum = ();
    my @ids = split (',', $id);
    foreach $pmcid (@ids) {
      if ( $pmcid !~ /^PMC/ ) {
        # add PMC prefix if not already in argument
        $pmcid = "PMC" . $pmcid;
      }
      push (@accum, $pmcid);
    }
    $id = join (',', @accum);

    $url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocxml";
    $arg = "pmcids=$id";

    $output = bioc_do_get ($url, $arg);

    if ( $output !~ /\n$/ ) {
      $output .= "\n";
    }

    print $output;

    return;
  }

  if ( $dbase ne "" and $id ne "" ) {

    if ( $id eq "0" ) {

      # id "0" returns a live UID for any database

      $id = get_zero_uid ($dbase);

      if ( $id eq "0" ) {

        # id "0" is an unrecognized accession

        return;
      }
    }

    $url = $base . $efetch;

    if ( $dbase eq "pubmed" ) {
      $id =~ s/(\d+)\.\d+/$1/g;
    }

    $arg = "db=$dbase&id=$id";

    if ( $type eq "gb" or $type eq "gbc" ) {
      if ( $style eq "withparts" or $style eq "master" ) {
        $arg .= "&rettype=$type";
        $arg .= "&retmode=$mode";
        $arg .= "&style=$style";
      } elsif ( $style eq "conwithfeat" or $style eq "withfeat" or $style eq "contigwithfeat" ) {
        $arg .= "&rettype=$type";
        $arg .= "&retmode=$mode";
        $arg .= "&gbconwithfeat=1";
      } else {
        $arg .= "&rettype=$type";
        $arg .= "&retmode=$mode";
      }
    } else {
      $arg .= "&rettype=$type";
      $arg .= "&retmode=$mode";
    }

    # -chr_start and -chr_stop are for 0-based sequence coordinates from EntrezGene
    if ( $chr_start > -1 && $chr_stop > -1 ) {
      $seq_start = $chr_start + 1;
      $seq_stop = $chr_stop + 1;
    }

    # if -seq_start > -seq_stop, swap values to normalize, indicate minus strand with -strand 2
    if ( $seq_start > 0 && $seq_stop > 0 ) {
      if ( $seq_start > $seq_stop ) {
        my $tmp = $seq_start;
        $seq_start = $seq_stop;
        $seq_stop = $tmp;
        $strand = "2";
      }
    }

    # option to show GI number (undocumented)
    if ( $showgi ) {
      $arg .= "&showgi=1";
    }

    # optionally extend retrieved sequence range in both directions
    if ( $extend > 0 ) {
      $seq_start -= $extend;
      $seq_stop += $extend;
    }

    if ( $strand ne "" ) {
      $arg .= "&strand=$strand";
    }
    if ( $seq_start > 0 ) {
      $arg .= "&seq_start=$seq_start";
    }
    if ( $seq_stop > 0 ) {
      $arg .= "&seq_stop=$seq_stop";
    }
    if ( $complexity > 0 ) {
      $arg .= "&complexity=$complexity";
    }
    if ( $extrafeat > -1 ) {
      $arg .= "&extrafeat=$extrafeat";
    }

    $data = do_post_yielding_ref ($url, $arg, $tool, $email, true);

    if ($$data =~ /<eFetchResult>/i and $$data =~ /<ERROR>(.+?)<\/ERROR>/i) {
      $err = $1;
      if ( $err ne "" ) {
        if ( ! $silent ) {
          print STDERR "ERROR in efetch: $err\n";
        }
      }
    }

    Encode::_utf8_on($$data);

    if (! $raw) {

      if ( $dbase eq "sra" and $type eq "full" and $mode eq "xml" ) {
        $$data = fix_sra_xml_encoding($$data);
      }

      if ( $dbase eq "pubmed" and $type eq "full" and $mode eq "xml" ) {
        $$data = fix_pubmed_xml_encoding($$data);
      }

      if ( $type eq "fasta" or $type eq "fasta_cds_aa" or $type eq "fasta_cds_na" or $type eq "gene_fasta" ) {
        # remove blank lines in FASTA format
        $$data =~ s/\n+/\n/g;
      }
    }

    if ( $$data !~ /\n$/ ) {
      $$data .= "\n";
    }

    if ( $json ) {
      $$data = xml_to_json($$data);
    }

    print $$data;

    return;
  }

  if ( $dbase ne "" and $web ne "" and $key eq "" and $num eq "0" ) {
    die "QueryKey value not found in fetch input\n";
    return;
  }

  test_edirect ( $dbase, $web, $key, $num, "fetch" );

  if ( $stp ne "" ) {
    $stpminusone = $stp - 1;
  }

  if ( $max == 0 ) {
    if ( $silent ) {
      return;
    }
  }

  if ( $dbase eq "pubmed" and $type eq "bioc" and $id eq "" ) {

    # echo "28483577" | epost -db pubmed -format uid | efetch -format bioc

    # use chunk of 100 for PubTator BioC GET
    $chunk = 100;

    for ( $start = $min; $start < $max; $start += $chunk ) {

      my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $max, $tool, $email );
      $id = join (',', @ids);

      $url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocxml";
      $arg = "pmids=$id";

      $output = bioc_do_get ($url, $arg);

      if ( $output !~ /\n$/ ) {
        $output .= "\n";
      }

      print $output;
    }

    return;
  }

  if ( $dbase eq "pmc" and $type eq "bioc" and $id eq "" ) {

    # echo -e "6207735\n3681088" | epost -db pmc -format uid | efetch -format bioc

    # use chunk of 100 for PubTator BioC GET
    $chunk = 100;

    for ( $start = $min; $start < $max; $start += $chunk ) {

      my @accum = ();
      my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $max, $tool, $email );
      foreach $pmcid (@ids) {
        if ( $pmcid !~ /^PMC/ ) {
          # add PMC prefix if not already in argument
          $pmcid = "PMC" . $pmcid;
        }
        push (@accum, $pmcid);
      }
      $id = join (',', @accum);

      $url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocxml";
      $arg = "pmcids=$id";

      $output = bioc_do_get ($url, $arg);

      if ( $output !~ /\n$/ ) {
        $output .= "\n";
      }

      print $output;
    }

    return;
  }

  # use small chunk because fetched records could be quite large
  $chunk = 500;

  # use larger chunk for accessions
  if ( $dbase eq "nucleotide" or
       $dbase eq "nuccore" or
       $dbase eq "est" or
       $dbase eq "gss" or
       $dbase eq "protein" ) {
    if ( $type eq "ACCN" or $type eq "accn" or $type eq "ACC" or $type eq "acc" ) {
      $chunk = 10000;
    }
  }

  # use smaller chunk for SNP JSON
  if ( $dbase eq "snp" ) {
    if ( $type eq "JSON" or $type eq "json" ) {
      # $chunk = 100;
      $chunk = 10;
    }
  }

  for ( $start = $min; $start < $max; $start += $chunk ) {
    $url = $base . $efetch;

    $chkx = $chunk;
    if ( $start + $chkx > $max ) {
      $chkx = $max - $start;
    }

    $arg = "db=$dbase&query_key=$key&WebEnv=$web";

    if ( $type eq "gb" or $type eq "gbc" ) {
      if ( $style eq "withparts" or $style eq "master" ) {
        $arg .= "&rettype=$type";
        $arg .= "&retmode=$mode";
        $arg .= "&style=$style";
      } elsif ( $style eq "conwithfeat" or $style eq "withfeat" or $style eq "contigwithfeat" ) {
        $arg .= "&rettype=$type";
        $arg .= "&retmode=$mode";
        $arg .= "&gbconwithfeat=1";
      } else {
        $arg .= "&rettype=$type";
        $arg .= "&retmode=$mode";
      }
    } else {
      $arg .= "&rettype=$type";
      $arg .= "&retmode=$mode";
    }

    $arg .= "&retstart=$start&retmax=$chkx";
    if ( $strand ne "" ) {
      $arg .= "&strand=$strand";
    }
    if ( $seq_start > 0 ) {
      $arg .= "&seq_start=$seq_start";
    }
    if ( $seq_stop > 0 ) {
      $arg .= "&seq_stop=$seq_stop";
    }
    if ( $complexity > 0 ) {
      $arg .= "&complexity=$complexity";
    }
    if ( $extrafeat > -1 ) {
      $arg .= "&extrafeat=$extrafeat";
    }

    $data = "";
    $retry = true;

    for ( $tries = 0; $tries < 3 && $retry; $tries++) {
      $data = do_post_yielding_ref ($url, $arg, $tool, $email, true);

      if ($$data =~ /<eFetchResult>/i and $$data =~ /<ERROR>(.+?)<\/ERROR>/i) {
        if ( ! $silent ) {
          print STDERR "Retrying efetch, step $stp: $err\n";
        }
        sleep 3;
      } else {
        $retry = false;
      }
    }

    if ($$data =~ /<eFetchResult>/i and $$data =~ /<ERROR>(.+?)<\/ERROR>/i) {
      $err = $1;
      if ( $err ne "" ) {
        if ( ! $silent ) {
          my $from = $start + 1;
          my $to = $start + $chunk;
          if ( $to > $num ) {
            $to = $num;
          }
          print STDERR "ERROR in efetch ($from-$to / $num): $err\n";
          print STDERR "Replicate for debugging with:\n";
          print STDERR "  edirutil -db $dbase -web $web -key $key -count $num";
          if ( $stpminusone > 0 ) {
            print STDERR " -step $stpminusone";
          }
          my $seconds = "";
          my $end_time = Time::HiRes::time();
          my $elapsed = $end_time - $begin_time;
          if ( $elapsed > 0.0005 ) {
            if ( $elapsed =~ /(\d+)\.(\d\d\d)\d+/ ) {
              $seconds = "$1.$2";
            }
          }
          if ( $seconds ne "" ) {
            print STDERR " -seconds $seconds";
          }
          print STDERR " | efetch -format $type";
          if ( $mode ne "text" ) {
            print STDERR " -mode $mode";
          }
          print STDERR " -start $from -stop $to\n";
        }
      }
    } else {
      if ( $verbose ) {
        my $from = $start + 1;
        my $to = $start + $chunk;
        if ( $to > $num ) {
          $to = $num;
        }
        print STDERR "( edirutil -db $dbase -web $web -key $key -count $num";
        if ( $stpminusone > 0 ) {
          print STDERR " -step $stpminusone";
        }
        my $seconds = "";
        my $end_time = Time::HiRes::time();
        my $elapsed = $end_time - $begin_time;
        if ( $elapsed > 0.0005 ) {
          if ( $elapsed =~ /(\d+)\.(\d\d\d)\d+/ ) {
            $seconds = "$1.$2";
          }
        }
        if ( $seconds ne "" ) {
          print STDERR " -seconds $seconds";
        }
        print STDERR " | efetch -format $type";
        if ( $mode ne "text" ) {
          print STDERR " -mode $mode";
        }
        print STDERR " -start $from -stop $to )\n";
        $begin_time = $end_time;
      }

      Encode::_utf8_on($$data);

      if (! $raw) {
        if ( $dbase eq "sra" and $type eq "full" and $mode eq "xml" ) {
          $$data = fix_sra_xml_encoding($$data);
        }

        if ( $dbase eq "pubmed" and $type eq "full" and $mode eq "xml" ) {
          $$data = fix_pubmed_xml_encoding($$data);
        }

        if ( $type eq "fasta" or $type eq "fasta_cds_aa" or $type eq "fasta_cds_na" or $type eq "gene_fasta" ) {
          # remove blank lines in FASTA format
          $$data =~ s/\n+/\n/g;
        }
      }

      if ( $$data !~ /\n$/ ) {
        $$data .= "\n";
      }

      if ( $json ) {
        $$data = xml_to_json($$data);
      }

      print $$data;
    }
  }
}

# einfo obtains names of databases or names of fields and links per database

my $info_help = qq{
Database Selection

  -db        Database name
  -dbs       Get all database names

Data Summaries

  -fields    Print field names
  -links     Print link names

Field Example

  <Field>
    <Name>ALL</Name>
    <FullName>All Fields</FullName>
    <Description>All terms from all searchable fields</Description>
    <TermCount>138982028</TermCount>
    <IsDate>N</IsDate>
    <IsNumerical>N</IsNumerical>
    <SingleToken>N</SingleToken>
    <Hierarchy>N</Hierarchy>
    <IsHidden>N</IsHidden>
    <IsTruncatable>Y</IsTruncatable>
    <IsRangable>N</IsRangable>
  </Field>

Link Example

  <Link>
    <Name>pubmed_protein</Name>
    <Menu>Protein Links</Menu>
    <Description>Published protein sequences</Description>
    <DbTo>protein</DbTo>
  </Link>
  <Link>
    <Name>pubmed_protein_refseq</Name>
    <Menu>Protein (RefSeq) Links</Menu>
    <Description>Link to Protein RefSeqs</Description>
    <DbTo>protein</DbTo>
  </Link>

};

sub einfo {

  # ... | edirect.pl -info -db pubmed | ...

  clearflags ();

  MyGetOptions(
    $info_help,
    "db=s" => \$db,
    "dbs" => \$dbs,
    "fields" => \$fields,
    "links" => \$links,
    "api_key=s" => \$api_key,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "internal" => \$internal,
    "log" => \$log,
    "compact" => \$compact,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "einfo $version\n";
    print $info_help;
    return;
  }

  read_aliases ();
  adjust_base ();

  if ( @ARGV  ||  ($^O =~ /^MSWin/ && ! -t STDIN) ) {
    while ( defined($thisline = <STDIN>) ) {
      $tool = $1 if ( $thisline =~ /<Tool>(.+?)<\/Tool>/i );
      $email = $1 if ( $thisline =~ /<Email>(.+?)<\/Email>/i );
    }
  }

  if ( $dbase eq "" ) {
    $dbase = $db;
  }

  $dbase = lc($dbase);

  if ( $dbase eq "" and (! $dbs) ) {
    die "Must supply -db or -dbs on command line\n";
  }

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }

  binmode STDOUT, ":utf8";

  $url = $base . $einfo;

  $prefix = "?";

  if ( $dbase ne "" ) {
    $url .= "$prefix" . "db=$dbase&version=2.0";
    $prefix = "&";
  }

  if ( $os ne "" ) {
    $url .= "$prefix" . "edirect_os=$os";
    $prefix = "&";
  }

  if ( $api_key ne "" ) {
    $url .= "$prefix" . "api_key=$api_key";
    $prefix = "&";
  }

  $url .= "$prefix" . "edirect=$version";
  $prefix = "&";

  if ( $tool eq "" ) {
    $tool = "edirect";
  }
  if ( $tool ne "" ) {
    $url .= "$prefix" . "tool=$tool";
    $prefix = "&";
  }

  if ( $email eq "" ) {
    $email = get_email ();
  }
  if ( $email ne "" ) {
    $url .= "$prefix" . "email=$email";
    $prefix = "&";
  }

  if ( $debug or $log ) {
    print STDERR "$url\n";
  }

  do_sleep ();

  $output = get ($url);

  if ( ! defined $output ) {
    print STDERR "Failure of '$url'\n";
    $result = 1;
    return;
  }

  if ( $output eq "" ) {
    print STDERR "No einfo output returned from '$url'\n";
    return;
  }

  if ($output =~ /IsTrunc/ or $output =~ /IsRang/) {
    # print STDERR "New server is out with IsTrunc/IsRang - update executable\n";
  }

  if ( $debug ) {
    print STDERR "$output\n";
  }

  if ( $dbs and $output =~ /<DbName>/ ) {

    # -dbs now extracts database names from XML

    $output =~ s/\r//g;
    $output =~ s/\n//g;
    $output =~ s/\t//g;
    $output =~ s/ +/ /g;
    $output =~ s/> +</></g;

    my @databases = ($output =~ /<DbName>(.+?)<\/DbName>/g);
    foreach $dtbs (@databases) {
      print "$dtbs\n";
    }

    return;
  }

  if ( ( $fields or $links ) and $output =~ /<DbInfo>/ ) {

    # -db can print information directly without need to process XML with xtract

    $output =~ s/\r//g;
    $output =~ s/\n//g;
    $output =~ s/\t//g;
    $output =~ s/ +/ /g;
    $output =~ s/> +</></g;

    my $name = "";
    my $full = "";
    my $menu = "";

    if ( $fields ) {
      my @unsorted = ();
      my @flds = ($output =~ /<Field>(.+?)<\/Field>/g);
      foreach $fld (@flds) {
        $name = "";
        $full = "_";
        if ( $fld =~ /<Name>(.+?)<\/Name>/ ) {
          $name = $1;
        }
        if ( $fld =~ /<FullName>(.+?)<\/FullName>/ ) {
          $full = $1;
        }
        if ( $name ne "" and $full ne "" ) {
          push (@unsorted, "$name\t$full\n");
        }
      }
      my @sorted = sort { "\U$a" cmp "\U$b" } @unsorted;
      foreach $itm (@sorted) {
        print "$itm";
      }
    }

    if ( $links ) {
      my @unsorted = ();
      my @lnks = ($output =~ /<Link>(.+?)<\/Link>/g);
      foreach $lnk (@lnks) {
        $name = "";
        $menu = "_";
        if ( $lnk =~ /<Name>(.+?)<\/Name>/ ) {
          $name = $1;
        }
        if ( $lnk =~ /<Menu>(.+?)<\/Menu>/ ) {
          $menu = $1;
        }
        if ( $name ne "" and $menu ne "" ) {
          push (@unsorted, "$name\t$menu\n");
        }
      }
      my @sorted = sort { "\U$a" cmp "\U$b" } @unsorted;
      foreach $itm (@sorted) {
        print "$itm";
      }
    }

    return;
  }

  print "$output";
}

# common link history processing

sub acheck_test {

  my $dbase = shift (@_);
  my $dbto = shift (@_);
  my $name = shift (@_);
  my $web = shift (@_);
  my $key = shift (@_);
  my $num = shift (@_);
  my $stp = shift (@_);
  my $err = shift (@_);
  my $tool = shift (@_);
  my $email = shift (@_);

  if ( $num > 0 ) {
    return;
  }

  # if failure, use acheck to confirm that there are no links

  $url = $base . $elink;
  $arg = "dbfrom=$dbase&query_key=$key&WebEnv=$web&cmd=acheck";
  $data = do_post ($url, $arg, $tool, $email, true);

  # remove newlines, tabs, space between tokens, compress runs of spaces,

  $data =~ s/\r//g;
  $data =~ s/\n//g;
  $data =~ s/\t//g;
  $data =~ s/ +/ /g;
  $data =~ s/> </></g;

  if ($data =~ /<Error>(.+?)<\/Error>/i) {
    $err = $1;
    if ( $err ne "" ) {
      write_edirect ( $dbto, $web, $key, $num, $stp, $err, $tool, $email );
      close (STDOUT);
      die "ERROR in acheck test: $err\n";
    }
  }

  if ( $data !~ /<LinkInfo>/ ) {
    write_edirect ( $dbto, $web, $key, $num, $stp, $err, $tool, $email );
    close (STDOUT);
    die "Elink acheck confirmation test failed\n";
  }

  if ( $data =~ /<LinkName>$name<\/LinkName>/ ) {
    write_edirect ( $dbto, $web, $key, $num, $stp, $err, $tool, $email );
    close (STDOUT);
    die "Elink acheck test indicates non-zero count expected\n";
  }
}

sub process_history_link {

  my $arg = shift (@_);
  my $output = shift (@_);
  my $dbase = shift (@_);
  my $dbto = shift (@_);
  my $name = shift (@_);
  my $wb = shift (@_);
  my $web = shift (@_);
  my $key = shift (@_);
  my $stp = shift (@_);
  my $err = shift (@_);
  my $lbl = shift (@_);
  my $tool = shift (@_);
  my $email = shift (@_);

  my $num = "";

  if ( $wb eq "" ) {
    $wb = $web;
  }

  if ( $err ne "" ) {
    write_edirect ( "", $wb, "", "", "", $err, "", "" );
    close (STDOUT);
    if ( ! $silent ) {
      die "ERROR in link output: $err\nWebEnv: $wb\nURL: $arg\nResult: $output\n\n";
    }
    return;
  }

  if ( $web eq "" ) {
    die "WebEnv value not found in link output - WebEnv1 $wb\n";
  }
  if ( $key eq "" ) {
    write_edirect ( $dbto, $wb, $key, "0", $stp, $err, $tool, $email );
    close (STDOUT);
    # no neighbors or links can be a normal response,
    # e.g., elink -db gene -id 496376 -target medgen
    # so suppress this message
    # die "QueryKey value not found in link output - WebEnv1 $wb\n";
    return;
  }

  if ( $web ne $wb ) {
    $err = "WebEnv mismatch in link-search output - WebEnv1 $wb, WebEnv2 $web";
    write_edirect ( "", $web, "", "", "", $err, "", "" );
    close (STDOUT);
    die "WebEnv value changed in link-search output - WebEnv1 $wb\nWebEnv2 $web\n";
  }

  ( $num, $key ) = get_count ( $dbto, $web, $key, $tool, $email );

  if ( $num eq "" ) {
    $err = "Missing count in link-search output - WebEnv $web";
    write_edirect ( "", $web, "", "", "", $err, "", "" );
    close (STDOUT);
    die "Count value not found in link-search output - WebEnv $web\n";
  }

  if ( $num == 0 ) {
    acheck_test ( $dbase, $dbto, $name, $web, $key, $num, $stp, $err, $tool, $email );
  }

  if ( $lbl ne "" and $key ne "" ) {
    $labels{"$lbl"} = "$key";
  }

  write_edirect ( $dbto, $web, $key, $num, $stp, $err, $tool, $email );
}

# for large lists, break into chunks, merge result on client, post to server

sub batch_elink {

  my $dbase = shift (@_);
  my $dbto = shift (@_);
  my $name = shift (@_);
  my $web = shift (@_);
  my $key = shift (@_);
  my $num = shift (@_);
  my $stp = shift (@_);
  my $err = shift (@_);
  my $lbl = shift (@_);
  my $tool = shift (@_);
  my $email = shift (@_);
  my $auto = shift (@_);

  my %seen = ();
  my @uniq = ();

  if ( $num == 0 ) {
    if ( $silent ) {
      write_edirect ( "", "", "", "", "", $err, "", "" );
      close (STDOUT);
      return;
    }
  }

  $chunk = 500;
  for ( $start = 0; $start < $num; $start += $chunk ) {

    do_sleep ();

    my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $num, $tool, $email );

    $url = $base . $elink;

    $arg = "dbfrom=$dbase&db=$dbto&cmd=neighbor&linkname=$name";
    if ( $mode ne "" ) {
      $arg .= "&retmode=$mode";
    }
    $arg .= "&id=";
    $arg .= join (',', @ids);

    $data = "";
    $retry = true;

    for ( $tries = 0; $tries < 3 && $retry; $tries++) {

      do_sleep ();

      $data = do_post ($url, $arg, $tool, $email, true);

      if ($data =~ /<Error>(.+?)<\/Error>/i) {
        if ( ! $silent ) {
          print STDERR "Retrying batch elink, step $stp: $err\n";
        }
        sleep 3;
      } else {
        $retry = false;
      }
    }

    # remove newlines, tabs, space between tokens, compress runs of spaces,

    $data =~ s/\r//g;
    $data =~ s/\n//g;
    $data =~ s/\t//g;
    $data =~ s/ +/ /g;
    $data =~ s/> </></g;

    if ($data =~ /<Error>(.+?)<\/Error>/i) {
      $err = $1;
      if ( $err ne "" ) {
        if ( ! $silent ) {
          my $from = $start + 1;
          my $to = $start + $chunk;
          if ( $to > $num ) {
            $to = $num;
          }
          print STDERR "ERROR in batch elink ($from-$to / $num): $err\n";
        }
      }
    }

    if ($data !~ /<LinkSetDb>/i) {
      if ( ! $silent ) {
        print STDERR "LinkSetDb missing in batch elink result\n";
      }
    }

    while ( $data =~ /<LinkSetDb>(.*?)<\/LinkSetDb>/g ) {
      $linkset = $1;
      while ( $linkset =~ /<Id>(\d+)<\/Id>/g ) {
        $uid = $1;
        if ( ! $seen{$uid}++ ) {
          push (@uniq, $uid);
        }
      }
    }
  }

  $dbase = $dbto;

  $dbase = lc($dbase);

  $num = scalar @uniq;

  if ( $num == 0 ) {

    do_sleep ();

    acheck_test ( $dbase, $dbto, $name, $web, $key, $num, $stp, $err, $tool, $email );
    write_edirect ( $dbase, $web, $key, "0", $stp, $err, $tool, $email );
    if ( $auto ) {
      close (STDOUT);
      die "Automatic batch link failed\n";
    }
    return;
  }

  # not certain if sort is necessary - need to experiment before final release

  @sorted = sort { $a <=> $b } @uniq;

  $url = $base . $epost;

  $arg = "db=$dbase";
  if ( $key ne "" ) {
    $arg .= "&query_key=$key";
  }
  if ( $web ne "" ) {
    $arg .= "&WebEnv=$web";
  }
  $ids = join (',', @sorted);

  $arg .= "&id=$ids";

  $output = do_post ($url, $arg, $tool, $email, true);

  if ( $debug ) {
    print STDERR "$output\n";
  }

  $wb = $web;

  $web = "";
  $key = "";
  $err = "";

  $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
  $key = $1 if ($output =~ /<QueryKey>(\S+)<\/QueryKey>/);
  $err = $1 if ($output =~ /<Error>(.+?)<\/Error>/i);

  process_history_link ( $arg, $output, $dbase, $dbto, $name, $wb, $web, $key, $stp, $err, $lbl, $tool, $email );
}

# elink without a target uses the source database for neighboring

my $link_help = qq{
Destination Database

  -related    Neighbors in same database
  -target     Links in different database
  -name       Link name (e.g., pubmed_protein_refseq, pubmed_pubmed_citedin, pubmed_pubmed_refs)

Direct Record Selection

  -db         Database name
  -id         Unique identifier(s)

Advanced Control

  -cmd        Command type (returns eLinkResult XML)
  -mode       "ref" uses LinkOut provider's web site
  -holding    Name of LinkOut provider

PubMed Citation Lookup

  -cited      References to this paper
  -cites      Publication reference list

Batch Processing

  -batch      Bypass Entrez history mechanism

Miscellaneous Arguments

  -label      Alias for query step

Command Option Examples

  -cmd              Result
  ____              ______

  neighbor          Neighbors or links

  neighbor_score    Neighbors with computed similarity scores

  acheck            All links available

  ncheck            Existence of neighbors

  lcheck            Existence of external links (LinkOuts)

  llinks            Non-library LinkOut providers

  llinkslib         All LinkOut providers

  prlinks           Primary LinkOut provider,
                    or URL for single UID with -mode ref

};

sub elink {

  # ... | edirect.pl -link -target structure | ...

  clearflags ();

  MyGetOptions(
    $link_help,
    "db=s" => \$db,
    "id=s" => \$id,
    "format=s" => \$type,
    "target=s" => \$dbto,
    "name=s" => \$name,
    "related" => \$related,
    "neighbor" => \$neighbor,
    "cmd=s" => \$cmd,
    "mode=s" => \$mode,
    "cited" => \$cited,
    "cites" => \$cites,
    "batch" => \$batch,
    "holding=s" => \$holding,
    "label=s" => \$lbl,
    "api_key=s" => \$api_key,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "internal" => \$internal,
    "log" => \$log,
    "compact" => \$compact,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "elink $version\n";
    print $link_help;
    return;
  }

  # convert spaces between UIDs to commas

  $id =~ s/ /,/g;
  $id =~ s/,,/,/g;

  if ( -t STDIN and not @ARGV ) {
  } elsif ( $db ne "" and $id ne "" ) {
  } else {
    ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();
  }

  read_aliases ();
  adjust_base ();

  if ( $err ne "" ) {
    write_edirect ( "", "", "", "", "", $err, "", "" );
    close (STDOUT);
    if ( ! $silent ) {
      die "ERROR in link input: $err\n\n";
    }
    return;
  }

  if ( $dbase eq "" ) {
    $dbase = $db;
  }

  $dbase = lc($dbase);

  my $adddbto = true;
  if ( $dbto eq "" ) {
    if ( $cmd eq "acheck" or
         $cmd eq "ncheck" or
         $cmd eq "lcheck" or
         $cmd eq "llinks" or
         $cmd eq "llinkslib" or
         $cmd eq "prlinks" ) {
      $dbto = $dbase;
      $adddbto = false;
    }
  }

  if ( $dbto eq "" and (! $related) and (! $neighbor) and (! $cited) and (! $cites) ) {
    die "Must supply -target or -related on command line\n";
  }

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }

  binmode STDOUT, ":utf8";

  if ( $dbto eq "" ) {
    $dbto = $dbase;
  }
  if ( $name eq "" ) {
    $name = $dbase . "_" . $dbto;
  }

  if ( $cmd eq "" ) {
    $cmd = "neighbor_history";
  }

  if ( $dbase eq "nlmcatalog" ) {
    die "Entrez Direct does not support links for the nlmcatalog database\n";
  }

  # experimental access to new icite service
  # equivalent of (PMC-requiring) -name pubmed_pubmed_citedin or -name pubmed_pubmed_refs

  if ( $dbase eq "pubmed" and ( $cited or $cites ) ) {

    my @unsorted = ();

    if ( $id ne "" and $id ne "0" ) {

      # elink -db pubmed -id 2539356 -cited
      # elink -db pubmed -id 2539356 -cites

      my @ids = split (',', $id);
      foreach $pmid (@ids) {
        push (@unsorted, $pmid);
      }

    } else {

      # esearch -db pubmed -query "NCBI [AFFL]" | elink -cited

      if ( $web eq "" ) {
        die "WebEnv value not found in link input\n";
      }

      $chunk = 5000;
      for ( $start = 0; $start < $num; $start += $chunk ) {

        my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $num, $tool, $email );
        foreach $pmid (@ids) {
          push (@unsorted, $pmid);
        }
      }
    }

    if ( scalar @unsorted < 1 ) {
      write_edirect ( $dbase, $web, $key, "0", $stp, $err, $tool, $email );
      return
    }

    @sorted = sort { $a <=> $b } @unsorted;

    my @uniqued = ();

    my $prev = "";
    foreach $pmid (@sorted) {
      if ( $pmid ne $prev ) {
        push (@uniqued, $pmid);
        $prev = $pmid;
      }
    }

    my @working = ();

    while ( my @chunk = splice @uniqued, 0, 100 ) {

      $id = join (',', @chunk);

      $url = "https://icite.od.nih.gov/api/pubs";
      $arg = "pmids=$id";

      my @idc = do_icite_link ($url, $arg);

      foreach $pmid (@idc) {
        push (@working, $pmid);
      }
    }

    if ( scalar @working < 1 ) {
      write_edirect ( $dbase, $web, $key, "0", $stp, $err, $tool, $email );
      return
    }

    @combined = sort { $a <=> $b } @working;

    my @resolved = ();

    my $last = "";
    foreach $pmid (@combined) {
      if ( $pmid ne $last ) {
        push (@resolved, $pmid);
        $last = $pmid;
      }
    }

    $idx = join (',', @resolved);

    ( $web, $key ) = post_chunk ( $dbase, $web, $key, $tool, $email, $idx, "" );

    ( $num, $key ) = get_count ( $dbase, $web, $key, $tool, $email );

    write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );

    return;
  }

  if ( $dbase ne "" and $id ne "" ) {

    # process db and id command-line arguments instead of getting from history

    $url = $base . $elink;

    $arg = "dbfrom=$dbase";
    if ( $adddbto ) {
      $arg .= "&db=$dbto";
    }
    $arg .= "&cmd=$cmd&linkname=$name";
    if ( $mode ne "" ) {
      $arg .= "&retmode=$mode";
    }
    if ( $dbase eq "pubmed" ) {
      $id =~ s/(\d+)\.\d+/$1/g;
    }
    $arg .= "&id=$id";
    if ( $type eq "acc" ) {
      $arg .= "&idtype=acc";
    }

    $data = do_post ($url, $arg, $tool, $email, true);

    if ( $cmd ne "neighbor_history" ) {

      # if not neighbor_history, write eLinkResult XML instead of ENTREZ_DIRECT

      print "$data";

      return;
    }

    $wb = $web;

    $web = "";
    $key = "";
    $err = "";

    $err = $1 if ($data =~ /<Error>(.+?)<\/Error>/i);
    if ( $err eq "" ) {
      $web = $1 if ($data =~ /<WebEnv>(\S+)<\/WebEnv>/);
      $key = $1 if ($data =~ /<QueryKey>(\S+)<\/QueryKey>/);
    }

    process_history_link ( $arg, $data, $dbase, $dbto, $name, $wb, $web, $key, $stp, $err, $lbl, $tool, $email );

    return;
  }

  test_edirect ( $dbase, $web, $key, $num, "link" );

  if ( $cmd ne "neighbor_history" ) {

    # if not neighbor_history, write eLinkResult XML instead of ENTREZ_DIRECT

    if ( $num == 0 ) {
      if ( $silent ) {
        return;
      }
    }

    $chunk = 500;
    for ( $start = 0; $start < $num; $start += $chunk ) {

      do_sleep ();

      my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $num, $tool, $email );

      $url = $base . $elink;

      $arg = "dbfrom=$dbase";
      if ( $adddbto ) {
        $arg .= "&db=$dbto";
      }
      $arg .= "&cmd=$cmd&linkname=$name";
      if ( $mode ne "" ) {
        $arg .= "&retmode=$mode";
      }
      $arg .= "&id=";
      $arg .= join ('&id=', @ids);
      if ( $type eq "acc" ) {
        $arg .= "&idtype=acc";
      }

      $data = do_post ($url, $arg, $tool, $email, true);

      print "$data";
    }

    return;
  }

  if ( $batch ) {

    # large list bypass

    batch_elink ( $dbase, $dbto, $name, $web, $key, $num, $stp, "", $lbl, $tool, $email, false );
    return;
  }

  # if not breaking large lists into chunks, use neighbor_history on server

  $url = $base . $elink;

  $arg = "dbfrom=$dbase";
  if ( $adddbto ) {
    $arg .= "&db=$dbto";
  }
  $arg .= "&query_key=$key&WebEnv=$web";
  $arg .= "&cmd=$cmd&linkname=$name";
  if ( $mode ne "" ) {
    $arg .= "&retmode=$mode";
  }
  if ( $type eq "acc" ) {
    $arg .= "&idtype=acc";
  }

  $wb = $web;
  $ky = $key;
  $nm = $num;

  $web = "";
  $key = "";
  $err = "";

  $output = "";

  for ( $tries = 0; $tries < 3 && $web eq ""; $tries++) {

    $output = do_post ($url, $arg, $tool, $email, true);

    $err = $1 if ($output =~ /<Error>(.+?)<\/Error>/i);
    if ( $err eq "" ) {
      $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
      $key = $1 if ($output =~ /<QueryKey>(\S+)<\/QueryKey>/);
    } else {
      if ( ! $silent ) {
        print STDERR "Retrying elink, step $stp: $err\n";
      }
      sleep 3;
    }
  }

  # automatically fail over to batch mode under certain failure conditions

  if ( $stp ne "" ) {
    $stpminusone = $stp - 1;
  }

  if ( $err =~ "^Query failed" or
       $err =~ "^Timeout waiting" or
       $err =~ "^Unable to obtain" or
       $err =~ "^The read request has timed out" ) {
    if ( ! $silent ) {
      print STDERR "$err\nReplicate for debugging with:\n";
      print STDERR "  edirutil -db $dbase -web $wb -key $ky -count $nm";
      if ( $stpminusone > 0 ) {
        print STDERR " -step $stpminusone";
      }
      my $seconds = "";
      my $end_time = Time::HiRes::time();
      my $elapsed = $end_time - $begin_time;
      if ( $elapsed > 0.0005 ) {
        if ( $elapsed =~ /(\d+)\.(\d\d\d)\d+/ ) {
          $seconds = "$1.$2";
        }
      }
      if ( $seconds ne "" ) {
        print STDERR " -seconds $seconds";
      }
      if ( $dbase eq $dbto ) {
        print STDERR " | elink -related\n";
      } else {
        print STDERR " | elink -target $dbto\n";
      }
      print STDERR "Automatically switching to -batch mode\n";
    }

    batch_elink ( $dbase, $dbto, $name, $wb, $ky, $nm, $stp, "", $lbl, $tool, $email, true );
    return;
  }

  # finish reality checks, get count and key, and write ENTREZ_DIRECT structure

  process_history_link ( $arg, $output, $dbase, $dbto, $name, $wb, $web, $key, $stp, $err, $lbl, $tool, $email );
}

# emmdb downloads "Ncbi-mime-asn1 ::= strucseq" ASN.1 from MMDB unique identifiers

sub emmdb {

  # ... | edirect.pl -fetch -format mmdb | ...

  $dbase = shift (@_);
  $web = shift (@_);
  $key = shift (@_);
  $num = shift (@_);
  $id = shift (@_);
  $tool = shift (@_);
  $email = shift (@_);

  if ( $err ne "" ) {
    die "ERROR in mmdb input: $err\n\n";
  }

  $mbase = "https://www.ncbi.nlm.nih.gov/Structure/mmdb/";
  $mprog = "mmdbsrv.cgi";

  if ( $id ne "" ) {

    my @ids = split (',', $id);
    foreach $uid (@ids) {
      $uid =~ s/(\d+)\.\d+/$1/g;
      $url = $mbase . $mprog;

      $arg = "uid=$uid";
      $arg .= "&save=asntext&form=6&db=t&Dopt=j&Complexity=Cn3D%20Subset";

      $mmdb = do_post ($url, $arg, $tool, $email, true);

      print "$mmdb\n";
    }

    return;
  }

  test_edirect ( $dbase, $web, $key, $num, "mmdb" );

  $chunk = 500;
  for ( $start = 0; $start < $num; $start += $chunk ) {

    my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $num, $tool, $email );

    foreach $uid (@ids) {
      $url = $mbase . $mprog;

      $arg = "uid=$uid";
      $arg .= "&save=asntext&form=6&db=t&Dopt=j&Complexity=Cn3D%20Subset";

      $mmdb = do_post ($url, $arg, $tool, $email, true);

      print "$mmdb\n";
    }
  }
}

# entfy sends e-mail

my $ntfy_help = qq{
  -email    Contact person's address
  -tool     Name of script or program

};

sub entfy {

  # ... | edirect.pl -notify

  clearflags ();

  MyGetOptions(
    $ntfy_help,
    "api_key=s" => \$api_key,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "internal" => \$internal,
    "log" => \$log,
    "compact" => \$compact,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "enotify $version\n";
    print $ntfy_help;
    return;
  }

  ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();

  read_aliases ();
  adjust_base ();

  if ( $err ne "" ) {
    die "ERROR in notify input: $err\n\n";
  }

  test_edirect ( $dbase, $web, $key, $num, "notify" );

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }
  if ( $email eq "" ) {
    $email = get_email ();
  }
  if ( $email eq "" ) {
    die "Email value not found in notify input\n";
  }

  binmode STDOUT, ":utf8";

  if ( $num > 0 ) {
    $chunk = 500;
    for ( $start = 0; $start < $num; $start += $chunk ) {

      my @ids = get_uids ( $dbase, $web, $key, $start, $chunk, $num, $tool, $email );

      foreach $uid (@ids) {
        $txt = "echo \"https://www.ncbi.nlm.nih.gov/$dbase/$uid\n\"";
        $str = "mail -s \"A new $dbase record is in Entrez\" $email";
        system "$txt | $str";
      }
    }
  }
}

# epost uploads UIDs or accessions

sub post_chunk {

  my $dbsx = shift (@_);
  my $webx = shift (@_);
  my $keyx = shift (@_);
  my $tulx = shift (@_);
  my $emlx = shift (@_);
  my $uids = shift (@_);
  my $qryx = shift (@_);

  $url = $base;
  $arg = "";
  $wb = $webx;

  if ( $uids ne "" ) {

    $url .= $epost;

    $arg = "db=$dbsx";
    if ( $web ne "" ) {
      $arg .= "&WebEnv=$webx";
    }
    $arg .= "&id=$uids";

  } elsif ( $qryx ne "" ) {

    $url .= $esearch;

    $qryx = map_labels ($qryx);
    $qryx = map_macros ($qryx);
    $enc = uri_escape($qryx);

    $arg = "db=$dbsx&term=$enc";
    if ( $web ne "" ) {
      $arg .= "&WebEnv=$webx";
    }
    $arg .= "&retmax=0&usehistory=y";
    if ( $rldate > 0 ) {
      $arg .= "&reldate=$rldate";
      if ( $dttype eq "" ) {
        $dttype = "PDAT";
      }
    }
    if ( $mndate ne "" and $mxdate ne "" ) {
      $arg .= "&mindate=$mndate&maxdate=$mxdate";
      if ( $dttype eq "" ) {
        $dttype = "PDAT";
      }
    }
    if ( $dttype ne "" ) {
      $arg .= "&datetype=$dttype";
    }
  }

  $output = do_post ($url, $arg, $tulx, $emlx, true);

  $webx = "";
  $keyx = "";
  $err = "";

  $webx = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
  $keyx = $1 if ($output =~ /<QueryKey>(\S+)<\/QueryKey>/);
  $err = $1 if ($output =~ /<Error>(.+?)<\/Error>/i);

  if ( $err ne "" ) {
    write_edirect ( "", "", "", "", "", $err, "", "" );
    close (STDOUT);
    die "ERROR in post output: $err\nURL: $arg\n\n";
  }

  if ( $webx eq "" ) {
    die "WebEnv value not found in post output\n";
  }
  if ( $keyx eq "" ) {
    die "QueryKey value not found in post output\n";
  }

  if ( $wb ne "" and $webx ne $wb ) {
    $err = "WebEnv mismatch in post output - WebEnv1 $wb, WebEnv2 $webx";
    write_edirect ( "", $wb, "", "", "", $err, "", "" );
    close (STDOUT);
    die "WebEnv value changed in post output - WebEnv1 $wb\nWebEnv2 $webx\n";
  }

  return $webx, $keyx;
}

my $post_help = qq{
  -db        Database name
  -id        Unique identifier(s) or accession number(s)
  -format    uid or acc
  -input     Read from file instead of stdin
  -label     Alias for query step

};

sub epost {

  # ... | edirect.pl -post -db nucleotide -format uid | ...

  clearflags ();

  MyGetOptions(
    $post_help,
    "db=s" => \$db,
    "id=s" => \$id,
    "format=s" => \$field,
    "input=s" => \$input,
    "label=s" => \$lbl,
    "api_key=s" => \$api_key,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "internal" => \$internal,
    "log" => \$log,
    "compact" => \$compact,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "epost $version\n";
    print $post_help;
    return;
  }

  # convert spaces between UIDs to commas

  $id =~ s/ /,/g;
  $id =~ s/,,/,/g;

  if ( -t STDIN and not @ARGV ) {
    if ( $id ne "" ) {
      if ( $id =~ /[,0-9.]/ ) {
        $just_num = true;
      }
    }
  } else {
    ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();
  }

  read_aliases ();
  adjust_base ();

  if ( $err ne "" ) {
    die "ERROR in post input: $err\n\n";
  }

  if ( $dbase eq "" ) {
    $dbase = $db;
  }

  $dbase = lc($dbase);

  if ( $dbase eq "" ) {
    die "Must supply -db database on command line\n";
  }

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }

  binmode STDOUT, ":utf8";

  if ( $field eq "" ) {
    $field = "UID";
  }

  # read data from input file instead of piping from stdin
  if ( $input ne "" ) {
    if (open (my $FILE_IN, $input)) {
      $has_num = false;
      $all_num = true;
      while ( $thisline = <$FILE_IN> ) {
        $thisline =~ s/\r//;
        $thisline =~ s/\n//;
        $thisline =~ s/^\s+//;
        $thisline =~ s/\s+$//;

        if ( $thisline =~ /^(\d+)$/ ) {
          push (@rest, $1);
          $has_num = true;
        } elsif ( $thisline =~ /^(.+)$/ ) {
          push (@rest, $1);
          $all_num = false;
        }
      }
      close ($FILE_IN);
      if ( $has_num && $all_num ) {
        $just_num = true;
      }
    } else {
      print STDERR "Unable to open input file '$input'\n";
    }
  }

  my $combo = "";
  my $pfx = "";
  my $loops = 0;

  my $accession_mode = false;
  if ( $field eq "ACCN" or $field eq "accn" or $field eq "ACC" or $field eq "acc" ) {
    $accession_mode = true;
  }

  if ( ! $just_num ) {
    if ( $dbase eq "nucleotide" or
         $dbase eq "nuccore" or
         $dbase eq "est" or
         $dbase eq "gss" or
         $dbase eq "protein" ) {
      $accession_mode = true;
      $field = "acc";
    }
  }

  if ( $field eq "UID" or $field eq "uid" ) {

    if ( $id ne "" ) {

      ( $web, $key ) = post_chunk ( $dbase, $web, $key, $tool, $email, $id, "" );

      $combo .= $pfx . "#" . $key;
      $pfx = " OR ";
      $loops++;

    } else {

      while ( @rest ) {
        my @chunk = splice(@rest, 0, 50000);

        $ids = join (',', @chunk);

        # newline to comma conversion for piped data

        $ids =~ s/\n/,/g;

        if ( $ids =~ /[a-zA-Z]/ ) {
          if ( $dbase eq "nucleotide" or
               $dbase eq "nuccore" or
               $dbase eq "est" or
               $dbase eq "gss" or
               $dbase eq "protein" ) {
            $accession_mode = true;
            $field = "acc";
          } else {
            die "Non-numeric value found in post input:$ids\n\n";
          }
        }

        ( $web, $key ) = post_chunk ( $dbase, $web, $key, $tool, $email, $ids, "" );

        $combo .= $pfx . "#" . $key;
        $pfx = " OR ";
        $loops++;
      }
    }

  } else {

    if ( $id ne "" ) {

      my @chunk = split (',', $id);

      if ( $accession_mode ) {
        $query = join (' [ACCN] OR ', @chunk);
        $query .= " [ACCN]";
        $query =~ s/\./_/g;
      } else {
        $query = join (' OR ', @chunk);
        $query .= " [$field]";
      }

      if ( $accession_mode and $dbase eq "assembly" ) {
        $query =~ s/\[ACCN\]/[ASAC]/g;
      }

      ( $web, $key ) = post_chunk ( $dbase, $web, $key, $tool, $email, "", $query );

      $combo .= $pfx . "#" . $key;
      $pfx = " OR ";
      $loops++;

    } else {

      while ( @rest ) {
        my @chunk = splice(@rest, 0, 50000);

        if ( $accession_mode ) {
          $query = join (' [ACCN] OR ', @chunk);
        } else {
          $query = join (' OR ', @chunk);
        }

        if ( $query eq "" ) {
          die "Must pipe data into stdin\n";
        }

        if ( $accession_mode ) {
        } elsif ( ! $just_num ) {
          die "Non-numeric value found in post input\n";
        }

        if ( $accession_mode ) {
          $query .= " [ACCN]";
          $query =~ s/\./_/g;
        } else {
          $query .= " [$field]";
        }

        if ( $accession_mode and $dbase eq "assembly" ) {
          $query =~ s/\[ACCN\]/[ASAC]/g;
        }

        ( $web, $key ) = post_chunk ( $dbase, $web, $key, $tool, $email, "", $query );

        $combo .= $pfx . "#" . $key;
        $pfx = " OR ";
        $loops++;
      }
    }
  }

  if ( $combo eq "" ) {
    $result = 1;
    die "Failure of post to find data to load\n";
  }

  if ( $loops > 1 ) {
    $url = $base . $esearch;

    $enc = uri_escape($combo);
    $arg = "db=$dbase&term=$enc&WebEnv=$web";
    $arg .= "&retmax=0&usehistory=y";

    $output = do_post ($url, $arg, $tool, $email, true);

    $web = "";
    $key = "";
    $err = "";

    $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
    $key = $1 if ($output =~ /<QueryKey>(\S+)<\/QueryKey>/);
    $err = $1 if ($output =~ /<Error>(.+?)<\/Error>/i);
  }

  ( $num, $key ) = get_count ( $dbase, $web, $key, $tool, $email );

  if ( $num eq "" ) {
    die "Count value not found in post output - WebEnv1 $web\n";
  }

  if ( $lbl ne "" and $key ne "" ) {
    $labels{"$lbl"} = "$key";
  }

  write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );
}

# espel performs an ESpell search

my $spell_help = qq{
  -db       Database name
  -query    Query string

};

sub espel {

  # ... | edirect.pl -spell -db pubmed -query "asthmaa OR alergies" | ...

  clearflags ();

  MyGetOptions(
    $spell_help,
    "db=s" => \$db,
    "query=s" => \$query,
    "q=s" => \$query,
    "api_key=s" => \$api_key,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "internal" => \$internal,
    "log" => \$log,
    "compact" => \$compact,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "espell $version\n";
    print $spell_help;
    return;
  }

  # convert spaces between UIDs to commas

  $id =~ s/ /,/g;
  $id =~ s/,,/,/g;

  if ( -t STDIN and not @ARGV ) {
  } else {
    ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();
  }

  read_aliases ();
  adjust_base ();

  if ( $err ne "" ) {
    die "ERROR in spell input: $err\n\n";
  }

  if ( $dbase eq "" ) {
    $dbase = $db;
  }

  $dbase = lc($dbase);

  if ( $dbase eq "" ) {
    die "Must supply -db database on command line\n";
  }

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }

  binmode STDOUT, ":utf8";

  if ( $query eq "" ) {
    die "Must supply -query term expression on command line\n";
  }

  $url = $base . $espell;

  $enc = uri_escape($query);
  $arg = "db=$dbase&term=$enc";

  $data = do_post ($url, $arg, $tool, $email, true);

  Encode::_utf8_on($data);

  print "$data";
}

# ecitmtch performs an ECitMatch search

my $citmatch_help = qq{
  -journal    Journal Title
  -year       Year
  -volume     Volume
  -page       First Page
  -author     Author Name

};

sub ecitmtch {

  # ... | edirect.pl -citmatch -journal "proc natl acad sci u s a" -year 2005 ...

  clearflags ();

  MyGetOptions(
    $citmatch_help,
    "journal=s" => \$journal,
    "year=s" => \$year,
    "volume=s" => \$volume,
    "page=s" => \$page,
    "author=s" => \$author,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "internal" => \$internal,
    "log" => \$log,
    "compact" => \$compact,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "ecitmatch $version\n";
    print $citmatch_help;
    return;
  }

  # convert spaces between UIDs to commas

  $id =~ s/ /,/g;
  $id =~ s/,,/,/g;

  if ( -t STDIN and not @ARGV ) {
  } else {
    ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();
  }

  read_aliases ();
  adjust_base ();

  if ( $err ne "" ) {
    die "ERROR in citation match input: $err\n\n";
  }

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }

  binmode STDOUT, ":utf8";

  $url = $base . $ecitmat;

  $query = "";

  if ( $journal ne "" ) {
    $query .= $journal;
  }
  $query .= "|";
  if ( $year ne "" ) {
    $query .= $year;
  }
  $query .= "|";
  if ( $volume ne "" ) {
    $query .= $volume;
  }
  $query .= "|";
  if ( $page ne "" ) {
    $query .= $page;
  }
  $query .= "|";
  if ( $author ne "" ) {
    $query .= $author;
  }
  $query .= "||";

  $enc = uri_escape($query);
  $arg = "db=pubmed&retmode=xml&bdata=$enc";

  $data = do_post ($url, $arg, $tool, $email, true);

  Encode::_utf8_on($data);

  if ( $data =~ "NOT_FOUND" ) {
    return;
  }

  if ( $data =~ /.+\|AMBIGUOUS (.+)$/ ) {
    my $my_uids = $1;
    my @ids = split (',', $my_uids);

    foreach $uid (@ids) {
      print "$uid\n";
    }

    return;
  }

  if ( $data =~ /.+\|(\d+)$/ ) {
    my $my_uid = $1;
    print "$my_uid\n";
  }
}

# eprxy reads a file of query proxies, can also pipe from stdin

my $prxy_help = qq{
  -alias    File of aliases
  -pipe     Read aliases from stdin

};

sub eprxy {

  # ... | edirect.pl -proxy -alias file_name | ...

  clearflags ();

  MyGetOptions(
    $prxy_help,
    "pipe" => \$pipe,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "internal" => \$internal,
    "log" => \$log,
    "compact" => \$compact,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "eproxy $version\n";
    print $prxy_help;
    return;
  }

  if ( -t STDIN and not @ARGV ) {
  } elsif ( $pipe ) {
  } else {
    ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();
  }

  read_aliases ();
  adjust_base ();

  if ( $pipe ) {
    while ( defined($thisline = <STDIN>) ) {
      $thisline =~ s/\r//;
      $thisline =~ s/\n//;
      $thisline =~ s/ +/ /g;
      $thisline =~ s/> </></g;
      if ( $thisline =~ /(.+)\t(.+)/ ) {
        $ky = $1;
        $vl = $2;
        $vl =~ s/\"//g;
        $macros{"$ky"} = "$vl";
      }
    }
  }

  write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );
}

# esrch performs a new EUtils search, but can read a previous web environment value

my $srch_help = qq{
Query Specification

  -db          Database name
  -query       Query string

Document Order

  -sort        Result presentation order

Date Constraint

  -days        Number of days in the past
  -datetype    Date field abbreviation
  -mindate     Start of date range
  -maxdate     End of date range

Limit by Field

  -field       Query words individually in field
  -pairs       Query overlapping word pairs

Spell Check

  -spell       Correct misspellings in query

Miscellaneous Arguments

  -label       Alias for query step

Sort Order Examples

  -db            -sort
  ___            _____

  gene
                 Chromosome
                 Gene Weight
                 Name
                 Relevance

  geoprofiles
                 Default Order
                 Deviation
                 Mean Value
                 Outliers
                 Subgroup Effect

  pubmed
                 First Author
                 Journal
                 Last Author
                 Pub Date
                 Recently Added
                 Relevance
                 Title

  (sequences)
                 Accession
                 Date Modified
                 Date Released
                 Default Order
                 Organism Name
                 Taxonomy ID

  snp
                 Chromosome Base Position
                 Default Order
                 Heterozygosity
                 Organism
                 SNP_ID
                 Success Rate

};

sub remove_punctuation {

  my $qury = shift (@_);

  $qury =~ s/[^a-zA-Z0-9]/ /g;
  $qury =~ s/ +/ /g;

  return $qury;
}

sub remove_stop_words {

  my $qury = shift (@_);

  # split to protect against regular expression artifacts
  $qury =~ s/[^a-zA-Z0-9]/ /g;
  $qury =~ s/ +/ /g;

  my @words = split (' ', $qury);
  my $kept = "";
  my $pfx = "";

  foreach $term (@words) {

    my $trm = lc($term);
    $trm = "#$trm#";

    if ($stop_words !~ /$trm/) {
      $kept .= "$pfx$term";
      $pfx = " ";
    }
  }

  if ( $kept ne "" ) {
    $qury = $kept;
  }

  return $qury;
}

sub field_each_word {

  my $fld = shift (@_);
  my $qury = shift (@_);

  $qury =~ s/,/ /g;
  $qury =~ s/ +/ /g;

  my @words = split (' ', $qury);
  $qury = "";
  my $pfx = "";

  foreach $term (@words) {
    $qury .= "$pfx$term [$fld]";
    $pfx = " AND ";
  }

  return $qury;
}

sub merge_each_word {

  my $fld = shift (@_);
  my $qury = shift (@_);

  $qury =~ s/,/ /g;
  $qury =~ s/ +/ /g;

  my @words = split (' ', $qury);
  $qury = "";
  my $pfx = "";

  foreach $term (@words) {
    $qury .= "$pfx$term [$fld]";
    $pfx = " OR ";
  }

  return $qury;
}

sub field_each_pair {

  my $fld = shift (@_);
  my $qury = shift (@_);

  my @words = split (' ', $qury);
  $qury = "";
  my $pfx = "";

  my $prev = "";
  foreach $term (@words) {
    my $trm = lc($term);
    $trm = "#$trm#";
    if ($stop_words =~ /$trm/) {
      $prev = "";
      $term = "";
    }
    if ( $prev ne "" ) {
      $qury .= "$pfx\"$prev $term\" [$fld]";
      $pfx = " AND ";
    }
    $prev = $term;
  }

  return $qury;
}

sub esrch {

  # ... | edirect.pl -search -db nucleotide -query "M6506* [ACCN] OR 1322283 [UID]" -days 365 | ...

  clearflags ();

  MyGetOptions(
    $srch_help,
    "db=s" => \$db,
    "query=s" => \$query,
    "q=s" => \$query,
    "sort=s" => \$sort,
    "days=i" => \$rldate,
    "mindate=s" => \$mndate,
    "maxdate=s" => \$mxdate,
    "datetype=s" => \$dttype,
    "label=s" => \$lbl,
    "journal=s" => \$journal,
    "pub=s" => \$pub,
    "released=s" => \$released,
    "country=s" => \$country,
    "feature=s" => \$feature,
    "location=s" => \$location,
    "molecule=s" => \$molecule,
    "organism=s" => \$organism,
    "source=s" => \$source,
    "status=s" => \$status,
    "type=s" => \$gtype,
    "class=s" => \$class,
    "kind=s" => \$kind,
    "pathway=s" => \$pathway,
    "clean" => \$clean,
    "field=s" => \$field,
    "word" => \$word,
    "drop" => \$drop,
    "trim" => \$trim,
    "trunc" => \$trunc,
    "spell" => \$spell,
    "split=s" => \$split,
    "merge=s" => \$meadow,
    "pairs=s" => \$pair,
    "api_key=s" => \$api_key,
    "email=s" => \$emaddr,
    "tool=s" => \$tuul,
    "help" => \$help,
    "silent" => \$silent,
    "verbose" => \$verbose,
    "debug" => \$debug,
    "internal" => \$internal,
    "log" => \$log,
    "compact" => \$compact,
    "http=s" => \$http,
    "https=s" => \$http,
    "alias=s" => \$alias,
    "base=s" => \$basx
  );

  if ( $help ) {
    print "esearch $version\n";
    print $srch_help;
    return;
  }

  if ( -t STDIN and not @ARGV ) {
  } else {
    ( $dbase, $web, $key, $num, $stp, $err, $tool, $email, $just_num, @rest ) = read_edirect ();
  }

  read_aliases ();
  adjust_base ();

  if ( $err ne "" ) {
    die "ERROR in search input: $err\n\n";
  }

  if ( $db ne "" ) {
    $dbase = $db;
  }

  $dbase = lc($dbase);

  if ( $dbase eq "" ) {
    die "Must supply -db database on command line\n";
  }

  if ( $tuul ne "" ) {
    $tool = $tuul;
  }
  if ( $emaddr ne "" ) {
    $email = $emaddr;
  }

  binmode STDOUT, ":utf8";

  # support all efilter shortcut flags in esearch (undocumented)
  $query = process_extras ( $query, $pub, $released, $journal, $country, $feature, $location, $molecule, $organism, $source, $status, $gtype, $class, $kind, $pathway );

  if ( $query eq "" ) {
    die "Must supply -query search expression on command line\n";
  }

  $url = $base . $esearch;

  $query = map_labels ($query);
  $query = map_macros ($query);

  # multi-step query cleaning (undocumented)
  if ( $clean ) {
    if ( $query =~ /^(.*)\(.+\)(.*)$/ ) {
      $query = "$1 $2";
    }
    if ( $query =~ /^ +(.+)$/ ) {
      $query = $1;
    }
    if ( $query =~ /^(.+) +$/ ) {
      $query = $1;
    }
    $query =~ s/ +/ /g;
    $query = remove_stop_words ($query);
    $query = spell_check_query ($dbase, $query);
  }

  # remove punctuation from query (undocumented)
  if ( $word ) {
    $query = remove_punctuation ($query);
  }

  # drop stop words from query (undocumented)
  if ( $drop ) {
    $query = remove_stop_words ($query);
  }

  # trim words within parentheses (undocumented)
  if ( $trim ) {
    if ( $query =~ /^(.*)\(.+\)(.*)$/ ) {
      $query = "$1 $2";
    }
  }

  # truncate words at first parenthesis (undocumented)
  if ( $trunc ) {
    if ( $query =~ /^(.+)\(.*$/ ) {
      $query = $1;
    }
  }

  # remove leading, trailing, and multiple spaces
  if ( $query =~ /^ +(.+)$/ ) {
    $query = $1;
  }
  if ( $query =~ /^(.+) +$/ ) {
    $query = $1;
  }
  $query =~ s/ +/ /g;

  # spell check each query word
  if ( $spell ) {
    $query = spell_check_query ($dbase, $query);
  }

  # force each query word to be separately fielded, combined with AND (undocumented)
  if ( $split ne "" ) {
    $query = field_each_word ($split, $query);
  }

  # force each query word to be separately fielded, combined with OR (undocumented)
  if ( $meadow ne "" ) {
    $query = merge_each_word ($meadow, $query);
  }

  # -field combines -drop and -split (-field TITL produces same behavior as Web PubMed)
  if ( $field ne "" ) {
    $query = remove_stop_words ($query);
    $query = field_each_word ($field, $query);
  }

  # -pairs separately fields query word pairs, breaking chain at stop words
  if ( $pair ne "" ) {
    $query = remove_punctuation ($query);
    if ( $query =~ /^ +(.+)$/ ) {
      $query = $1;
    }
    if ( $query =~ /^(.+) +$/ ) {
      $query = $1;
    }
    $query =~ s/ +/ /g;
    $query = field_each_pair ($pair, $query);
  }

  $enc = uri_escape($query);
  $arg = "db=$dbase&term=$enc";
  if ( $web ne "" ) {
    $arg .= "&WebEnv=$web";
  }

  if ( $sort ne "" ) {
    if ( $sort eq "Relevance" ) {
      $sort = "relevance";
    }
    $arg .= "&sort=$sort";
  }
  $arg .= "&retmax=0&usehistory=y";

  if ( $rldate > 0 ) {
    $arg .= "&reldate=$rldate";
    if ( $dttype eq "" ) {
      $dttype = "PDAT";
    }
  }
  if ( $mndate ne "" and $mxdate ne "" ) {
    $arg .= "&mindate=$mndate&maxdate=$mxdate";
    if ( $dttype eq "" ) {
      $dttype = "PDAT";
    }
  }
  if ( $dttype ne "" ) {
    $arg .= "&datetype=$dttype";
  }

  $wb = $web;

  $web = "";
  $key = "";
  $num = "";
  $err = "";
  my $trn = "";

  $output = "";

  for ( $tries = 0; $tries < 3 && $web eq ""; $tries++) {
    $output = do_post ($url, $arg, $tool, $email, true);

    $err = $1 if ($output =~ /<Error>(.+?)<\/Error>/i);
    if ( $err eq "" ) {
      $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
      $key = $1 if ($output =~ /<QueryKey>(\S+)<\/QueryKey>/);
      $num = $1 if ($output =~ /<Count>(\S+)<\/Count>/);
      $trn = $1 if ($output =~ /<QueryTranslation>(.+?)<\/QueryTranslation>/i);
    } else {
      if ( ! $silent ) {
        print STDERR "Retrying esearch, step $stp: $err\n";
      }
      sleep 3;
    }
  }

  if ( $err ne "" ) {
    write_edirect ( "", "", "", "", "", $err, "", "" );
    close (STDOUT);
    die "ERROR in search output: $err\nURL: $arg\n\n";
  }

  if ( $web eq "" ) {
    die "WebEnv value not found in search output - WebEnv1 $wb\n";
  }
  if ( $key eq "" ) {
    die "QueryKey value not found in search output - WebEnv1 $wb\n";
  }

  if ( $wb ne "" and $web ne $wb ) {
    $err = "WebEnv mismatch in search output - WebEnv1 $wb, WebEnv2 $web";
    write_edirect ( "", $wb, "", "", "", $err, "", "" );
    close (STDOUT);
    die "WebEnv value changed in search output - WebEnv1 $wb\nWebEnv2 $web\n";
  }

  if ( $num eq "" ) {
    die "Count value not found in search output - WebEnv1 $web\n";
  }

  if ( $lbl ne "" and $key ne "" ) {
    $labels{"$lbl"} = "$key";
  }

  write_edirect ( $dbase, $web, $key, $num, $stp, $err, $tool, $email );

  if ( $verbose ) {
    my $seconds = "";
    my $end_time = Time::HiRes::time();
    my $elapsed = $end_time - $begin_time;
    if ( $elapsed > 0.0005 ) {
      if ( $elapsed =~ /(\d+)\.(\d\d\d)\d+/ ) {
        $seconds = "$1.$2";
      }
    }
    if ( $seconds ne "" ) {
      print STDERR "Elapsed time is $seconds seconds\n";
    }
  }

  if ( $log ) {
    if ( $trn ne "" ) {
      print STDERR "$trn\n";
    }
  }
}

#  eaddr returns the current user's e-mail address

sub eaddr {
  my $addr = get_email ();
  print "$addr\n";
}

#  eblst front-end to Blast API must be fed chunks of FASTA sequences from stdin

sub eblst {

  clearflags ();

  $encoded_query = "";

  while ( defined($thisline = <STDIN>) ) {
    $encoded_query = $encoded_query . uri_escape($thisline);
  }

  $url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi";
  $arg = "CMD=Put&PROGRAM=blastp&DATABASE=nr&QUERY=" . $encoded_query;

  $ua = new LWP::UserAgent (timeout => 300);
  $ua->env_proxy;

  $req = new HTTP::Request POST => "$url";
  $req->content_type('application/x-www-form-urlencoded');
  $req->content("$arg");

  do_sleep ();

  $response = $ua->request ( $req );

  # parse out the request id
  $response->content =~ /^    RID = (.*$)/m;
  $rid=$1;

  if ( $rid eq "" ) {
    print STDERR "Unable to create BLAST request.\n";
    exit 1;
  }

  print STDERR "RID:  $rid\n";

  # parse out the estimated time to completion
  $response->content =~ /^    RTOE = (.*$)/m;
  $rtoe=$1;

  # wait for search to complete
  if ( $rtoe ne "" ) {
      sleep $rtoe;
  }

  # poll for results
  while (true)
    {
    $req = new HTTP::Request GET => "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid";

    do_sleep ();

    $response = $ua->request($req);

    if ($response->content =~ /\s+Status=WAITING/m)
        {
        print STDERR "Searching...\n";
        sleep 20;
        next;
        }

    if ($response->content =~ /\s+Status=FAILED/m)
        {
        print STDERR "Search $rid failed; please report to blast-help\@ncbi.nlm.nih.gov.\n";
        exit 4;
        }

    if ($response->content =~ /\s+Status=UNKNOWN/m)
        {
        print STDERR "Search $rid expired.\n";
        exit 3;
        }

    if ($response->content =~ /\s+Status=READY/m) 
        {
        if ($response->content =~ /\s+ThereAreHits=yes/m)
            {
            print STDERR "Search complete, retrieving results...\n";
            last;
            }
        else
            {
            print STDERR "No hits found.\n";
            exit 2;
            }
        }

        # if we get here, something unexpected happened.
        print STDERR "Unexpected situation.\n";
        exit 5;
    } # end poll loop

  # retrieve and display results
  $req = new HTTP::Request GET => "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID=$rid";

  do_sleep ();

  $response = $ua->request($req);

  print $response->content;
}

sub ftls {

  @args = @ARGV;
  $max = scalar @args;
  if ( $max < 2 ) {
    die "Needs SERVER and PATH arguments\n";
  }
  my $server = $args[0];
  my $dir    = $args[1];

  if ( $server =~ /^ftp:\/\/(.+)/ ) {
    $server = $1;
  }

  my $ftp    = new Net::FTP($server, Passive => 1)
    or die "Unable to connect to FTP server: $!";

  $ftp->login or die "Unable to log in to FTP server: ", $ftp->message;
  $ftp->cwd($dir) or die "Unable to change to $dir: ", $ftp->message;
  my $contents = $ftp->dir;
  die "Unable to list contents" unless defined $contents;

  for (@$contents) {
    if (/^-.*?(\S*)$/) {
        print "$1\n";
    } elsif (/^d.*?(\S*)$/) {
        print "$1/\n";
    } elsif (/^l.*?(\S*) -> \S*$/) {
        print "$1@\n";
    }
  }
}

sub asls {

  @args = @ARGV;
  $max = scalar @args;
  if ( $max < 1 ) {
    die "Needs PATH argument\n";
  }
  my $server = "ftp.ncbi.nlm.nih.gov";
  my $dir    = $args[0];

  my $ftp    = new Net::FTP($server, Passive => 1)
    or die "Unable to connect to FTP server: $!";

  $ftp->login or die "Unable to log in to FTP server: ", $ftp->message;
  $ftp->cwd($dir) or die "Unable to change to $dir: ", $ftp->message;
  my $contents = $ftp->dir;
  die "Unable to list contents" unless defined $contents;

  for (@$contents) {
    if (/^-.*?(\S*)$/) {
        print "$1\n";
    } elsif (/^d.*?(\S*)$/) {
        print "$1/\n";
    } elsif (/^l.*?(\S*) -> \S*$/) {
        print "$1@\n";
    }
  }
}

sub ftcp {

  @args = @ARGV;
  $max = scalar @args;
  if ( $max < 2 ) {
    die "Needs SERVER and PATH arguments\n";
  }
  my $server = $args[0];
  my $dir    = $args[1];

  if ( $server =~ /^ftp:\/\/(.+)/ ) {
    $server = $1;
  }

  my $ftp    = new Net::FTP($server, Passive => 1)
    or die "Unable to connect to FTP server: $!";

  my @failed = ();

  $ftp->login or die "Unable to log in to FTP server: ", $ftp->message;
  $ftp->cwd($dir) or die "Unable to change to $dir: ", $ftp->message;
  $ftp->binary or warn "Unable to set binary mode: ", $ftp->message;

  if ($max > 2) {
  # file names on command line
    for ( $i = 2; $i < $max; $i++) {
      my $fl = $args[$i];
      if (! -e $fl) {
        if (! $ftp->get($fl) ) {
          my $msg = $ftp->message;
          chomp $msg;
          push (@failed, "$fl ($msg)");
        }
      }
    }
  } elsif ( -t STDIN ) {
    print STDERR "\nNO INPUT PIPED FROM STDIN\n\n";
  } else {
  # read file names from stdin
    while ( <STDIN> ) {
      chomp;
      $_ =~ s/\r$//;
      print "$_\n";
      my $fl = $_;
      if (! -e $fl) {
        if (! $ftp->get($fl) ) {
          my $msg = $ftp->message;
          chomp $msg;
          push (@failed, "$fl ($msg)");
        }
      }
    }
  }

  if (@failed) {
    my $errs = join ("\n", @failed);
    print STDERR "\nFAILED TO DOWNLOAD:\n\n$errs\n";
    exit 1;
  }
}

# convert_json based on convert_bool from:
#
# https://stackoverflow.com/questions/41039792correct-and-easy-way-convert-jsonppboolean-to-0-1-with-perl/41040316#41040316
#
# with hash section improvements by Aaron Ucko to convert non-XML-legal element characters to underscore.

sub convert_json {
  my %unrecognized;

  local *_convert_json = sub {
    my $ref_type = ref($_[0]);
    if (!$ref_type) {
      # Nothing.
    }
    elsif ($ref_type eq 'HASH') {
      my %converted = ();
      while (my ($k, $v) = each %{ $_[0] }) {
        $k =~ s/[^[:alnum:]_:.-]|^[0-9.-]/_/g;
        _convert_json($v);
        $converted{$k} = $v;
      }
      %{ $_[0] } = %converted;
    }
    elsif ($ref_type eq 'ARRAY') {
      _convert_json($_) for @{ $_[0] };
    }
    elsif (
       $ref_type eq 'JSON::PP::Boolean' || $ref_type eq 'Types::Serialiser::Boolean'
    ) {
      $_[0] = $_[0] ? 1 : 0;
    }
    else {
      ++$unrecognized{$ref_type};
    }
  };

  &_convert_json;
}

my $transmute_help = qq{
  [escape|unescape|encode64|decode64|plain|simple|script|pretty|docsum|pubmed]

};

sub tmut {

  @args = @ARGV;
  $max = scalar @args;
  if ( $max < 1 ) {
    die "Must supply conversion type on command line\n";
  }
  # read required function argument
  my $type = $args[0];
  my $obj = "";
  my $spt = "";
  if ( $max > 1 ) {
    # read optional parent object name
    $obj = $args[1];
  }
  if ( $max > 2 ) {
    # read optional json tag for splitting
    $spt = $args[2];
  }

  if ( $type eq "-help" ) {
    print "transmute $version\n";
      print $transmute_help;
    return;
  }

  # read entire XML input stream into a single string

  my $holdTerminator = $/;
  undef $/;
  my $data = <STDIN>;
  $/ = $holdTerminator;

  # exit on empty data
  if ( $data eq "" ) {
    exit 1;
  }

  # perform specific conversions

  if ( $type eq "unescape" || $type eq "-unescape" ) {

    $data = uri_unescape($data);

    # convert plus signs to spaces
    $data =~ s/\+/ /g;

    # compress runs of spaces
    $data =~ s/ +/ /g;

    print "$data";
  }

  if ( $type eq "escape" || $type eq "-escape" ) {

    # compress runs of spaces
    $data =~ s/ +/ /g;

    $data = uri_escape($data);

    print "$data";
  }

  if ( $type eq "decode64" || $type eq "-decode64" ) {

    $data = decode_base64($data);

    print "$data";
  }

  if ( $type eq "encode64" || $type eq "-encode64" ) {

    $data = encode_base64($data);

    print "$data";
  }

  if ( $type eq "plain" || $type eq "-plain" ) {

    # remove embedded mixed-content tags
    $data =~ s/<[^>]*>//g;

    # compress runs of spaces
    $data =~ s/ +/ /g;

    print "$data";
  }

  if ( $type eq "simple" || $type eq "-simple" ) {

    # remove embedded mixed-content tags and everything in between
    $data =~ s,<[^>]*/>,,g;
    $data =~ s,<(\S+)[^>]*>.*?</\1>,,g;

    # compress runs of spaces
    $data =~ s/ +/ /g;

    print "$data";
  }

  if ( $type eq "script" || $type eq "-script" ) {

    # remove newlines, tabs, space between tokens, compress runs of spaces
    $data =~ s/\r/ /g;
    $data =~ s/\n/ /g;
    $data =~ s/\t//g;
    $data =~ s/ +/ /g;
    $data =~ s/> +</></g;

    # remove embedded script tags
    $data =~ s|<script.*?</script>||g;

    # compress runs of spaces
    $data =~ s/ +/ /g;

    # restore newlines between objects
    $data =~ s/> *?</>\n</g;

    print "$data";
  }

  if ( $type eq "pubmed" || $type eq "-pubmed" ) {

    # remove newlines, tabs, space between tokens, compress runs of spaces
    $data =~ s/\r/ /g;
    $data =~ s/\n/ /g;
    $data =~ s/\t//g;
    $data =~ s/ +/ /g;
    $data =~ s/> +</></g;

    # my $markup = '(?:[biu]|su[bp])';
    # my $attrs = ' ?';
    my $markup = '(?:[\w.:_-]*:)?[[:lower:]-]+';
    my $attrs = '(?:\s[\w.:_-]+=[^>]*)?';

    # check for possible newline artifact
    $data =~ s|</$markup>\n||g;
    $data =~ s|\n<$markup$attrs>||g;

    # remove mixed content tags
    $data =~ s|</$markup>||g;
    $data =~ s|<$markup$attrs/?>||g;
    $data =~ s|</?DispFormula$attrs>| |g;

    # check for encoded tags
    if ( $data =~ /\&amp\;/ || $data =~ /\&lt\;/ || $data =~ /\&gt\;/ ) {
      # remove runs of amp
      $data =~ s|&amp;(?:amp;)+|&amp;|g;
      # fix secondary encoding
      $data =~ s|&amp;lt;|&lt;|g;
      $data =~ s|&amp;gt;|&gt;|g;
      $data =~ s|&amp;#(\d+);|&#$1;|g;
      # temporarily protect encoded scientific symbols, e.g., PMID 9698410 and 21892341
      $data =~ s|(?<= )(&lt;)(=*$markup&gt;)(?= )|$1=$2|g;
      # remove encoded markup
      $data =~ s|&lt;/$markup&gt;||g;
      $data =~ s|&lt;$markup$attrs/?&gt;||g;
      # undo temporary protection of scientific symbols adjacent to space
      $data =~ s|(?<= )(&lt;)=(=*$markup&gt;)(?= )|$1$2|g;
    }

    # compress runs of horizontal whitespace
    $data =~ s/\h+/ /g;

    # remove lines with just space
    $data =~ s/\n \n/\n/g;

    # remove spaces just outside of angle brackets
    $data =~ s|> |>|g;
    $data =~ s| <|<|g;

    # remove spaces just inside of parentheses
    $data =~ s|\( |\(|g;
    $data =~ s| \)|\)|g;

    # remove newlines flanking spaces
    $data =~ s|\n ||g;
    $data =~ s| \n| |g;

    # restore newlines between objects
    $data =~ s/> *?</>\n</g;

    print "$data\n";
  }

  if ( $type eq "docsum" || $type eq "-docsum" ) {

    # remove newlines, tabs, space between tokens, compress runs of spaces
    $data =~ s/\r/ /g;
    $data =~ s/\n/ /g;
    $data =~ s/\t//g;
    $data =~ s/ +/ /g;
    $data =~ s/> +</></g;

    # move UID from attribute to object
    if ($data !~ /<Id>\d+<\/Id>/) {
      $data =~ s/<DocumentSummary uid=\"(\d+)\">/<DocumentSummary><Id>$1<\/Id>/g;
    }
    $data =~ s/<DocumentSummary uid=\"\d+\">/<DocumentSummary>/g;

    # fix bad encoding
    my @accum = ();
    my @working = ();
    my $prefix = "";
    my $suffix = "";
    my $docsumset_attrs = '';

    if ( $data =~ /(.+?)<DocumentSummarySet(\s+.+?)?>(.+)<\/DocumentSummarySet>(.+)/s ) {
      $prefix = $1;
      $docsumset_attrs = $2;
      my $docset = $3;
      $suffix = $4;

      my @vals = ($docset =~ /<DocumentSummary>(.+?)<\/DocumentSummary>/sg);
      foreach $val (@vals) {
        push (@working, "<DocumentSummary>");
        if ( $val =~ /<Title>(.+?)<\/Title>/ ) {
          my $x = $1;
          if ( $x =~ /\&amp\;/ || $x =~ /\&lt\;/ || $x =~ /\&gt\;/ || $x =~ /\</ || $x =~ /\>/ ) {
            while ( $x =~ /\&amp\;/ || $x =~ /\&lt\;/ || $x =~ /\&gt\;/ ) {
              HTML::Entities::decode_entities($x);
            }
            # removed mixed content tags
            $x =~ s|<b>||g;
            $x =~ s|<i>||g;
            $x =~ s|<u>||g;
            $x =~ s|<sup>||g;
            $x =~ s|<sub>||g;
            $x =~ s|</b>||g;
            $x =~ s|</i>||g;
            $x =~ s|</u>||g;
            $x =~ s|</sup>||g;
            $x =~ s|</sub>||g;
            $x =~ s|<b/>||g;
            $x =~ s|<i/>||g;
            $x =~ s|<u/>||g;
            $x =~ s|<sup/>||g;
            $x =~ s|<sub/>||g;
            # Reencode any resulting less-than or greater-than entities to avoid breaking the XML.
            $x =~ s/</&lt;/g;
            $x =~ s/>/&gt;/g;
            $val =~ s/<Title>(.+?)<\/Title>/<Title>$x<\/Title>/;
          }
        }
        if ( $val =~ /<Summary>(.+?)<\/Summary>/ ) {
          my $x = $1;
          if ( $x =~ /\&amp\;/ ) {
            HTML::Entities::decode_entities($x);
            # Reencode any resulting less-than or greater-than entities to avoid breaking the XML.
            $x =~ s/</&lt;/g;
            $x =~ s/>/&gt;/g;
            $val =~ s/<Summary>(.+?)<\/Summary>/<Summary>$x<\/Summary>/;
          }
        }
        push (@working, $val );
        push (@working, "</DocumentSummary>");
      }
    }

    if ( scalar @working > 0 ) {
      push (@accum, $prefix);
      push (@accum, "<DocumentSummarySet$docsumset_attrs>");
      push (@accum, @working);
      push (@accum, "</DocumentSummarySet>");
      push (@accum, $suffix);
      $data = join ("\n", @accum);
      $data =~ s/\n\n/\n/g;
    }

    # restore newlines between objects
    $data =~ s/> *?</>\n</g;

    print "$data\n";
  }

  if ( $type eq "pretty" || $type eq "-pretty" ) {

    # remove newlines, tabs, space between tokens, compress runs of spaces
    $data =~ s/\r/ /g;
    $data =~ s/\n/ /g;
    $data =~ s/\t//g;
    $data =~ s/ +/ /g;
    $data =~ s/> +</></g;

    # restore newlines between objects
    $data =~ s/> *?</>\n</g;

    print "$data\n";
  }

  if ( $type eq "json2xml" || $type eq "-json2xml" || $type eq "j2x" || $type eq "-j2x" ) {

    # convert JSON to XML

    my @item = ();

    if ( defined($spt) && $spt ne "" ) {
      @items = split(/(?<=\})\s*(?=\{"\Q$spt\E":)/, $data);
    } else {
      push (@items, $data);
    }

    binmode(STDOUT, ":utf8");

    if ( defined($obj) && $obj ne "" ) {
      print "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
      print "<!DOCTYPE $obj>\n";
      print "<$obj>\n";
    }

    foreach $item (@items) {
      my $jc = JSON::PP->new->ascii->pretty->allow_nonref;
      my $conv = $jc->decode($item);
      convert_json($conv);
      my $result = XMLout($conv, SuppressEmpty => undef);

      # remove newlines, tabs, space between tokens, compress runs of spaces
      $result =~ s/\r/ /g;
      $result =~ s/\n/ /g;
      $result =~ s/\t//g;
      $result =~ s/ +/ /g;
      $result =~ s/> +</></g;

      # remove <opt> flanking object
      if ( $result =~ /<opt>\s*?</ and $result =~ />\s*?<\/opt>/ ) {
        $result =~ s/<opt>\s*?</</g;
        $result =~ s/>\s*?<\/opt>/>/g;
      }

      # restore newlines between objects
      $result =~ s/> *?</>\n</g;

      print "$result\n";
    }

    if ( defined($obj) && $obj ne "" ) {
      print "</$obj>\n";
    }
  }

  if ( $type eq "xml2json" || $type eq "-xml2json" || $type eq "x2j" || $type eq "-x2j" ) {

    # convert XML to JSON

    my $xc = new XML::Simple(KeepRoot => 1);
    my $conv = $xc->XMLin($data);
    convert_json($conv);
    my $jc = JSON::PP->new->ascii->pretty->allow_nonref;
    my $result = $jc->encode($conv);

    $data = "$result";

    print "$data\n";
  }
}

# send actual query

sub do_nquire_post {

  $urlx = shift (@_);
  $argx = shift (@_);

  $rslt = "";

  if ( $debug ) {
    if ( $argx ne "" ) {
      print STDERR "URL: $urlx?$argx\n\n";
    } else {
      print STDERR "URL: $urlx\n\n";
    }
  }

  if ( $http eq "get" or $http eq "GET" ) {
    if ( $argx ne "" ) {
      $urlx .= "?";
      $urlx .= "$argx";
    }

    $usragnt = new LWP::UserAgent (timeout => 300);
    $usragnt->agent( "$agent" );

    $res = $usragnt->get ( $urlx );

    if ( $res->is_success) {
      $rslt = $res->content;
    } elsif ( $debug ) {
      print STDERR "STATUS: " . $res->status_line . "\n";
    }

    if ( $rslt eq "" and $debug ) {
      print STDERR "No do_get output returned from '$urlx'\n";
    }

    if ( $debug ) {
      print STDERR "$rslt\n";
    }

    return $rslt;
  }

  $usragnt = new LWP::UserAgent (timeout => 300);
  $usragnt->agent( "$agent" );

  $req = new HTTP::Request POST => "$urlx";
  $req->content_type('application/x-www-form-urlencoded');
  $req->content("$argx");

  $res = $usragnt->request ( $req );

  if ( $res->is_success) {
    $rslt = $res->content;
  } elsif ( $debug ) {
    print STDERR "STATUS: " . $res->status_line . "\n";
  }

  if ( $rslt eq "" && $debug ) {
    if ( $argx ne "" ) {
      $urlx .= "?";
      $urlx .= "$argx";
    }
    print STDERR "No do_post output returned from '$urlx'\n";
  }

  if ( $debug ) {
    print STDERR "$rslt\n";
  }

  return $rslt;
}

# uri_escape with backslash exceptions

sub do_uri_escape {

  $patx = shift (@_);

  $rslt = "";

  while ( $patx ne "" ) {
    if ( $patx =~ /^\\\\(.+)/ ) {
      $rslt .= "\\";
      $patx = $1;
    } elsif ( $patx =~ /^\\(.)(.+)/ ) {
      $rslt .= $1;
      $patx = $2;
    } elsif ( $patx =~ /^(.)(.+)/ ) {
      $rslt .= uri_escape ($1);
      $patx = $2;
    } elsif ( $patx =~ /^(.)/ ) {
      $rslt .= uri_escape ($1);
      $patx = "";
    }
  }

  return $rslt;
}

# nquire executes an external URL query from command line arguments

my $nquire_help = qq{
Query Commands

  -ftp         Uses FTP instead of HTTP
  -get         Uses HTTP GET instead of POST
  -url         Base URL for external search

Documentation

  -help        Print this document
  -examples    Examples of advanced queries
  -version     Print version number

Examples

  nquire -get "http://collections.mnh.si.edu/services/resolver/resolver.php" \\
    -voucher "Birds:625456" |
  xtract -pattern Result -element ScientificName Country

  nquire -get http://w1.weather.gov/xml/current_obs/KSFO.xml |
  xtract -pattern current_observation -tab "\\n" \\
    -element weather temp_f wind_dir wind_mph

  nquire -url "https://eutils.ncbi.nlm.nih.gov/entrez/eutils" elink.fcgi \\
    -dbfrom pubmed -db pubmed -cmd neighbor -linkname pubmed_pubmed -id 2539356

  nquire -eutils efetch.fcgi -db pubmed -id 2539356 -rettype medline -retmode text

  nquire -eutils esummary.fcgi -db pubmed -id 2539356 -version 2.0

  nquire -eutils esearch.fcgi -db pubmed -term "tn3 transposition immunity" |
  xtract -pattern eSearchResult -element QueryTranslation

  nquire -bioc-pubmed 17299597 |
  xtract -pattern collection -block passage -if infon -equals title -element text

  nquire -bioc-pmc 1790863 |
  xtract -pattern collection -block passage -if infon -equals abstract -tab "\\n" -element text

  nquire -ftp ftp.ncbi.nlm.nih.gov pub/gdp ideogram_9606_GCF_000001305.14_850_V1 |
  grep acen | cut -f 1,2,6,7 | grep "^X\\t"

};

my $nquire_examples = qq{
Medical Subject Headings

  nquire -get "http://id.nlm.nih.gov/mesh/sparql" \\
    -query "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> \\
      SELECT DISTINCT ?class FROM <http://id.nlm.nih.gov/mesh> \\
      WHERE { ?s rdf:type ?class } ORDER BY ?class" |
  xtract -pattern result -pfx "meshv:" -first "uri[http://id.nlm.nih.gov/mesh/vocab#|]"

  meshv:AllowedDescriptorQualifierPair
  meshv:CheckTag
  meshv:Concept
  meshv:DisallowedDescriptorQualifierPair
  meshv:GeographicalDescriptor
  meshv:PublicationType
  meshv:Qualifier
  meshv:SCR_Chemical
  meshv:SCR_Disease
  meshv:SCR_Organism
  meshv:SCR_Protocol
  meshv:Term
  meshv:TopicalDescriptor
  meshv:TreeNumber

MeSH Predicates

  nquire -get "http://id.nlm.nih.gov/mesh/sparql" \\
    -query "SELECT DISTINCT ?p FROM <http://id.nlm.nih.gov/mesh> WHERE { ?s ?p ?o } ORDER BY ?p" |
  xtract -pattern result -pfx "meshv:" -first "uri[http://id.nlm.nih.gov/mesh/vocab#|]"

  meshv:abbreviation
  meshv:active
  meshv:allowableQualifier
  meshv:altLabel
  meshv:annotation
  meshv:broaderConcept
  meshv:broaderDescriptor
  meshv:broaderQualifier
  meshv:casn1_label
  meshv:concept
  meshv:considerAlso
  meshv:dateCreated
  meshv:dateEstablished
  meshv:dateRevised
  meshv:entryVersion
  meshv:frequency
  meshv:hasDescriptor
  meshv:hasQualifier
  meshv:historyNote
  meshv:identifier
  meshv:indexerConsiderAlso
  meshv:lastActiveYear
  meshv:lexicalTag
  meshv:mappedTo
  meshv:narrowerConcept
  meshv:nlmClassificationNumber
  meshv:note
  meshv:onlineNote
  meshv:parentTreeNumber
  meshv:pharmacologicalAction
  meshv:prefLabel
  meshv:preferredConcept
  meshv:preferredMappedTo
  meshv:preferredTerm
  meshv:previousIndexing
  meshv:publicMeSHNote
  meshv:registryNumber
  meshv:relatedConcept
  meshv:relatedRegistryNumber
  meshv:scopeNote
  meshv:seeAlso
  meshv:sortVersion
  meshv:source
  meshv:term
  meshv:thesaurusID
  meshv:treeNumber
  meshv:useInstead

WikiData Predicate List

  nquire -url "https://query.wikidata.org/sparql" \\
    -query "SELECT ?property ?propertyType ?propertyLabel \\
      ?propertyDescription ?propertyAltLabel WHERE { \\
      ?property wikibase:propertyType ?propertyType . SERVICE wikibase:label \\
      { bd:serviceParam wikibase:language '[AUTO_LANGUAGE],en'. } } \\
      ORDER BY ASC(xsd:integer(STRAFTER(STR(?property), 'P')))" |
  xtract -pattern result -first "uri[http://www.wikidata.org/entity/|]" -first literal

Selected WikiData Predicates

  P6       head of government
  P16      highway system
  P17      country
  P19      place of birth
  P21      sex or gender
  P22      father
  P25      mother
  P26      spouse
  P30      continent
  P31      instance of
  P35      head of state
  P36      capital
  P40      child
  P105     taxon rank
  P660     EC enzyme classification
  P672     MeSH Code
  P680     molecular function
  P681     cell component
  P682     biological process
  P685     NCBI Taxonomy ID
  P698     PubMed ID
  P699     Disease Ontology ID
  P932     PMCID
  P1340    eye color
  P2067    mass
  P2410    WikiPathways ID
  P2888    exact match

Vitamin Binding Site

  nquire -get "http://www.wikidata.org/entity/Q22679758" |
  transmute -j2x |
  xtract -pattern entities -group claims -block P527 -element "value\@id"

Children of JS Bach

  nquire -url "https://query.wikidata.org/sparql" \\
    -query "SELECT ?child ?childLabel WHERE \\
      { ?child wdt:P22 wd:Q1339. SERVICE wikibase:label \\
        { bd:serviceParam wikibase:language '[AUTO_LANGUAGE],en'. } }" |
  xtract -pattern result -block binding -if "\@name" -equals childLabel -element literal

Eye Color Frequency

  nquire -url "https://query.wikidata.org/sparql" \\
    -query "SELECT ?eyeColorLabel WHERE \\
      { ?human wdt:P31 wd:Q5. ?human wdt:P1340 ?eyeColor. SERVICE wikibase:label \\
        { bd:serviceParam wikibase:language '[AUTO_LANGUAGE],en'. } }" |
  xtract -pattern result -element literal |
  sort-uniq-count-rank

Federated Query

  nquire -url "https://query.wikidata.org/sparql" \\
    -query " \\
      PREFIX wp:      <http://vocabularies.wikipathways.org/wp#> \\
      PREFIX dcterms:  <http://purl.org/dc/terms/> \\
      PREFIX dc:      <http://purl.org/dc/elements/1.1/> \\
      SELECT DISTINCT ?metabolite1Label ?metabolite2Label ?mass1 ?mass2 WITH { \\
        SELECT ?metabolite1 ?metabolite2 WHERE { \\
          ?pathwayItem wdt:P2410 'WP706'; \\
                       wdt:P2888 ?pwIri. \\
          SERVICE <http://sparql.wikipathways.org/> { \\
            ?pathway dc:identifier ?pwIri. \\
            ?interaction rdf:type wp:Interaction; \\
                         wp:participants ?wpmb1, ?wpmb2; \\
                         dcterms:isPartOf ?pathway. \\
            FILTER (?wpmb1 != ?wpmb2) \\
            ?wpmb1 wp:bdbWikidata ?metabolite1. \\
            ?wpmb2 wp:bdbWikidata ?metabolite2. \\
          } \\
        } \\
      } AS %metabolites WHERE { \\
        INCLUDE %metabolites. \\
        ?metabolite1 wdt:P2067 ?mass1. \\
        ?metabolite2 wdt:P2067 ?mass2. \\
        SERVICE wikibase:label { bd:serviceParam wikibase:language '[AUTO_LANGUAGE],en'. } \\
      }" |
  xtract -pattern result -block binding -element "binding\@name" literal

BioThings Queries

  nquire -get -j2x -url "http://mygene.info/v3" query -q "pathway.wikipathways.id:WP455" -fetch_all true |
  xtract -pattern hits -element \@_id

  nquire -j2x -url "http://mygene.info/v3" query -q "WP455" -scopes "pathway.wikipathways.id" -size 300 |
  xtract -pattern anon -element \@_id

  nquire -variant variant "chr6:g.26093141G>A" -fields dbsnp.gene |
  xtract -pattern gene -element \@geneid

  nquire -gene query -q "symbol:OPN1MW" -species 9606 |
  xtract -pattern hits -element \@_id

  nquire -gene query -q "symbol:OPN1MW AND taxid:9606" |
  xtract -pattern hits -element \@_id

  nquire -gene gene 2652 -fields pathway.wikipathways |
  xtract -pattern pathway -element \@id

  nquire -gene query -q "pathway.wikipathways.id:WP455" -size 300 |
  xtract -pattern hits -element \@_id

  nquire -chem query -q "drugbank.targets.uniprot:P05231 AND drugbank.targets.actions:inhibitor" -fields hgvs |
  xtract -pattern hits -element \@_id

  nquire -get "http://myvariant.info/v1/variant/chr6:g.26093141G>A" \\
    -fields clinvar.rcv.conditions.identifiers \\
    -always_list clinvar.rcv.conditions.identifiers |
  xtract -j2x |
  xtract -biopath opt clinvar.rcv.conditions.identifiers.omim

  cat uniprots.xml |
  xtract -pattern ENTREZ_EXTEND -sep "\n" -element Id |
  while read uid
  do
    echo "UID \$uid"
    echo "\$uid" | xplore -load uniprot | xplore -link inchikey
    echo ""
  done

EDirect Expansion

  ExtractIDs() {
    xtract -pattern BIO_THINGS -block Id -tab "\\n" -element Id
  }

  WrapIDs() {
    xtract -wrp BIO_THINGS -pattern opt -wrp Type -lbl "\$1" \\
      -wrp Count -num "\$2" -block "\$2" -wrp Id -element "\$3" |
    xtract -format
  }

  nquire -gene query -q "symbol:OPN1MW AND taxid:9606" |
  WrapIDs entrezgene hits "\@entrezgene" |

  ExtractIDs |
  while read geneid
  do
    nquire -gene gene "\$geneid" -fields pathway.wikipathways
  done |
  WrapIDs pathway.wikipathways.id pathway "\@id" |

  ExtractIDs |
  while read pathid
  do
    nquire -gene query -q "pathway.wikipathways.id:\$pathid" -size 300
  done |
  WrapIDs entrezgene hits "\@entrezgene" |

  ExtractIDs |
  sort -n

};

my @pubchem_properties = qw(
  MolecularFormula
  MolecularWeight
  CanonicalSMILES
  IsomericSMILES
  InChI
  InChIKey
  IUPACName
  XLogP
  ExactMass
  MonoisotopicMass
  TPSA
  Complexity
  Charge
  HBondDonorCount
  HBondAcceptorCount
  RotatableBondCount
  HeavyAtomCount
  IsotopeAtomCount
  AtomStereoCount
  DefinedAtomStereoCount
  UndefinedAtomStereoCount
  BondStereoCount
  DefinedBondStereoCount
  UndefinedBondStereoCount
  CovalentUnitCount
  Volume3D
  XStericQuadrupole3D
  YStericQuadrupole3D
  ZStericQuadrupole3D
  FeatureCount3D
  FeatureAcceptorCount3D
  FeatureDonorCount3D
  FeatureAnionCount3D
  FeatureCationCount3D
  FeatureRingCount3D
  FeatureHydrophobeCount3D
  ConformerModelRMSD3D
  EffectiveRotorCount3D
  ConformerCount3D
  Fingerprint2D
);

sub nqir {

  %macros = ();
  $agent = "Nquire/1.0";
  $alias = "";
  $debug = false;
  $http = "";
  $j2x = false;
  $x2j = false;
  $output = "";

  # nquire -url http://... -tag value -tag value | ...

  $url = "";
  $arg = "";
  $pfx = "";
  $amp = "";
  $pat = "";
  $sfx = "";

  @args = @ARGV;
  $max = scalar @args;

  %biothingsHash = (
    '-gene'     =>  'http://mygene.info/v3',
    '-variant'  =>  'http://myvariant.info/v1',
    '-chem'     =>  'http://mychem.info/v1',
  );

  if ( $max < 1 ) {
    return;
  }

  if ( $ARGV[0] eq "-version" ) {
    print "$version\n";
    return;
  }

  if ( $ARGV[0] eq "-help" ) {
    print "nquire $version\n";
    print $nquire_help;
    return;
  }

  # -examples prints advanced sparql queries (undocumented)

  if ( $ARGV[0] eq "-examples" or $ARGV[0] eq "-example" or $ARGV[0] eq "-extras" or $ARGV[0] eq "-extra" ) {
    print "nquire $version\n";
    print $nquire_examples;
    return;
  }

  if ( $max < 2 ) {
    return;
  }

  $i = 0;

  # if present, -debug must be first argument, only prints generated URL (undocumented)

  if ( $i < $max ) {
    $pat = $args[$i];
    if ( $pat eq "-debug" ) {
      $i++;
      $debug = true;
    }
  }

  # if present, -delay must be next (if not using -eutils shortcut)

  if ( $i < $max ) {
    $pat = $args[$i];
    if ( $pat eq "-delay" ) {
      $i++;
      # add 1/3 second sleep to avoid server limit (undocumented)
      Time::HiRes::usleep(350000);
    }
  }

  # if present, -ftp must be next

  if ( $i < $max ) {
    $pat = $args[$i];
    if ( $pat eq "-ftp" ) {
      $i++;
      if ( $i < $max ) {
        my $server = $args[$i];
        if ( $server =~ /^ftp:\/\/(.+)/ ) {
          $server = $1;
        }
        $i++;
        if ( $i < $max ) {
          my $dir = $args[$i];
          $i++;
          if ( $i < $max ) {
            my $fl = $args[$i];

            my $ftp = new Net::FTP($server, Passive => 1)
              or die "Unable to connect to FTP server: $!";

            $ftp->login or die "Unable to log in to FTP server: ", $ftp->message;
            $ftp->cwd($dir) or die "Unable to change to $dir: ", $ftp->message;
            $ftp->binary or warn "Unable to set binary mode: ", $ftp->message;

            if (! $ftp->get($fl, "/dev/stdout") ) {
              my $msg = $ftp->message;
              chomp $msg;
              print STDERR "\nFAILED TO DOWNLOAD:\n\n$fl ($msg\n";
            }
          }
        }
      }
      return;
    }
  }

  # if present, -http get or -get must be next (now also allow -http post or -post)

  # nquire -get -url "http://collections.mnh.si.edu/services/resolver/resolver.php" -voucher "Birds:625456"

  if ( $i < $max ) {
    $pat = $args[$i];
    if ( $pat eq "-http" ) {
      $i++;
      if ( $i < $max ) {
        $http = $args[$i];
        $i++;
      }
    } elsif ( $pat eq "-get" ) {
      $http = "get";
      $pat = $args[$i + 1];
      # allow URL argument immediately after -get
      if ( $pat =~ /^-(.+)/ ) {
        $i++;
      }
    } elsif ( $pat eq "-post" ) {
      $http = "post";
      $pat = $args[$i + 1];
      # allow URL argument immediately after -post
      if ( $pat =~ /^-(.+)/ ) {
        $i++;
      }
    }
  }

  # if present, -agent must be next argument (undocumented)

  if ( $i < $max ) {
    $pat = $args[$i];
    if ( $pat eq "-agent" ) {
      $i++;
      if ( $i < $max ) {
        $agent = $args[$i];
        $i++;
      }
    }
  }

  # if present, -j2x or -x2j must be next argument (undocumented)

  if ( $i < $max ) {
    $pat = $args[$i];
    if ( $pat eq "-j2x" ) {
      $i++;
      $j2x = true;
    } elsif ( $pat eq "-x2j" ) {
      $i++;
      $x2j = true;
    }
  }

  # read file of keyword shortcuts for URL expansion

  if ( $i < $max ) {
    $pat = $args[$i];
    if ( $pat eq "-alias" ) {
      $i++;
      if ( $i < $max ) {
        $alias = $args[$i];
        if ( $alias ne "" ) {
          read_aliases ();
        }
        $i++;
      }
    }
  }

  # read URL

  # -get or -post can now be followed immediately by the URL, without a -url argument

  # nquire -get "http://collections.mnh.si.edu/services/resolver/resolver.php" -voucher "Birds:625456"

  if ( $i < $max ) {
    $pat = $args[$i];
    if ( $pat eq "-url" or $pat eq "-get" or $pat eq "-post" ) {
      $i++;
      if ( $i < $max ) {
        $url = $args[$i];
        $url = map_macros ($url);
        $i++;
      }
    } elsif ( $pat eq "-ncbi" ) {
      # shortcut for ncbi base (undocumented)
      $i++;
      if ( $i < $max ) {
        $url = "https://www.ncbi.nlm.nih.gov";
      }
    } elsif ( $pat eq "-eutils" ) {
      # shortcut for eutils base (undocumented)
      $i++;
      if ( $i < $max ) {
        $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";
        # add 1/3 second sleep to avoid server limit (undocumented)
        Time::HiRes::usleep(350000);
      }
    } elsif ( $pat eq "-test" ) {
      # shortcut for eutilstest base (undocumented)
      $i++;
      if ( $i < $max ) {
        $url = "https://eutilstest.ncbi.nlm.nih.gov/entrez/eutils";
      }
    } elsif ( $pat eq "-qa" ) {
      # shortcut for eutils QA base (undocumented)
      $i++;
      if ( $i < $max ) {
        $url = "http://qa.ncbi.nlm.nih.gov/entrez/eutils";
      }

    } elsif ( $pat eq "-hydra" ) {
      # internal citation match request (undocumented)
      $i++;
      if ( $i < $max ) {
        $url = "https://www.ncbi.nlm.nih.gov/projects/hydra/hydra_search.cgi";
        $pat = $args[$i];
        $pat = map_macros ($pat);
        $enc = do_uri_escape ($pat);
        $arg="search=pubmed_search_citation_top_20.1&query=$enc";
        $amp = "&";
        $i++;
      }

    } elsif ( $pat eq "-revhist" ) {
      # internal sequence revision history request (undocumented)
      $i++;
      if ( $i < $max ) {
        $url = "https://www.ncbi.nlm.nih.gov/sviewer/girevhist.cgi";
        $pat = $args[$i];
        $arg="cmd=seqid&txt=on&seqid=asntext&os=PUBSEQ_OS&val=$pat";
        $amp = "&";
        $i++;
      }

    } elsif ( $pat eq "-pubchem" ) {
      # shortcut for PubChem Power User Gateway REST service base (undocumented)
      # nquire -pubchem "compound/name/creatine/property" "IUPACName,MolecularWeight,MolecularFormula" "XML"
      $i++;
      if ( $i < $max ) {
        $url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug";
        if ( $i + 2 == $max && $args[$i] eq "compound" ) {
          # even shorter shortcut
          # nquire -pubchem compound creatine
          $pat = $args[$i + 1];
          if ( $pat =~ /^-(.+)/ ) {
          } elsif ( $pat !~ /\// ) {
            $i = $i + 2;
            $url .= "/compound/name/";
            $pat = map_macros ($pat);
            $url .= $pat;
            $url .= "/property/";
            $sfx = join(",", @pubchem_properties);
            $url .= $sfx;
            $url .= "/XML";
          }
        }
      }

    } elsif ( $pat eq "-bioc-pubmed" ) {
      # shortcut for BioC on PubMed (undocumented)
      $i++;
      if ( $i < $max ) {
        $id = $args[$i];
        $http = "get";
        $url = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pubmed.cgi/BioC_xml/$id/unicode";
        $amp = "&";
        $i++;
      }
    } elsif ( $pat eq "-bioc-pmc" ) {
      # shortcut for BioC on PMC (undocumented)
      $i++;
      if ( $i < $max ) {
        $id = $args[$i];
        if ( $id !~ /^PMC/ ) {
          # add PMC prefix if not already in argument
          $id = "PMC" . $id;
        }
        $http = "get";
        $url = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/$id/unicode";
        $amp = "&";
        $i++;
      }

    } elsif ( $pat eq "-biocxml-pubmed" ) {
      # shortcut for annotated PubMed BioC on PubTator (undocumented)
      $i++;
      if ( $i < $max ) {
        $id = $args[$i];
        $http = "get";
        $url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocxml?pmids=$id";
        $amp = "&";
        $i++;
      }
    } elsif ( $pat eq "-biocxml-pmc" ) {
      # shortcut for annotated PMC BioC on PubTator (undocumented)
      $i++;
      if ( $i < $max ) {
        $id = $args[$i];
        if ( $id !~ /^PMC/ ) {
          # add PMC prefix if not already in argument
          $id = "PMC" . $id;
        }
        $http = "get";
        $url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocxml?pmcids=$id";
        $amp = "&";
        $i++;
      }

    } elsif ( defined $biothingsHash{$pat} ) {
      # shortcuts for biothings services (undocumented)
      $i++;
      $url = $biothingsHash{$pat};
      if ( $http eq "" ) {
          $http = "get";
      }
      $j2x = true;

    } elsif ( $pat eq "-wikipathways" ) {
      # shortcut for webservice.wikipathways.org (undocumented)
      $i++;
      if ( $i < $max ) {
        $url = "http://webservice.wikipathways.org";
      }

    } elsif ( $pat eq "-biosample" ) {
      # internal biosample_chk request on live database (undocumented)
      $i++;
      if ( $i < $max ) {
        $http = "get";
        $url = "https://api-int.ncbi.nlm.nih.gov/biosample/fetch";
        $bid = $args[$i];
        $arg="format=source&id=$bid";
        $amp = "&";
        $i++;
      }
    } elsif ( $pat eq "-biosample-dev" ) {
      # internal biosample_chk request on development database (undocumented)
      $i++;
      if ( $i < $max ) {
        $http = "get";
        $url = "https://dev-api-int.ncbi.nlm.nih.gov/biosample/fetch";
        $bid = $args[$i];
        $arg="format=source&id=$bid";
        $amp = "&";
        $i++;
      }
    }
  }

  if ( $url eq "" ) {
    return;
  }

  # hard-coded URL aliases for common NCBI web sites

  if ( $url =~ /\(#/ ) {

    $ky = "ncbi_url";
    if ( $url =~ /\(#$ky\)/ ) {
      $vl = "https://www.ncbi.nlm.nih.gov";
      $url =~ s/\((\#$ky)\)/$vl/g;
    }

    $ky = "eutils_url";
    if ( $url =~ /\(#$ky\)/ ) {
      $vl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";
      $url =~ s/\((\#$ky)\)/$vl/g;
    }
  }

  # arguments before next minus are added to base URL as /value

  $go_on = true;
  while ( $i < $max and $go_on ) {
    $pat = $args[$i];
    if ( $pat =~ /^-(.+)/ ) {
      $go_on = false;
    } else {
      $pat = map_macros ($pat);
      $url .= "/" . $pat;
      $i++;
    }
  }

  # now expect tag with minus and value[s] without, add as &tag=value[,value]

  while ( $i < $max ) {
    $pat = $args[$i];
    if ( $pat =~ /^-(.+)/ ) {
      $pat = $1;
      $pfx = $amp . "$pat=";
      $amp = "";
    } else {
      $pat =~ s/^\\-/-/g;
      $pat = map_macros ($pat);
      $enc = do_uri_escape ($pat);
      $arg .= $pfx . $enc;
      $pfx = ",";
      $amp = "&";
    }
    $i++;
  }

  # perform query
  $output = do_nquire_post ($url, $arg);

  if ( $j2x ) {
    my $jc = JSON::PP->new->ascii->pretty->allow_nonref;
    my $conv = $jc->decode($output);
    convert_json($conv);
    my $result = XMLout($conv, SuppressEmpty => undef);

    # remove newlines, tabs, space between tokens, compress runs of spaces
    $result =~ s/\r/ /g;
    $result =~ s/\n/ /g;
    $result =~ s/\t//g;
    $result =~ s/ +/ /g;
    $result =~ s/> +</></g;

    # remove <opt> flanking object
    if ( $result =~ /<opt>\s*?</ and $result =~ />\s*?<\/opt>/ ) {
      $result =~ s/<opt>\s*?</</g;
      $result =~ s/>\s*?<\/opt>/>/g;
    }

    $output = "$result";

    # restore newlines between objects
    $output =~ s/> *?</>\n</g;

    binmode(STDOUT, ":utf8");
  }

  if ( $x2j ) {
    my $xc = new XML::Simple(KeepRoot => 1);
    my $conv = $xc->XMLin($output);
    convert_json($conv);
    my $jc = JSON::PP->new->ascii->pretty->allow_nonref;
    my $result = $jc->encode($conv);

    $output = "$result";
  }

  print "$output";
}

#  etest is an unadvertised function for development

sub etest {
#  $addr = get_email ();
#  print "e-mail:  $addr\n";
}

# main block dispatches control to appropriate subroutine

if ( scalar @ARGV > 0 and $ARGV[0] eq "-version" ) {
  print "$version\n";
} elsif ( $fnc eq "-search" ) {
  esrch ();
} elsif ( $fnc eq "-link" ) {
  elink ();
} elsif ( $fnc eq "-filter" ) {
  efilt ();
} elsif ( $fnc eq "-summary" ) {
  eftch ();
} elsif ( $fnc eq "-fetch" ) {
  eftch ();
} elsif ( $fnc eq "-info" ) {
  einfo ();
} elsif ( $fnc eq "-post" ) {
  epost ();
} elsif ( $fnc eq "-spell" ) {
  espel ();
} elsif ( $fnc eq "-citmatch" ) {
  ecitmtch ();
} elsif ( $fnc eq "-proxy" ) {
  eprxy ();
} elsif ( $fnc eq "-contact" ) {
  ecntc ();
} elsif ( $fnc eq "-notify" ) {
  entfy ();
} elsif ( $fnc eq "-address" ) {
  eaddr ();
} elsif ( $fnc eq "-blast" ) {
  eblst ();
} elsif ( $fnc eq "-ftpls" ) {
  ftls ();
} elsif ( $fnc eq "-aspls" ) {
  asls ();
} elsif ( $fnc eq "-ftpcp" ) {
  ftcp ();
} elsif ( $fnc eq "-tmute" ) {
  tmut ();
} elsif ( $fnc eq "-nquir" ) {
  nqir ();
} elsif ( $fnc eq "-test" ) {
  etest ();
} else {
  die "Function name '$fnc' not recognized\n";
}

# close input and output files

close (STDIN);
close (STDOUT);
close (STDERR);

exit $result;
