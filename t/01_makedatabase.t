#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests=>3;
use FindBin qw/$RealBin/;
use lib "$RealBin/../lib";
use Data::Dumper;

use_ok "Mashpit";

my $dbPath = "$RealBin/testDb";
my $sqlitePath = "$dbPath/db.sqlite";
my $mp = Mashpit->new($dbPath);

is(-d $dbPath, 1, "Database is a directory at $dbPath");
is(-f $sqlitePath, 1, "SQLite database found at $sqlitePath");

