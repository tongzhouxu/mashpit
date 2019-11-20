#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw/$RealBin/;
use lib "$RealBin/../lib";
use Mashpit;
use Data::Dumper;

use Test::More tests=>3;

__END__

here=$BATS_TEST_DIRNAME
db=$here/test.sqlite
script=$here/../scripts/addFromNcbi.sh

# Just test one concurrently
@test "Add SAMN02182870" {
  $script $db SAMN02182870
}

