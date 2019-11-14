#!/usr/bin/env bats

here=$BATS_TEST_DIRNAME
db=$here/test.sqlite
script=$here/../scripts/addFromNcbi.sh

# Just test one concurrently
@test "Add SAMN02182870" {
  $script $db SAMN02182870
}

