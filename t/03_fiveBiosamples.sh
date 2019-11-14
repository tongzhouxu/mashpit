#!/usr/bin/env bats

here=$BATS_TEST_DIRNAME
db=$here/test.sqlite
script=$here/../scripts/addFromNcbi.sh

# Download five at once, test database locking, etc
@test "Add SAMN02182865" {
  $script $db SAMN02182865
}
@test "Add SAMN02182866" {
  $script $db SAMN02182866
}
@test "Add SAMN02182867" {
  $script $db SAMN02182867
}
@test "Add SAMN02182868" {
  $script $db SAMN02182868
}
@test "Add SAMN02182869" {
  $script $db SAMN02182869
}

