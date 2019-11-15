#!/usr/bin/env bats

here=$BATS_TEST_DIRNAME
db=$here/test.sqlite
script=$here/../scripts/addFromNcbi.sh
lockfile=$db.sketches/bats.lock

function rmlockfile(){
  if [ -e "$lockfile" ]; then
    rm -f $lockfile
    echo "# removed $lockfile" >&3
  fi
}
function setup(){
  rmlockfile
}
function teardown(){
  rmlockfile
}

# Download five at once, test database locking, etc
@test "Add SAMN02182865 concurrently" {
  lockfile $lockfile
  $script $db SAMN02182865
  rm -f $lockfile
}
@test "Add SAMN02182866 concurrently" {
  lockfile $lockfile
  $script $db SAMN02182866
  rm -f $lockfile
}
@test "Add SAMN02182867 concurrently" {
  lockfile $lockfile
  $script $db SAMN02182867
  rm -f $lockfile
}
@test "Add SAMN02182868 concurrently" {
  lockfile $lockfile
  $script $db SAMN02182868
  rm -f $lockfile
}
@test "Add SAMN02182869 concurrently" {
  lockfile $lockfile
  $script $db SAMN02182869
  rm -f $lockfile
}

