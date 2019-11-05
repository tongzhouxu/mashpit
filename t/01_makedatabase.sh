#!/usr/bin/env bats

here=$BATS_TEST_DIRNAME
db=$here/test.sqlite

@test "Create database" {
  run sqlite3 $db < $here/../scripts/mashpit.sql

  # Some filesize
  filesize=$(wc -c < $db)
  [[ "$filesize" -gt 1 ]]
  
  # Some tables to find
  dump=$(sqlite3 $db .dump)

  table_biosample=$(grep -o -i "CREATE TABLE BIOSAMPLE" <<< $dump || true)
  [[ "$table_biosample" == "CREATE TABLE BIOSAMPLE" ]]

  table_taxonomy=$(grep -o -i "CREATE TABLE TAXONOMY" <<< $dump || true)
  [[ "$table_taxonomy" == "CREATE TABLE TAXONOMY" ]]

  table_sketch=$(grep -o -i "CREATE TABLE SKETCH" <<< $dump || true)
  [[ "$table_sketch" == "CREATE TABLE SKETCH" ]]
}

