#!/usr/bin/env bats

here=$BATS_TEST_DIRNAME
db=$here/test.sqlite

# Test for six biosamples
@test "Number of biosamples" {
  run sqlite3 $db "SELECT COUNT(*) FROM BIOSAMPLE"
  [[ "${lines[0]}" -eq 6 ]]
}
# Number of sketches
@test "Number of sketches" {
  run sqlite3 $db "SELECT COUNT(*) FROM SKETCH"
  [[ "${lines[0]}" -eq 6 ]]
}
# Number of SRA entries
@test "Number of sketches" {
  run sqlite3 $db "SELECT COUNT(*) FROM SRA"
  [[ "${lines[0]}" -eq 6 ]]
}
# Two entries on taxonomy
@test "Taxonomy IDs" {
  run sqlite3 $db "SELECT taxid FROM TAXONOMY ORDER BY taxid ASC"
  [[ "${lines[0]}" -eq 28901 ]]
  [[ "${lines[1]}" -eq 59201 ]]
}
