#!/bin/bash

# Run an sql statement with the right pragmas
function sql(){
  sqlite3 "pragma foreign_keys='on';
    $@
  "

