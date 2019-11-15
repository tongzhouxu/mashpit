#!/bin/sh

cd "$GOPATH"

go get -u github.com/fatih/color
go get -u github.com/fiam/gounidecode/unidecode
go get -u github.com/gedex/inflector
go get -u github.com/klauspost/cpuid
go get -u github.com/pbnjay/memory
go get -u github.com/surgebase/porter2
go get -u golang.org/x/text/runes
go get -u golang.org/x/text/transform
go get -u golang.org/x/text/unicode/norm

cd "$GOPATH/src/xtract"

go build -o xtract xtract.go common.go
go build -o rchive rchive.go common.go
