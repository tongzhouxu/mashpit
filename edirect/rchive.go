// ===========================================================================
//
//                            PUBLIC DOMAIN NOTICE
//            National Center for Biotechnology Information (NCBI)
//
//  This software/database is a "United States Government Work" under the
//  terms of the United States Copyright Act. It was written as part of
//  the author's official duties as a United States Government employee and
//  thus cannot be copyrighted. This software/database is freely available
//  to the public for use. The National Library of Medicine and the U.S.
//  Government do not place any restriction on its use or reproduction.
//  We would, however, appreciate having the NCBI and the author cited in
//  any work or product based on this material.
//
//  Although all reasonable efforts have been taken to ensure the accuracy
//  and reliability of the software and data, the NLM and the U.S.
//  Government do not and cannot warrant the performance or results that
//  may be obtained by using this software or data. The NLM and the U.S.
//  Government disclaim all warranties, express or implied, including
//  warranties of performance, merchantability or fitness for any particular
//  purpose.
//
// ===========================================================================
//
// File Name:  rchive.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package main

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"container/heap"
	"encoding/binary"
	"fmt"
	"github.com/fiam/gounidecode/unidecode"
	"github.com/klauspost/cpuid"
	"github.com/pbnjay/memory"
	"github.com/surgebase/porter2"
	"hash/crc32"
	"html"
	"io"
	"io/ioutil"
	"os"
	"os/user"
	"path"
	"regexp"
	"runtime"
	"runtime/debug"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
	"unicode"
)

// RCHIVE VERSION AND HELP MESSAGE TEXT

const rchiveVersion = "12.5"

const rchiveHelp = `
Processing Flags

  -strict     Remove HTML and MathML tags
  -mixed      Allow mixed content XML

Data Source

  -input      Read XML from file instead of stdin

Local Record Cache

  -archive    Base path for saving individual XML files
  -index      Use [parent/element@attribute^version] for identifier

  -fetch      Base path for retrieving XML files
  -stream     Path for retrieving compressed XML

  -flag       [strict|mixed|none]
  -gzip       Use compression for local XML files
  -hash       Print UIDs and checksum values to stdout

  -trie       Print archive trie

Local Record Index

  -e2index    Create Entrez index XML (in xtract)
  -invert     Generate inverted index
  -join       Collect subsets of inverted index files
  -fuse       Combine subsets of inverted index files
  -merge      Combine inverted indices, divide by term prefix
  -promote    Create term lists and posting files

  -path       Path to postings directory

  -query      Search on words or phrases in Boolean formulas
  -exact      Strict search for article title round-tripping

  -count      Print terms and counts, merging wildcards
  -counts     Expand wildcards, print individual term counts

Documentation

  -help       Print this document
  -version    Print version number

Sample File Download

  ftp-cp ftp.ncbi.nlm.nih.gov /entrez/entrezdirect/samples carotene.xml.zip
  unzip carotene.xml.zip
  rm carotene.xml.zip

Mammalian Sequence Download

  download-sequence gbmam gbpri gbrod

Human Subset Extraction

  #!/bin/sh

  for fl in gbpri?.aso.gz gbpri??.aso.gz
  do
    run-ncbi-converter asn2all -i "$fl" -a t -b -c -O 9606 -f s > ${fl%.aso.gz}.xml
  done

Populate PubMed Archive and Positional Term Index

  export EDIRECT_PUBMED_MASTER=/Volumes/cachet
  export EDIRECT_PUBMED_WORKING=/Volumes/scratch

  archive-pubmed

  index-pubmed

Retrieve from PubMed Archive

  cat subset.uid | fetch-pubmed > subset.xml

Entrez Indexing

  cat carotene.xml | xtract -strict -e2index > carotene.e2x

Index Inversion

  cat carotene.e2x | rchive -invert > carotene.inv

Merge Indices

  rchive -merge "$MASTER/Merged" carotene.inv

Create Postings

  rchive -promote "$MASTER/Postings" NORM carotene.mrg

Record Counts

  phrase-search -count "catabolite repress*"

Wildcard Expansion

  phrase-search -counts "catabolite repress*"

Query Processing

  phrase-search -query "selective serotonin reuptake inhibitor [STEM]"

  phrase-search -query "(literacy AND numeracy) NOT (adolescent OR child)"

  phrase-search -query "vitamin c + + common cold"

  phrase-search -query "vitamin c ~ ~ common cold"

  phrase-search -exact "Genetic Control of Biochemical Reactions in Neurospora."

Large-Scale Record Retrieval

  esearch -db pubmed -query "DNA Repair [MESH]" |
  efetch -format uid |
  fetch-pubmed |
  xtract -pattern PubmedArticle -num Author |
  sort-uniq-count -n |
  reorder-columns 2 1 |
  head -n 25 |
  tee /dev/tty |
  xy-plot auth.png

XML Data Transformation

  seconds_start=$(date "+%s")
  esearch -db pubmed -query "PNAS [JOUR]" -pub abstract |
  efetch -format uid | stream-pubmed | gunzip -c |
  xtract -stops -wrp Set,Rec -pattern PubmedArticle \
    -wrp "Year" -year "PubDate/*" \
    -wrp "Abst" -words Abstract/AbstractText |
  xtract -wrp Set,Pub -pattern Rec \
    -wrp "Year" -element Year \
    -wrp "Num" -num Abst > countsByYear.xml
  for yr in {1960..2020}
  do
    cat countsByYear.xml |
    xtract -wrp Raw -pattern Pub -select Year -eq "$yr" |
    xtract -pattern Raw -lbl "$yr" -avg Num
  done |
  tee /dev/tty |
  xy-plot verbosity.png
  rm countsByYear.xml
  seconds_end=$(date "+%s")
  seconds=$((seconds_end - seconds_start))
  echo "$seconds seconds"

Query Automation

  ascend_mesh_tree() {
    var="${1%\*}"
    while :
    do
      phrase-search -count "$var* [TREE]"
      case "$var" in
        *.* ) var="${var%????}" ;;
        *   ) break             ;;
      esac
    done
  }

  ascend_mesh_tree "C14.907.617.812"

  6584       c14 907 617 812*
  50722      c14 907 617*
  1567114    c14 907*
  2232414    c14*

Medical Subject Heading Code Viewer

  https://meshb.nlm.nih.gov/treeView

DISABLE ANTI-VIRUS FILE SCANNING FOR LOCAL ARCHIVES OR DESIGNATE AS TRUSTED FILES

DISABLE SPOTLIGHT INDEXING FOR EXTERNAL DISKS CONTAINING LOCAL ARCHIVES
`

const rchiveExtras = `
Maintenance Commands

  -prepare    [release|report] Compare daily update to archive
  -ignore     Ignore contents of object in -prepare comparisons
  -damaged    Report UIDs containing damaged embedded HTML tags
  -missing    Print list of missing identifiers
  -unique     File of UIDs for skipping all but last version

Miscellaneous

  -head       Print before everything else
  -tail       Print after everything else
  -hd         Print before each record
  -tl         Print after each record

Update Candidate Report

  cd "$MASTER/Pubmed"
  gunzip -c *.xml.gz | xtract -strict -compress -format flush |
  rchive -prepare report -ignore DateRevised -archive "$MASTER/Archive" \
    -index MedlineCitation/PMID -pattern PubmedArticle

Unnecessary Update Removal

  cd "$MASTER/Pubmed"
  gunzip -c *.xml.gz | xtract -strict -compress -format flush |
  rchive -prepare release -ignore DateRevised -archive "$MASTER/Archive" -index MedlineCitation/PMID \
    -head "<PubmedArticleSet>" -tail "</PubmedArticleSet>" -pattern PubmedArticle |
  xtract -format indent -xml '<?xml version="1.0" encoding="UTF-8"?>' \
    -doctype '<!DOCTYPE PubmedArticleSet SYSTEM "http://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_180101.dtd">' |
  gzip > newupdate.xml.gz

Get Archive UID List

  pm-uids "$MASTER/Archive" > complete.uid

Reconstruct List of Versioned PMIDs

  cd "$MASTER/Pubmed"
  rm -f "$MASTER/Archive/versioned.uid"
  gunzip -c *.xml.gz |
  xtract -strict -pattern PubmedArticle -if MedlineCitation/PMID@Version -gt 1 \
    -element MedlineCitation/PMID > "$MASTER/Archive/versioned.uid"

Reconstruct Release Files

  split -a 3 -l 30000 release.uid uids-
  n=1
  for x in uids-???
  do
    xmlfile=$(printf "pubmed18n%04d.xml.gz" "$n")
    n=$((n+1))
    echo "$xmlfile"
    cat "$x" |
    rchive -fetch "$MASTER/Archive" -head "<PubmedArticleSet>" -tail "</PubmedArticleSet>" |
    xtract -strict -format indent -xml '<?xml version="1.0" encoding="UTF-8"?>' \
      -doctype '<!DOCTYPE PubmedArticleSet SYSTEM "http://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_180101.dtd">' |
    gzip > "$xmlfile"
  done
  rm -rf uids-???

Damaged Embedded HTML Tag Search

  for fl in *.xml.gz
  do
    echo "$fl"
    gunzip -c "$fl" | rchive -mixed -damaged -index MedlineCitation/PMID^Version -pattern PubmedArticle
  done

  grep -v pubmed18n | grep AMPER | cut -f 1,6

Reconstruct Term List Keys

  rm -f "$MASTER/Postings/sections.txt"
  find "$MASTER/Postings" -name "*.mst" |
  sed -e 's,.*/\(.*\)\.mst,\1,' |
  sort | uniq > "$MASTER/Postings/sections.txt"

Generate Term List Paths

  find "$MASTER/Postings" -name "*.trm" |
  sed -e 's,\(.*/\)\(.*\.trm\),\1 \2,' |
  sort -k 2 | uniq | tr -d ' '
`

const rchiveInternal = `
Performance Default Overrides

  -proc     Number of CPU processors used
  -cons     Ratio of parsers to processors
  -serv     Concurrent parser instances
  -chan     Communication channel depth
  -heap     Order restoration heap size
  -farm     Node allocation buffer length
  -gogc     Garbage collection tuning knob

Debugging

  -debug    Display run-time parameter summary
  -stats    Print performance tuning values
  -timer    Report processing duration and rate

Entrez Invert Performance Measurement

  InvertTrials() {
    echo -e "<Trials>"
    for tries in {1..5}
    do
      cat "$1" | rchive -debug -proc "$2" -invert
    done
    echo -e "</Trials>"
  }

  for proc in {1..8}
  do
    InvertTrials "carotene.e2x" "$proc" |
    xtract -pattern Trials -lbl "$proc" -avg Rate -dev Rate
  done

Execution Profiling

  cat carotene.e2x | ./rchive -profile -invert > /dev/null
  go tool pprof --pdf ./rchive ./cpu.pprof > ./callgraph.pdf
`

// GLOBAL VARIABLES

// DATA OBJECTS

type Master struct {
	TermOffset int32
	PostOffset int32
}

type Arrays struct {
	Data []int32
	Ofst [][]int16
	Dist int
}

// UTILITIES

func ReportEncodedMarkup(typ, id, str string) {

	var buffer strings.Builder

	max := len(str)

	lookAhead := func(txt string, to int) string {

		mx := len(txt)
		if to > mx {
			to = mx
		}
		pos := strings.Index(txt[:to], "gt;")
		if pos > 0 {
			to = pos + 3
		}
		return txt[:to]
	}

	findContext := func(fr, to int) string {

		numSpaces := 0

		for fr > 0 {
			ch := str[fr]
			if ch == ' ' {
				numSpaces++
				if numSpaces > 1 {
					fr++
					break
				}
			} else if ch == '\n' || ch == '>' {
				fr++
				break
			}
			fr--
		}

		numSpaces = 0

		for to < max {
			ch := str[to]
			if ch == ' ' {
				numSpaces++
				if numSpaces > 1 {
					break
				}
			} else if ch == '\n' || ch == '<' {
				break
			}
			to++
		}

		return str[fr:to]
	}

	reportMarkup := func(lbl string, fr, to int, txt string) {

		if lbl == typ || typ == "ALL" {
			// extract XML of SELF, SINGLE, DOUBLE, or AMPER types, or ALL
			buffer.WriteString(str)
			buffer.WriteString("\n")
		} else if typ == "" {
			// print report
			buffer.WriteString(id)
			buffer.WriteString("\t")
			buffer.WriteString(lbl)
			buffer.WriteString("\t")
			buffer.WriteString(txt)
			buffer.WriteString("\t| ")
			ctx := findContext(fr, to)
			buffer.WriteString(ctx)
			if HasUnicodeMarkup(ctx) {
				ctx = RepairUnicodeMarkup(ctx, SPACE)
			}
			ctx = RepairEncodedMarkup(ctx)
			buffer.WriteString("\t| ")
			buffer.WriteString(ctx)
			if HasAmpOrNotASCII(ctx) {
				ctx = html.UnescapeString(ctx)
			}
			buffer.WriteString("\t| ")
			buffer.WriteString(ctx)
			buffer.WriteString("\n")
		}
	}

	/*
		badTags := [10]string{
			"<i/>",
			"<i />",
			"<b/>",
			"<b />",
			"<u/>",
			"<u />",
			"<sup/>",
			"<sup />",
			"<sub/>",
			"<sub />",
		}
	*/

	skip := 0

	/*
		var prev rune
	*/

	for i, ch := range str {
		if skip > 0 {
			skip--
			continue
		}
		/*
			if ch > 127 {
				if IsUnicodeSuper(ch) {
					if IsUnicodeSubsc(prev) {
						// reportMarkup("UNIUP", i, i+2, string(ch))
					}
				} else if IsUnicodeSubsc(ch) {
					if IsUnicodeSuper(prev) {
						// reportMarkup("UNIDN", i, i+2, string(ch))
					}
				} else if ch == '\u0038' || ch == '\u0039' {
					// reportMarkup("ANGLE", i, i+2, string(ch))
				}
				prev = ch
				continue
			} else {
				prev = ' '
			}
		*/
		if ch == '<' {
			/*
				j := i + 1
				if j < max {
					nxt := str[j]
					if nxt == 'i' || nxt == 'b' || nxt == 'u' || nxt == 's' {
						for _, tag := range badTags {
							if strings.HasPrefix(str, tag) {
								k := len(tag)
								reportMarkup("SELF", i, i+k, tag)
								break
							}
						}
					}
				}
				if strings.HasPrefix(str[i:], "</sup><sub>") {
					// reportMarkup("SUPSUB", i, i+11, "</sup><sub>")
				} else if strings.HasPrefix(str[i:], "</sub><sup>") {
					// reportMarkup("SUBSUP", i, i+11, "</sub><sup>")
				}
			*/
			continue
		} else if ch != '&' {
			continue
		} else if strings.HasPrefix(str[i:], "&lt;") {
			sub := lookAhead(str[i:], 14)
			_, ok := htmlRepair[sub]
			if ok {
				skip = len(sub) - 1
				reportMarkup("SINGLE", i, i+skip+1, sub)
				continue
			}
		} else if strings.HasPrefix(str[i:], "&amp;lt;") {
			sub := lookAhead(str[i:], 22)
			_, ok := htmlRepair[sub]
			if ok {
				skip = len(sub) - 1
				reportMarkup("DOUBLE", i, i+skip+1, sub)
				continue
			}
		} else if strings.HasPrefix(str[i:], "&amp;amp;") {
			reportMarkup("AMPER", i, i+9, "&amp;amp;")
			skip = 8
			continue
		}
	}

	res := buffer.String()

	os.Stdout.WriteString(res)
}

// DIRECTORY PATH UTILITIES

// MakeArchiveTrie allows a short prefix of letters with an optional underscore, and splits the remainder into character pairs
func MakeArchiveTrie(str string, arry [132]rune) string {

	if len(str) > 64 {
		return ""
	}

	if IsAllDigits(str) {

		// pad numeric identifier to 8 characters with leading zeros
		ln := len(str)
		if ln < 8 {
			zeros := "00000000"
			str = zeros[ln:] + str
		}
	}

	if IsAllDigitsOrPeriod(str) {

		// limit trie to first 6 characters
		if len(str) > 6 {
			str = str[:6]
		}
	}

	max := 4
	k := 0
	for _, ch := range str {
		if unicode.IsLetter(ch) {
			k++
			continue
		}
		if ch == '_' {
			k++
			max = 6
		}
		break
	}

	// prefix is up to three letters if followed by digits, or up to four letters if followed by an underscore
	pfx := str[:k]
	if len(pfx) < max {
		str = str[k:]
	} else {
		pfx = ""
	}

	i := 0

	if pfx != "" {
		for _, ch := range pfx {
			arry[i] = ch
			i++
		}
		arry[i] = '/'
		i++
	}

	between := 0
	doSlash := false

	// remainder is divided in character pairs, e.g., NP_/06/00/51 for NP_060051.2
	for _, ch := range str {
		// break at period separating accession from version
		if ch == '.' {
			break
		}
		if doSlash {
			arry[i] = '/'
			i++
			doSlash = false
		}
		if ch == ' ' {
			ch = '_'
		}
		if !unicode.IsLetter(ch) && !unicode.IsDigit(ch) {
			ch = '_'
		}
		arry[i] = ch
		i++
		between++
		if between > 1 {
			doSlash = true
			between = 0
		}
	}

	res := string(arry[:i])

	if !strings.HasSuffix(res, "/") {
		arry[i] = '/'
		i++
		res = string(arry[:i])
	}

	return strings.ToUpper(res)
}

// MakePostingsTrie splits a string into characters, separated by path delimiting slashes
func MakePostingsTrie(str string, arry [516]rune) string {

	if len(str) > 256 {
		return ""
	}

	// expand Greek letters, anglicize characters in other alphabets
	if IsNotASCII(str) {
		if HasGreek(str) {
			str = SpellGreek(str)
			str = CompressRunsOfSpaces(str)
		}
		str = unidecode.Unidecode(str)
		str = strings.TrimSpace(str)
	}

	i := 0
	doSlash := false

	for _, ch := range str {
		if doSlash {
			arry[i] = '/'
			i++
		}
		if ch == ' ' {
			ch = '_'
		}
		if !unicode.IsLetter(ch) && !unicode.IsDigit(ch) {
			ch = '_'
		}
		arry[i] = ch
		i++
		doSlash = true
	}

	return strings.ToLower(string(arry[:i]))
}

// POSTINGS FILE UTILITIES

// trieLen directory depth parameters are based on the observed size distribution of PubMed indices
var trieLen = map[string]int{
	"19": 4,
	"20": 4,
	"ac": 4,
	"af": 4,
	"an": 4,
	"c0": 3,
	"ca": 4,
	"ce": 4,
	"cl": 4,
	"co": 4,
	"d0": 3,
	"de": 3,
	"di": 4,
	"ex": 4,
	"ge": 4,
	"gr": 4,
	"he": 4,
	"hi": 4,
	"in": 4,
	"me": 4,
	"mo": 4,
	"no": 4,
	"pa": 4,
	"pe": 4,
	"pl": 4,
	"po": 4,
	"pr": 4,
	"re": 4,
	"si": 4,
	"sp": 4,
	"st": 4,
	"su": 4,
	"tr": 4,
	"tw": 4,
	"un": 3,
	"va": 3,
	"ve": 3,
	"vi": 3,
	"wh": 3,
}

func PostingDir(term string) string {

	if len(term) < 3 {
		return term
	}

	key := term[:2]

	num, ok := trieLen[key]
	if ok && len(term) >= num {
		return term[:num]
	}

	switch term[0] {
	case 'u', 'v', 'w', 'x', 'y', 'z':
		return term[:2]
	}

	return term[:3]
}

func IdentifierKey(term string) string {

	// remove punctuation from term
	key := strings.Map(func(c rune) rune {
		if !unicode.IsLetter(c) && !unicode.IsDigit(c) && c != ' ' && c != '-' && c != '_' {
			return -1
		}
		return c
	}, term)

	key = strings.Replace(key, " ", "_", -1)
	key = strings.Replace(key, "-", "_", -1)

	// use first 2, 3, or 4 characters of identifier for directory
	key = PostingDir(key)

	return key
}

func PostingPath(prom, field, term string, arry [516]rune) (string, string) {

	// use first few characters of identifier for directory
	dir := IdentifierKey(term)

	trie := MakePostingsTrie(dir, arry)
	if trie == "" {
		return "", ""
	}

	dpath := path.Join(prom, field, trie)

	return dpath, dir
}

func CommonOpenFile(dpath, fname string) (*os.File, int64) {

	fpath := path.Join(dpath, fname)
	if fpath == "" {
		return nil, 0
	}

	inFile, err := os.Open(fpath)
	if err != nil && os.IsNotExist(err) {
		return nil, 0
	}
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil, 0
	}

	fi, err := inFile.Stat()
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil, 0
	}

	size := fi.Size()

	return inFile, size
}

func ReadMasterIndex(dpath, key, field string) []Master {

	inFile, size := CommonOpenFile(dpath, key+"."+field+".mst")
	if inFile == nil {
		return nil
	}

	defer inFile.Close()

	data := make([]Master, size/8)
	if data == nil || len(data) < 1 {
		return nil
	}

	err := binary.Read(inFile, binary.LittleEndian, &data)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	return data
}

func ReadTermList(dpath, key, field string) []byte {

	inFile, size := CommonOpenFile(dpath, key+"."+field+".trm")
	if inFile == nil {
		return nil
	}

	defer inFile.Close()

	data := make([]byte, size)
	if data == nil || len(data) < 1 {
		return nil
	}

	err := binary.Read(inFile, binary.LittleEndian, &data)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	return data
}

func ReadPostingData(dpath, key, field string, offset int32, size int32) []int32 {

	inFile, _ := CommonOpenFile(dpath, key+"."+field+".pst")
	if inFile == nil {
		return nil
	}

	defer inFile.Close()

	data := make([]int32, size/4)
	if data == nil || len(data) < 1 {
		return nil
	}

	_, err := inFile.Seek(int64(offset), io.SeekStart)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	err = binary.Read(inFile, binary.LittleEndian, data)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	return data
}

func ReadPositionIndex(dpath, key, field string, offset int32, size int32) []int32 {

	inFile, _ := CommonOpenFile(dpath, key+"."+field+".uqi")
	if inFile == nil {
		return nil
	}

	defer inFile.Close()

	data := make([]int32, size/4)
	if data == nil || len(data) < 1 {
		return nil
	}

	_, err := inFile.Seek(int64(offset), io.SeekStart)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	err = binary.Read(inFile, binary.LittleEndian, data)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	return data
}

func ReadOffsetData(dpath, key, field string, offset int32, size int32) []int16 {

	inFile, _ := CommonOpenFile(dpath, key+"."+field+".ofs")
	if inFile == nil {
		return nil
	}

	defer inFile.Close()

	data := make([]int16, size/2)
	if data == nil || len(data) < 1 {
		return nil
	}

	_, err := inFile.Seek(int64(offset), io.SeekStart)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	err = binary.Read(inFile, binary.LittleEndian, data)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	return data
}

func ReadMasterIndexFuture(dpath, key, field string) <-chan []Master {

	out := make(chan []Master, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create master index channel\n")
		os.Exit(1)
	}

	// masterIndexFuture asynchronously gets the master file and sends results through channel
	masterIndexFuture := func(dpath, key, field string, out chan<- []Master) {

		data := ReadMasterIndex(dpath, key, field)

		out <- data

		close(out)
	}

	// launch single future goroutine
	go masterIndexFuture(dpath, key, field, out)

	return out
}

func ReadTermListFuture(dpath, key, field string) <-chan []byte {

	out := make(chan []byte, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create term list channel\n")
		os.Exit(1)
	}

	// termListFuture asynchronously gets posting IDs and sends results through channel
	termListFuture := func(dpath, key, field string, out chan<- []byte) {

		data := ReadTermList(dpath, key, field)

		out <- data

		close(out)
	}

	// launch single future goroutine
	go termListFuture(dpath, key, field, out)

	return out
}

func GetPostingIDs(prom, term, field string, simple bool) ([]int32, [][]int16) {

	var (
		arry [516]rune
	)

	dpath, key := PostingPath(prom, field, term, arry)
	if dpath == "" {
		return nil, nil
	}

	// schedule asynchronous fetching
	mi := ReadMasterIndexFuture(dpath, key, field)

	tl := ReadTermListFuture(dpath, key, field)

	// fetch master index and term list
	indx := <-mi

	trms := <-tl

	if indx == nil || len(indx) < 1 {
		return nil, nil
	}

	if trms == nil || len(trms) < 1 {
		return nil, nil
	}

	// master index is padded with phantom term and postings position
	numTerms := len(indx) - 1

	strs := make([]string, numTerms)
	if strs == nil || len(strs) < 1 {
		return nil, nil
	}

	retlength := int32(len("\n"))

	// populate array of strings from term list
	for i, j := 0, 1; i < numTerms; i++ {
		from := indx[i].TermOffset
		to := indx[j].TermOffset - retlength
		j++
		txt := string(trms[from:to])
		strs[i] = txt
	}

	// change protecting underscore to space
	term = strings.Replace(term, "_", " ", -1)

	// if term ends with dollar sign, use porter2 stemming, then add asterisk
	if strings.HasSuffix(term, "$") && term != "$" {
		term = strings.TrimSuffix(term, "$")
		term = porter2.Stem(term)
		term += "*"
	}

	isWildCard := false
	if strings.HasSuffix(term, "*") && term != "*" {
		tlen := len(term)
		isWildCard = true
		term = strings.TrimSuffix(term, "*")
		pdlen := len(PostingDir(term))
		if tlen < pdlen {
			fmt.Fprintf(os.Stderr, "Wildcard term '%s' must be at least %d characters long - ignoring this word\n", term, pdlen)
			return nil, nil
		}
	}

	// binary search in term list
	L, R := 0, numTerms-1
	for L < R {
		mid := (L + R) / 2
		if strs[mid] < term {
			L = mid + 1
		} else {
			R = mid
		}
	}

	// wild card search scans term lists, fuses adjacent postings lists
	if isWildCard {
		if R < numTerms && strings.HasPrefix(strs[R], term) {
			offset := indx[R].PostOffset
			for R < numTerms && strings.HasPrefix(strs[R], term) {
				R++
			}
			size := indx[R].PostOffset - offset

			// read relevant postings list section
			data := ReadPostingData(dpath, key, field, offset, size)
			if data == nil || len(data) < 1 {
				return nil, nil
			}

			if simple {

				merged := make(map[int32]bool)

				// combine all postings in term range
				for _, val := range data {
					merged[val] = true
				}

				fused := make([]int32, len(merged))

				// convert map to slice
				i := 0
				for num := range merged {
					fused[i] = num
					i++
				}

				sort.Slice(fused, func(i, j int) bool { return fused[i] < fused[j] })

				return fused, nil
			}

			// read relevant word position section, includes phantom offset at end
			uqis := ReadPositionIndex(dpath, key, field, offset, size+4)
			if uqis == nil {
				return nil, nil
			}
			ulen := len(uqis)
			if ulen < 1 {
				return nil, nil
			}

			from := uqis[0]
			to := uqis[ulen-1]

			// read offset section
			ofst := ReadOffsetData(dpath, key, field, from, to-from)
			if ofst == nil {
				return nil, nil
			}

			combo := make(map[int32][]int16)

			addPositions := func(uid int32, pos int16) {

				arrs, ok := combo[uid]
				if !ok {
					arrs = make([]int16, 0, 1)
				}
				arrs = append(arrs, pos)
				combo[uid] = arrs
			}

			// populate array of positions per UID
			for i, j, k := 0, 1, int32(0); i < ulen-1; i++ {
				uid := data[i]
				num := (uqis[j] - uqis[i]) / 2
				j++
				for q := k; q < k+num; q++ {
					addPositions(uid, ofst[q])
				}
				k += num
			}

			fused := make([]int32, len(combo))

			// convert map to slice
			i := 0
			for num := range combo {
				fused[i] = num
				i++
			}

			sort.Slice(fused, func(i, j int) bool { return fused[i] < fused[j] })

			// make array of int16 arrays, populate for each UID
			arrs := make([][]int16, ulen-1)
			if arrs == nil {
				return nil, nil
			}

			for j, uid := range fused {
				posn := combo[uid]

				if len(posn) > 1 {
					sort.Slice(posn, func(i, j int) bool { return posn[i] < posn[j] })
				}

				arrs[j] = posn
			}

			return fused, arrs
		}

		return nil, nil
	}

	// regular search requires exact match from binary search
	if R < numTerms && strs[R] == term {

		offset := indx[R].PostOffset
		size := indx[R+1].PostOffset - offset

		// read relevant postings list section
		data := ReadPostingData(dpath, key, field, offset, size)
		if data == nil || len(data) < 1 {
			return nil, nil
		}

		if simple {
			return data, nil
		}

		// read relevant word position section, includes phantom offset at end
		uqis := ReadPositionIndex(dpath, key, field, offset, size+4)
		if uqis == nil {
			return nil, nil
		}
		ulen := len(uqis)
		if ulen < 1 {
			return nil, nil
		}

		from := uqis[0]
		to := uqis[ulen-1]

		// read offset section
		ofst := ReadOffsetData(dpath, key, field, from, to-from)
		if ofst == nil {
			return nil, nil
		}

		// make array of int16 arrays, populate for each UID
		arrs := make([][]int16, ulen)
		if arrs == nil || len(arrs) < 1 {
			return nil, nil
		}

		// populate array of positions per UID
		for i, j, k := 0, 1, int32(0); i < ulen-1; i++ {
			num := (uqis[j] - uqis[i]) / 2
			j++
			arrs[i] = ofst[k : k+num]
			k += num
		}

		return data, arrs
	}

	return nil, nil
}

func PrintTermCount(base, term, field string) int {

	data, _ := GetPostingIDs(base, term, field, true)
	size := len(data)
	fmt.Fprintf(os.Stdout, "%d\t%s\n", size, term)

	return size
}

func PrintTermCounts(base, term, field string) int {

	pdlen := len(PostingDir(term))

	if len(term) < pdlen {
		fmt.Fprintf(os.Stderr, "\nERROR: Term count argument must be at least %d characters\n", pdlen)
		os.Exit(1)
	}

	if strings.Contains(term[:pdlen], "*") {
		fmt.Fprintf(os.Stderr, "\nERROR: Wildcard asterisk must not be in first %d characters\n", pdlen)
		os.Exit(1)
	}

	var arry [516]rune
	dpath, key := PostingPath(base, field, term, arry)
	if dpath == "" {
		return 0
	}

	// schedule asynchronous fetching
	mi := ReadMasterIndexFuture(dpath, key, field)

	tl := ReadTermListFuture(dpath, key, field)

	// fetch master index and term list
	indx := <-mi

	trms := <-tl

	if indx == nil || len(indx) < 1 {
		return 0
	}

	if trms == nil || len(trms) < 1 {
		return 0
	}

	// master index is padded with phantom term and postings position
	numTerms := len(indx) - 1

	strs := make([]string, numTerms)
	if strs == nil || len(strs) < 1 {
		return 0
	}

	retlength := int32(len("\n"))

	// populate array of strings from term list
	for i, j := 0, 1; i < numTerms; i++ {
		from := indx[i].TermOffset
		to := indx[j].TermOffset - retlength
		j++
		txt := string(trms[from:to])
		strs[i] = txt
	}

	// change protecting underscore to space
	term = strings.Replace(term, "_", " ", -1)

	// flank pattern with start-of-string and end-of-string symbols
	pat := "^" + term + "$"

	// change asterisk in query to dot + star for regular expression
	pat = strings.Replace(pat, "*", ".*", -1)

	re, err := regexp.Compile(pat)

	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return 0
	}

	count := 0

	for R, str := range strs {
		if re.MatchString(str) {
			offset := indx[R].PostOffset
			size := indx[R+1].PostOffset - offset
			fmt.Fprintf(os.Stdout, "%d\t%s\n", size/4, str)
			count++
		}
	}

	return count
}

func PrintTermPositions(base, term, field string) int {

	data, ofst := GetPostingIDs(base, term, field, false)
	size := len(data)
	fmt.Fprintf(os.Stdout, "\n%d\t%s\n\n", size, term)

	for i := 0; i < len(data); i++ {
		fmt.Fprintf(os.Stdout, "%12d\t", data[i])
		pos := ofst[i]
		sep := ""
		for j := 0; j < len(pos); j++ {
			fmt.Fprintf(os.Stdout, "%s%d", sep, pos[j])
			sep = ","
		}
		fmt.Fprintf(os.Stdout, "\n")
	}

	return size
}

// BOOLEAN OPERATIONS FOR POSTINGS LISTS

func ExtendPositionalIDs(N []int32, np [][]int16, M []int32, mp [][]int16, delta int, proc func(pn, pm []int16, dlt int16) []int16) ([]int32, [][]int16) {

	if proc == nil {
		return nil, nil
	}

	if N == nil || len(N) < 1 || np == nil || len(np) < 1 {
		return M, mp
	}
	if M == nil || len(M) < 1 || mp == nil || len(mp) < 1 {
		return N, np
	}

	n, m := len(N), len(M)

	// order matters when extending phrase or testing proximity, do not swap lists based on size

	sz := n
	if sz > m {
		sz = m
	}

	if sz < 1 {
		return N, np
	}

	res := make([]int32, sz)
	ofs := make([][]int16, sz)

	if res == nil || len(res) < 1 || ofs == nil || len(ofs) < 1 {
		return nil, nil
	}

	i, j, k := 0, 0, 0

	// use local variables for speed
	en, em := N[i], M[j]

	for {
		// do inequality tests first
		if en < em {
			i++
			if i == n {
				break
			}
			en = N[i]
		} else if en > em {
			j++
			if j == m {
				break
			}
			em = M[j]
		} else {
			// specific callbacks test position arrays to match terms by adjacency or phrases by proximity
			adj := proc(np[i], mp[j], int16(delta))
			if adj != nil && len(adj) > 0 {
				res[k] = en
				ofs[k] = adj
				k++
			}
			i++
			j++
			if i == n || j == m {
				break
			}
			en = N[i]
			em = M[j]
		}
	}

	// truncate output arrays to actual size of intersection
	res = res[:k]
	ofs = ofs[:k]

	return res, ofs
}

func IntersectIDs(N, M []int32) []int32 {

	if N == nil {
		return M
	}
	if M == nil {
		return N
	}

	n, m := len(N), len(M)

	// swap to make M the smaller list
	if n < m {
		N, M = M, N
		n, m = m, n
	}

	if m < 1 {
		return N
	}

	res := make([]int32, m)

	i, j, k := 0, 0, 0

	// use local variables for speed
	en, em := N[i], M[j]

	for {
		// do inequality tests first
		if en < em {
			// index to larger list most likely to be advanced
			i++
			if i == n {
				break
			}
			en = N[i]
		} else if en > em {
			j++
			if j == m {
				break
			}
			em = M[j]
		} else {
			// equality (intersection match) least likely
			res[k] = en
			k++
			i++
			j++
			if i == n || j == m {
				break
			}
			en = N[i]
			em = M[j]
		}
	}

	// truncate output array to actual size of intersection
	res = res[:k]

	return res
}

// if m * log(n) < m + n, binary search has fewer comparisons, but processor memory caches make linear algorithm faster
/*
func IntersectBinary(N, M []int32) []int32 {

	if N == nil {
		return M
	}
	if M == nil {
		return N
	}

	n, m := len(N), len(M)

	// swap to make M the smaller list
	if n < m {
		N, M = M, N
		n, m = m, n
	}

	if m < 1 {
		return N
	}

	k := 0

	res := make([]int32, m)

	for _, uid := range M {
		// inline binary search is faster than sort.Search
		L, R := 0, n-1
		for L < R {
			mid := (L + R) / 2
			if N[mid] < uid {
				L = mid + 1
			} else {
				R = mid
			}
		}
		// R := sort.Search(len(N), func(i int) bool { return N[i] >= uid })
		if R < n && N[R] == uid {
			res[k] = uid
			k++
			// remove leading part of N for slight speed gain
			N = N[R:]
			n = len(N)
		}
	}

	res = res[:k]

	return res
}
*/

func CombineIDs(N, M []int32) []int32 {

	if N == nil {
		return M
	}
	if M == nil {
		return N
	}

	n, m := len(N), len(M)

	// swap to make M the smaller list
	if n < m {
		N, M = M, N
		n, m = m, n
	}

	if m < 1 {
		return N
	}

	i, j, k := 0, 0, 0

	res := make([]int32, n+m)

	for i < n && j < m {
		if N[i] < M[j] {
			res[k] = N[i]
			k++
			i++
		} else if N[i] > M[j] {
			res[k] = M[j]
			k++
			j++
		} else {
			res[k] = N[i]
			k++
			i++
			j++
		}
	}
	for i < n {
		res[k] = N[i]
		k++
		i++
	}
	for j < m {
		res[k] = M[j]
		k++
		j++
	}

	res = res[:k]

	return res
}

func ExcludeIDs(N, M []int32) []int32 {

	if N == nil {
		return nil
	}
	if M == nil {
		return N
	}

	n, m := len(N), len(M)

	if m < 1 {
		return N
	}

	res := make([]int32, n)

	i, j, k := 0, 0, 0

	// use local variables for speed
	en, em := N[i], M[j]

	for {
		// do inequality tests first
		if en < em {
			// item is not excluded
			res[k] = en
			k++
			i++
			if i == n {
				break
			}
			en = N[i]
		} else if en > em {
			// advance second list pointer
			j++
			if j == m {
				break
			}
			em = M[j]
		} else {
			// exclude
			i++
			j++
			if i == n || j == m {
				break
			}
			en = N[i]
			em = M[j]
		}
	}

	// truncate output array to actual size of result
	res = res[:k]

	return res
}

// QUERY EVALUATION FUNCTION

func PostingIDsFuture(base, term, field string, dist int) <-chan Arrays {

	out := make(chan Arrays, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create postings channel\n")
		os.Exit(1)
	}

	// postingFuture asynchronously gets posting IDs and sends results through channel
	postingFuture := func(base, term, field string, dist int, out chan<- Arrays) {

		data, ofst := GetPostingIDs(base, term, field, false)

		out <- Arrays{Data: data, Ofst: ofst, Dist: dist}

		close(out)
	}

	// launch single future goroutine
	go postingFuture(base, term, field, dist, out)

	return out
}

func EvaluateQuery(base string, clauses []string) int {

	if clauses == nil || clauses[0] == "" {
		return 0
	}

	count := 0

	// flag set if no tildes, indicates no proximity tests in query
	noProx := true
	for _, tkn := range clauses {
		if strings.HasPrefix(tkn, "~") {
			noProx = false
		}
	}

	phrasePositions := func(pn, pm []int16, dlt int16) []int16 {

		var arry []int16

		ln, lm := len(pn), len(pm)

		q, r := 0, 0

		vn, vm := pn[q], pm[r]
		vnd := vn + dlt

		for {
			if vnd > vm {
				r++
				if r == lm {
					break
				}
				vm = pm[r]
			} else if vnd < vm {
				q++
				if q == ln {
					break
				}
				vn = pn[q]
				vnd = vn + dlt
			} else {
				// store position of first word in current growing phrase
				arry = append(arry, vn)
				q++
				r++
				if q == ln || r == lm {
					break
				}
				vn = pn[q]
				vm = pm[r]
				vnd = vn + dlt
			}
		}

		return arry
	}

	proximityPositions := func(pn, pm []int16, dlt int16) []int16 {

		var arry []int16

		ln, lm := len(pn), len(pm)

		q, r := 0, 0

		vn, vm := pn[q], pm[r]
		vnd := vn + dlt

		for {
			if vnd < vm {
				q++
				if q == ln {
					break
				}
				vn = pn[q]
				vnd = vn + dlt
			} else if vn < vm {
				// store position of first word in downstream phrase that passes proximity test
				arry = append(arry, vm)
				q++
				r++
				if q == ln || r == lm {
					break
				}
				vn = pn[q]
				vm = pm[r]
				vnd = vn + dlt
			} else {
				r++
				if r == lm {
					break
				}
				vm = pm[r]
			}
		}

		return arry
	}

	eval := func(str string) ([]int32, [][]int16, int) {

		// extract optional [FIELD] qualifier
		field := "NORM"

		slen := len(str)

		for _, key := range IdxFields {
			if strings.HasSuffix(str, " ["+key+"]") {
				klen := len(key) + 3
				str = str[:slen-klen]
				field = key
				switch key {
				case "CHEM", "CODE", "DISZ", "GENE", "PATH", "THME", "TREE":
					str = strings.Replace(str, " ", "_", -1)
				case "PIPE":
					// esearch -db pubmed -query "complement system proteins [MESH]" -pub clinical |
					// efetch -format uid | phrase-search -query "[PIPE] AND L [THME]"
					var data []int32
					// read UIDs from stdin
					uidq := CreateUIDReader(os.Stdin)
					for ext := range uidq {

						val, err := strconv.Atoi(ext.Text)
						if err != nil {
							fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized UID %s\n", ext.Text)
							os.Exit(1)
						}

						data = append(data, int32(val))
					}
					// sort UIDs before returning
					sort.Slice(data, func(i, j int) bool { return data[i] < data[j] })
					return data, nil, 0
				default:
				}
				break
			}
		}

		words := strings.Fields(str)

		if words == nil || len(words) < 1 {
			return nil, nil, 0
		}

		// if no tilde proximity tests, and not building up phrase from multiple words,
		// no need to use more expensive position tests when calculating intersection
		if noProx && len(words) == 1 {
			term := words[0]
			if strings.HasPrefix(term, "+") {
				return nil, nil, 0
			}
			term = strings.Replace(term, "_", " ", -1)
			data, _ := GetPostingIDs(base, term, field, true)
			count++
			return data, nil, 1
		}

		dist := 0

		var intersect []Arrays

		var futures []<-chan Arrays

		// schedule asynchronous fetching
		for _, term := range words {

			term = strings.Replace(term, "_", " ", -1)

			if strings.HasPrefix(term, "+") {
				dist += strings.Count(term, "+")
				// run of stop words or explicit plus signs skip past one or more words in phrase
				continue
			}

			fetch := PostingIDsFuture(base, term, field, dist)

			futures = append(futures, fetch)

			dist++
		}

		runtime.Gosched()

		for _, chn := range futures {

			// fetch postings data
			fut := <-chn

			if len(fut.Data) < 1 {
				// bail if word not present
				return nil, nil, 0
			}

			// append posting and positions
			intersect = append(intersect, fut)

			runtime.Gosched()
		}

		if len(intersect) < 1 {
			return nil, nil, 0
		}

		// start phrase with first word
		data, ofst, dist := intersect[0].Data, intersect[0].Ofst, intersect[0].Dist+1

		if len(intersect) == 1 {
			return data, ofst, dist
		}

		for i := 1; i < len(intersect); i++ {

			// add subsequent words, keep starting positions of phrases that contain all words in proper position
			data, ofst = ExtendPositionalIDs(data, ofst, intersect[i].Data, intersect[i].Ofst, intersect[i].Dist, phrasePositions)
			if len(data) < 1 {
				// bail if phrase not present
				return nil, nil, 0
			}
			dist = intersect[i].Dist + 1
		}

		count += len(intersect)

		// return UIDs and all positions of current phrase
		return data, ofst, dist
	}

	prevTkn := ""

	nextToken := func() string {

		if len(clauses) < 1 {
			return ""
		}

		// remove next token from slice
		tkn := clauses[0]
		clauses = clauses[1:]

		if tkn == "(" && prevTkn != "" && prevTkn != "&" && prevTkn != "|" && prevTkn != "!" {
			fmt.Fprintf(os.Stderr, "\nERROR: Tokens '%s' and '%s' should be separated by AND, OR, or NOT\n", prevTkn, tkn)
			os.Exit(1)
		}

		if prevTkn == ")" && tkn != "" && tkn != "&" && tkn != "|" && tkn != "!" && tkn != ")" {
			fmt.Fprintf(os.Stderr, "\nERROR: Tokens '%s' and '%s' should be separated by AND, OR, or NOT\n", prevTkn, tkn)
			os.Exit(1)
		}

		prevTkn = tkn

		return tkn
	}

	// recursive definitions
	var fact func() ([]int32, [][]int16, int, string)
	var prox func() ([]int32, string)
	var excl func() ([]int32, string)
	var term func() ([]int32, string)
	var expr func() ([]int32, string)

	fact = func() ([]int32, [][]int16, int, string) {

		var (
			data  []int32
			ofst  [][]int16
			delta int
			tkn   string
		)

		tkn = nextToken()

		if tkn == "(" {
			// recursively process expression in parentheses
			data, tkn = expr()
			if tkn == ")" {
				tkn = nextToken()
			} else {
				fmt.Fprintf(os.Stderr, "\nERROR: Expected ')' but received '%s'\n", tkn)
				os.Exit(1)
			}
		} else if tkn == ")" {
			fmt.Fprintf(os.Stderr, "\nERROR: Unexpected ')' token\n")
			os.Exit(1)
		} else if tkn == "&" || tkn == "|" || tkn == "!" {
			fmt.Fprintf(os.Stderr, "\nERROR: Unexpected operator '%s' in expression\n", tkn)
			os.Exit(1)
		} else if tkn == "" {
			fmt.Fprintf(os.Stderr, "\nERROR: Unexpected end of expression\n")
			os.Exit(1)
		} else {
			// evaluate current phrase
			data, ofst, delta = eval(tkn)
			tkn = nextToken()
		}

		return data, ofst, delta, tkn
	}

	prox = func() ([]int32, string) {

		var (
			next []int32
			noff [][]int16
			ndlt int
		)

		data, ofst, delta, tkn := fact()
		if len(data) < 1 {
			return nil, ""
		}

		for strings.HasPrefix(tkn, "~") {
			dist := strings.Count(tkn, "~")
			next, noff, ndlt, tkn = fact()
			if len(next) < 1 {
				return nil, ""
			}
			// next phrase must be within specified distance after the previous phrase
			data, ofst = ExtendPositionalIDs(data, ofst, next, noff, delta+dist, proximityPositions)
			if len(data) < 1 {
				return nil, ""
			}
			delta = ndlt
		}

		return data, tkn
	}

	excl = func() ([]int32, string) {

		var next []int32

		data, tkn := prox()
		for tkn == "!" {
			next, tkn = prox()
			data = ExcludeIDs(data, next)
		}

		return data, tkn
	}

	term = func() ([]int32, string) {

		var next []int32

		data, tkn := excl()
		for tkn == "&" {
			next, tkn = excl()
			data = IntersectIDs(data, next)
		}

		return data, tkn
	}

	expr = func() ([]int32, string) {

		var next []int32

		data, tkn := term()
		for tkn == "|" {
			next, tkn = term()
			data = CombineIDs(data, next)
		}

		return data, tkn
	}

	// enter recursive descent parser
	result, tkn := expr()

	if tkn != "" {
		fmt.Fprintf(os.Stderr, "\nERROR: Unexpected token '%s' at end of expression\n", tkn)
		os.Exit(1)
	}

	// sort final result
	sort.Slice(result, func(i, j int) bool { return result[i] < result[j] })

	// use buffers to speed up uid printing
	var buffer strings.Builder

	wrtr := bufio.NewWriter(os.Stdout)

	for _, pmid := range result {
		val := strconv.Itoa(int(pmid))
		buffer.WriteString(val[:])
		buffer.WriteString("\n")
	}

	txt := buffer.String()
	if txt != "" {
		// print buffer
		wrtr.WriteString(txt[:])
	}

	wrtr.Flush()

	runtime.Gosched()

	return count
}

// QUERY PARSING FUNCTIONS

func PrepareQuery(str string) string {

	if str == "" {
		return ""
	}

	if strings.HasPrefix(str, "[PIPE]") {
		str = "stdin " + str
	}

	// cleanup string
	if IsNotASCII(str) {
		str = DoAccentTransform(str)
		if HasUnicodeMarkup(str) {
			str = RepairUnicodeMarkup(str, SPACE)
		}
	}

	if HasBadSpace(str) {
		str = CleanupBadSpaces(str)
	}
	if HasAngleBracket(str) {
		str = RepairEncodedMarkup(str)
		str = RepairScriptMarkup(str, SPACE)
		str = RepairMathMLMarkup(str, SPACE)
		str = RemoveEmbeddedMarkup(str)
	}

	if HasAmpOrNotASCII(str) {
		str = html.UnescapeString(str)
	}

	if IsNotASCII(str) {
		if HasGreek(str) {
			str = SpellGreek(str)
			str = CompressRunsOfSpaces(str)
		}
		str = UnicodeToASCII(str)
	}

	str = strings.Replace(str, "~ ~", "~~", -1)
	str = strings.Replace(str, "~ ~", "~~", -1)

	str = strings.TrimSpace(str)

	// temporarily flank with spaces to detect misplaced operators at ends
	str = " " + str + " "

	str = strings.Replace(str, " AND ", " & ", -1)
	str = strings.Replace(str, " OR ", " | ", -1)
	str = strings.Replace(str, " NOT ", " ! ", -1)

	str = strings.Replace(str, "(", " ( ", -1)
	str = strings.Replace(str, ")", " ) ", -1)
	str = strings.Replace(str, "&", " & ", -1)
	str = strings.Replace(str, "|", " | ", -1)
	str = strings.Replace(str, "!", " ! ", -1)

	// ensure that bracketed fields are flanked by spaces
	str = strings.Replace(str, "[", " [", -1)
	str = strings.Replace(str, "]", "] ", -1)

	// remove temporary flanking spaces
	str = strings.TrimSpace(str)

	str = strings.ToLower(str)

	str = strings.Replace(str, "_", " ", -1)

	if HasPlusOrMinus(str) {
		str = FixThemeCases(str)
	}

	if HasHyphenOrApostrophe(str) {
		str = FixSpecialCases(str)
	}

	str = strings.Replace(str, "-", " ", -1)

	// convert bracketed fields to capitalized words
	str = DecodeFields(str)

	// break terms at punctuation, and at non-ASCII characters, allowing Boolean control symbols, along with
	// underscore for protected terms, asterisk to indicate truncation wildcard, tilde for maximum proximity,
	// and plus sign for exactly one wildcard word
	terms := strings.FieldsFunc(str, func(c rune) bool {
		return (!unicode.IsLetter(c) && !unicode.IsDigit(c) && c != '_' && c != '*' && c != '~' &&
			c != '+' && c != '$' && c != '&' && c != '|' && c != '!' && c != '(' && c != ')') || c > 127
	})

	// rejoin into processed sentence
	tmp := strings.Join(terms, " ")

	tmp = CompressRunsOfSpaces(tmp)
	tmp = strings.TrimSpace(tmp)

	return tmp
}

func PrepareExact(str string) string {

	if str == "" {
		return ""
	}

	if str == "[Not Available]." || str == "Health." {
		return ""
	}

	str = html.EscapeString(str)

	if IsNotASCII(str) {
		str = DoAccentTransform(str)
		if HasUnicodeMarkup(str) {
			str = RepairUnicodeMarkup(str, SPACE)
		}
	}

	str = strings.ToLower(str)

	if HasBadSpace(str) {
		str = CleanupBadSpaces(str)
	}
	if HasAngleBracket(str) {
		str = RepairEncodedMarkup(str)
		str = RepairScriptMarkup(str, SPACE)
		str = RepairMathMLMarkup(str, SPACE)
		// RemoveEmbeddedMarkup must be called before UnescapeString, which was suppressed in ExploreElements
		str = RemoveEmbeddedMarkup(str)
	}

	if HasAmpOrNotASCII(str) {
		str = html.UnescapeString(str)
	}

	if IsNotASCII(str) {
		if HasGreek(str) {
			str = SpellGreek(str)
			str = CompressRunsOfSpaces(str)
		}
	}

	str = strings.Replace(str, "(", " ", -1)
	str = strings.Replace(str, ")", " ", -1)

	str = strings.Replace(str, "_", " ", -1)

	if HasHyphenOrApostrophe(str) {
		str = FixSpecialCases(str)
	}

	str = strings.Replace(str, "-", " ", -1)

	// remove trailing punctuation from each word
	var arry []string

	terms := strings.Fields(str)
	for _, item := range terms {
		max := len(item)
		for max > 1 {
			ch := item[max-1]
			if ch != '.' && ch != ',' && ch != ':' && ch != ';' {
				break
			}
			// trim trailing period, comma, colon, and semicolon
			item = item[:max-1]
			// continue checking for runs of punctuation at end
			max--
		}
		if item == "" {
			continue
		}
		arry = append(arry, item)
	}

	// rejoin into string
	cleaned := strings.Join(arry, " ")

	// break clauses at punctuation other than space or underscore, and at non-ASCII characters
	clauses := strings.FieldsFunc(cleaned, func(c rune) bool {
		return (!unicode.IsLetter(c) && !unicode.IsDigit(c)) && c != ' ' && c != '_' || c > 127
	})

	// space replaces plus sign to separate runs of unpunctuated words
	phrases := strings.Join(clauses, " ")

	var chain []string

	// break phrases into individual words
	words := strings.Fields(phrases)

	for _, item := range words {

		// skip at site of punctuation break
		if item == "+" {
			chain = append(chain, "+")
			continue
		}

		// skip terms that are all digits
		if IsAllDigitsOrPeriod(item) {
			chain = append(chain, "+")
			continue
		}

		// optional stop word removal
		if DeStop && IsStopWord(item) {
			chain = append(chain, "+")
			continue
		}

		// index single normalized term
		chain = append(chain, item)
	}

	// rejoin into processed sentence
	tmp := strings.Join(chain, " ")

	tmp = strings.Replace(tmp, "+ +", "++", -1)
	tmp = strings.Replace(tmp, "+ +", "++", -1)

	tmp = CompressRunsOfSpaces(tmp)
	tmp = strings.TrimSpace(tmp)

	return tmp
}

func ProcessStopWords(str string) string {

	if str == "" {
		return ""
	}

	var chain []string

	terms := strings.Fields(str)

	nextField := func(terms []string) (string, int) {

		for j, item := range terms {
			for _, key := range IdxFields {
				if item == key {
					return item, j + 1
				}
			}
		}

		return "", 0
	}

	// replace unwanted and stop words with plus sign
	for len(terms) > 0 {

		item := terms[0]
		terms = terms[1:]

		fld, j := nextField(terms)

		switch fld {
		case "CHEM", "CODE", "DISZ", "GENE", "PATH", "THME", "TREE":
			for j > 0 {
				// do not treat non-TIAB terms as stop words
				chain = append(chain, item)
				j--
				item = terms[0]
				terms = terms[1:]
			}
			chain = append(chain, fld)
			continue
		case "PIPE":
			for j > 0 {
				// add leading stdin term
				chain = append(chain, item)
				j--
				item = terms[0]
				terms = terms[1:]
			}
			chain = append(chain, fld)
			continue
		default:
		}

		// skip if stop word, breaking phrase chain
		if DeStop && IsStopWord(item) {
			chain = append(chain, "+")
			continue
		}

		// index single normalized term
		chain = append(chain, item)
	}

	// rejoin into processed sentence
	tmp := strings.Join(chain, " ")

	tmp = strings.Replace(tmp, "+ +", "++", -1)
	tmp = strings.Replace(tmp, "+ +", "++", -1)

	tmp = strings.Replace(tmp, "~ +", "~+", -1)
	tmp = strings.Replace(tmp, "+ ~", "+~", -1)

	for strings.Contains(tmp, "~+") {
		tmp = strings.Replace(tmp, "~+", "~~", -1)
	}
	for strings.Contains(tmp, "+~") {
		tmp = strings.Replace(tmp, "+~", "~~", -1)
	}

	tmp = CompressRunsOfSpaces(tmp)
	tmp = strings.TrimSpace(tmp)

	return tmp
}

func PartitionQuery(str string) []string {

	if str == "" {
		return nil
	}

	str = CompressRunsOfSpaces(str)
	str = strings.TrimSpace(str)

	str = " " + str + " "

	// flank all operators with caret
	str = strings.Replace(str, " ( ", " ^ ( ^ ", -1)
	str = strings.Replace(str, " ) ", " ^ ) ^ ", -1)
	str = strings.Replace(str, " & ", " ^ & ^ ", -1)
	str = strings.Replace(str, " | ", " ^ | ^ ", -1)
	str = strings.Replace(str, " ! ", " ^ ! ^ ", -1)
	str = strings.Replace(str, " ~", " ^ ~", -1)
	str = strings.Replace(str, "~ ", "~ ^ ", -1)

	str = CompressRunsOfSpaces(str)
	str = strings.TrimSpace(str)

	str = strings.Replace(str, "^ ^", "^", -1)

	if strings.HasPrefix(str, "^ ") {
		str = str[2:]
	}
	if strings.HasSuffix(str, " ^") {
		max := len(str)
		str = str[:max-2]
	}

	str = strings.Replace(str, "~ ^ +", "~+", -1)
	str = strings.Replace(str, "+ ^ ~", "+~", -1)

	str = strings.Replace(str, "~ +", "~+", -1)
	str = strings.Replace(str, "+ ~", "+~", -1)

	for strings.Contains(str, "~+") {
		str = strings.Replace(str, "~+", "~~", -1)
	}
	for strings.Contains(str, "+~") {
		str = strings.Replace(str, "+~", "~~", -1)
	}

	// split into non-broken phrase segments or operator symbols
	tmp := strings.Split(str, " ^ ")

	return tmp
}

func SetFieldQualifiers(clauses []string, rlxd bool) []string {

	var res []string

	if clauses == nil {
		return nil
	}

	performStemming := func(str string, allowStemming bool) string {

		if str == "" {
			return ""
		}

		var chain []string

		terms := strings.Fields(str)

		// replace unwanted and stop words with plus sign
		for _, item := range terms {

			// allow tilde proximity indicator
			if item == "~" {
				chain = append(chain, item)
				continue
			}

			// skip terms that are all digits
			if IsAllDigitsOrPeriod(item) {
				chain = append(chain, "+")
				continue
			}

			// skip if stop word, breaking phrase chain
			if DeStop && IsStopWord(item) {
				chain = append(chain, "+")
				continue
			}

			// apply stemming algorithm
			if allowStemming {
				isWildCard := strings.HasSuffix(item, "*")
				if isWildCard {
					// temporarily remove trailing asterisk
					item = strings.TrimSuffix(item, "*")
				}

				item = porter2.Stem(item)
				item = strings.TrimSpace(item)

				if isWildCard {
					// do wildcard search in stemmed term list
					item += "*"
				}
			}

			// record single term
			chain = append(chain, item)
		}

		// rejoin into processed sentence
		tmp := strings.Join(chain, " ")

		tmp = strings.Replace(tmp, "+ +", "++", -1)
		tmp = strings.Replace(tmp, "+ +", "++", -1)

		tmp = CompressRunsOfSpaces(tmp)
		tmp = strings.TrimSpace(tmp)

		return tmp
	}

	for _, str := range clauses {

		// pass control symbols unchanged
		if str == "(" || str == ")" || str == "&" || str == "|" || str == "!" || strings.HasPrefix(str, "~") {
			res = append(res, str)
			continue
		}

		// pass angle bracket content delimiters (for -phrase, -require, -exclude)
		if str == "<" || str == ">" {
			res = append(res, str)
			continue
		}

		slen := len(str)

		allowStemming := false
		explicitNorm := false

		if strings.HasSuffix(str, " YEAR") {

			str = str[:slen-5]

			// regular 4-digit year
			if len(str) == 4 && IsAllDigitsOrPeriod(str) {
				res = append(res, str+" [YEAR]")
				continue
			}

			// check for year wildcard
			if len(str) == 4 && str[3] == '*' && IsAllDigitsOrPeriod(str[:3]) {

				fmt.Fprintf(os.Stderr, "\nERROR: Wildcards not supported for years - use ####:#### range instead\n")
				os.Exit(1)
			}

			// check for year range
			if len(str) == 9 && str[4] == ' ' && IsAllDigitsOrPeriod(str[:4]) && IsAllDigitsOrPeriod(str[5:]) {
				start, err := strconv.Atoi(str[:4])
				if err != nil {
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to recognize starting year '%s'\n", str[:4])
					os.Exit(1)
				}
				stop, err := strconv.Atoi(str[5:])
				if err != nil {
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to recognize stopping year '%s'\n", str[5:])
					os.Exit(1)
				}
				if start > stop {
					continue
				}
				// expand year range into individual year-by-year queries
				pfx := "("
				sfx := ")"
				for start <= stop {
					res = append(res, pfx)
					pfx = "|"
					yr := strconv.Itoa(start)
					res = append(res, yr+" [YEAR]")
					start++
				}
				res = append(res, sfx)
				continue
			}

			fmt.Fprintf(os.Stderr, "\nERROR: Unable to recognize year expression '%s'\n", str)
			os.Exit(1)

		} else if strings.HasSuffix(str, " TREE") {

			str = str[:slen-5]
			str = strings.Replace(str, " ", ".", -1)
			tmp := str
			tmp = strings.TrimSuffix(tmp, "*")
			if len(tmp) > 2 && unicode.IsLower(rune(tmp[0])) && IsAllDigitsOrPeriod(tmp[1:]) {
				str = strings.Replace(str, ".", " ", -1)
				res = append(res, str+" [TREE]")
				continue
			}

			fmt.Fprintf(os.Stderr, "\nERROR: Unable to recognize mesh code expression '%s'\n", str)
			os.Exit(1)

		} else if strings.HasSuffix(str, " CODE") {

			str = str[:slen-5]
			res = append(res, str+" [CODE]")
			continue

		} else if strings.HasSuffix(str, " THME") {

			str = str[:slen-5]
			res = append(res, str+" [THME]")
			continue

		} else if strings.HasSuffix(str, " PATH") {

			str = str[:slen-5]
			res = append(res, str+" [PATH]")
			continue

		} else if strings.HasSuffix(str, " STEM") {

			allowStemming = true
			str = str[:slen-5]

		} else if strings.HasSuffix(str, " NORM") {

			explicitNorm = true
			str = str[:slen-5]

		} else {

			found := false
			for _, key := range IdxFields {
				if strings.HasSuffix(str, " "+key) {
					klen := len(key) + 1
					str = str[:slen-klen]
					res = append(res, str+" ["+key+"]")
					found = true
					break
				}
			}
			if found {
				continue
			}

			if rlxd {

				allowStemming = true
			}
		}

		// process clause, using plus sign to break runs of words
		tmp := performStemming(str, allowStemming)

		// remove leading and trailing plus signs and spaces
		for strings.HasPrefix(tmp, "+") || strings.HasPrefix(tmp, " ") {
			tmp = tmp[1:]
		}
		for strings.HasSuffix(tmp, "+") || strings.HasSuffix(tmp, " ") {
			tlen := len(tmp)
			tmp = tmp[:tlen-1]
		}

		if allowStemming {
			tmp += " [STEM]"
		} else if explicitNorm {
			tmp += " [NORM]"
		}

		res = append(res, tmp)
	}

	return res
}

// SEARCH TERM LISTS FOR PHRASES OR NORMALIZED TERMS, OR MATCH BY PATTERN

func ProcessSearch(base, phrase string, xact, rlxd bool) int {

	if phrase == "" {
		return 0
	}

	if xact {
		phrase = PrepareExact(phrase)
	} else {
		phrase = PrepareQuery(phrase)
	}

	phrase = ProcessStopWords(phrase)

	clauses := PartitionQuery(phrase)

	clauses = SetFieldQualifiers(clauses, rlxd)

	return EvaluateQuery(base, clauses)
}

func ProcessMock(base, phrase string, xact, rlxd bool) int {

	if phrase == "" {
		return 0
	}

	fmt.Fprintf(os.Stdout, "ProcessSearch:\n\n%s\n\n", phrase)

	if xact {
		phrase = PrepareExact(phrase)

		fmt.Fprintf(os.Stdout, "PrepareExact:\n\n%s\n\n", phrase)
	} else {
		phrase = PrepareQuery(phrase)

		fmt.Fprintf(os.Stdout, "PrepareQuery:\n\n%s\n\n", phrase)
	}

	phrase = ProcessStopWords(phrase)

	fmt.Fprintf(os.Stdout, "ProcessStopWords:\n\n%s\n\n", phrase)

	clauses := PartitionQuery(phrase)

	fmt.Fprintf(os.Stdout, "PartitionQuery:\n\n")
	for _, tkn := range clauses {
		fmt.Fprintf(os.Stdout, "%s\n", tkn)
	}
	fmt.Fprintf(os.Stdout, "\n")

	clauses = SetFieldQualifiers(clauses, rlxd)

	fmt.Fprintf(os.Stdout, "SetFieldQualifiers:\n\n")
	for _, tkn := range clauses {
		fmt.Fprintf(os.Stdout, "%s\n", tkn)
	}
	fmt.Fprintf(os.Stdout, "\n")

	return 0
}

func ProcessCount(base, phrase string, plrl, psns, rlxd bool) int {

	if phrase == "" {
		return 0
	}

	phrase = PrepareQuery(phrase)

	phrase = ProcessStopWords(phrase)

	clauses := PartitionQuery(phrase)

	clauses = SetFieldQualifiers(clauses, rlxd)

	if clauses == nil {
		return 0
	}

	count := 0

	splitIntoWords := func(str string) []string {

		if str == "" {
			return nil
		}

		var arry []string

		parts := strings.Split(str, "+")

		for _, segment := range parts {

			segment = strings.TrimSpace(segment)

			if segment == "" {
				continue
			}

			words := strings.Fields(segment)

			for _, item := range words {
				if strings.HasPrefix(item, "~") {
					continue
				}
				arry = append(arry, item)
			}
		}

		return arry
	}

	parseField := func(str string) (string, string) {

		field := "NORM"

		slen := len(str)

		for _, key := range IdxFields {
			if strings.HasSuffix(str, " ["+key+"]") {
				klen := len(key) + 3
				str = str[:slen-klen]
				field = key
				switch key {
				case "CHEM", "CODE", "DISZ", "GENE", "PATH", "THME", "TREE":
					str = strings.Replace(str, " ", "_", -1)
				case "PIPE":
				default:
				}
				break
			}
		}

		return field, str
	}

	checkTermCounts := func(txt string) {

		field, str := parseField(txt)

		var words []string

		words = splitIntoWords(str)

		if words == nil || len(words) < 1 {
			return
		}

		for _, term := range words {

			term = strings.Replace(term, "_", " ", -1)

			if psns {
				count += PrintTermPositions(base, term, field)
			} else if plrl {
				count += PrintTermCounts(base, term, field)
			} else {
				count += PrintTermCount(base, term, field)
			}
		}
	}

	for _, item := range clauses {

		// skip control symbols
		if item == "(" || item == ")" || item == "&" || item == "|" || item == "!" {
			continue
		}

		checkTermCounts(item)
	}

	runtime.Gosched()

	return count
}

// READ FILE IN ARCHIVE THAT STORES UIDS WITH NON-DEFAULT VERSIONS

// CONCURRENT GOROUTINE SERVERS

// processes with single goroutine call defer close(out) so consumer(s) can range over channel
// processes with multiple instances call defer wg.Done(), separate goroutine uses wg.Wait() to delay close(out)

func CreateUIDReader(in io.Reader) <-chan Extract {

	if in == nil {
		return nil
	}

	out := make(chan Extract, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create uid reader channel\n")
		os.Exit(1)
	}

	// uidReader reads uids from input stream and sends through channel
	uidReader := func(in io.Reader, out chan<- Extract) {

		// close channel when all records have been processed
		defer close(out)

		scanr := bufio.NewScanner(in)

		idx := 0
		for scanr.Scan() {

			// read lines of identifiers
			file := scanr.Text()
			idx++

			pos := strings.Index(file, ".")
			if pos >= 0 {
				// remove version suffix
				file = file[:pos]
			}

			out <- Extract{idx, "", file, nil}
		}
	}

	// launch single uid reader goroutine
	go uidReader(in, out)

	return out
}

func CreateStashers(stash, parent, indx, sfx string, hash, zipp bool, report int, inp <-chan Extract) <-chan string {

	if inp == nil {
		return nil
	}

	out := make(chan string, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create stasher channel\n")
		os.Exit(1)
	}

	find := ParseIndex(indx)

	if zipp {
		sfx += ".gz"
	}

	type StasherType int

	const (
		OKAY StasherType = iota
		WAIT
		BAIL
	)

	// mutex to protect access to inUse map
	var wlock sync.Mutex

	// map to track files currently being written
	inUse := make(map[string]int)

	// lockFile function prevents colliding writes
	lockFile := func(id string, index int) StasherType {

		// map is non-reentrant, protect with mutex
		wlock.Lock()

		// multiple return paths, schedule the unlock command up front
		defer wlock.Unlock()

		idx, ok := inUse[id]

		if ok {
			if idx < index {
				// later version is being written by another goroutine, skip this
				return BAIL
			}
			// earlier version is being written by another goroutine, wait
			return WAIT
		}

		// okay to write file, mark in use to prevent collision
		inUse[id] = index
		return OKAY
	}

	// freeFile function removes entry from inUse map
	freeFile := func(id string) {

		wlock.Lock()

		// free entry in map, later versions of same record can now be written
		delete(inUse, id)

		wlock.Unlock()
	}

	// mutex to protect access to rollingCount variable
	var tlock sync.Mutex

	rollingCount := 0

	countSuccess := func() {

		tlock.Lock()

		rollingCount++
		if rollingCount >= report {
			rollingCount = 0
			// print dot (progress monitor)
			fmt.Fprintf(os.Stderr, ".")
		}

		tlock.Unlock()
	}

	// stashRecord saves individual XML record to archive file accessed by trie
	stashRecord := func(str, id string, index int) string {

		pos := strings.Index(id, ".")
		if pos >= 0 {
			// remove version from UID
			id = id[:pos]
		}

		var arry [132]rune
		trie := MakeArchiveTrie(id, arry)
		if trie == "" {
			return ""
		}

		attempts := 5
		keepChecking := true

		for keepChecking {
			// check if file is not being written by another goroutine
			switch lockFile(id, index) {
			case OKAY:
				// okay to save this record now
				keepChecking = false
			case WAIT:
				// earlier version is being saved, wait one second and try again
				time.Sleep(time.Second)
				attempts--
				if attempts < 1 {
					// could not get lock after several attempts
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to save '%s'\n", id)
					return ""
				}
			case BAIL:
				// later version is being saved, skip this one
				return ""
			default:
			}
		}

		// delete lock after writing file
		defer freeFile(id)

		dpath := path.Join(stash, trie)
		if dpath == "" {
			return ""
		}
		err := os.MkdirAll(dpath, os.ModePerm)
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			return ""
		}
		fpath := path.Join(dpath, id+sfx)
		if fpath == "" {
			return ""
		}

		// overwrites and truncates existing file
		fl, err := os.Create(fpath)
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			return ""
		}

		res := ""

		if hash {
			// calculate hash code for verification table
			hsh := crc32.NewIEEE()
			hsh.Write([]byte(str))
			val := hsh.Sum32()
			res = strconv.FormatUint(uint64(val), 10)
		}

		if zipp {

			zpr, err := gzip.NewWriterLevel(fl, gzip.DefaultCompression)

			if err == nil {

				wrtr := bufio.NewWriter(zpr)

				// compress and copy record to file
				wrtr.WriteString(str)
				if !strings.HasSuffix(str, "\n") {
					wrtr.WriteString("\n")
				}

				err = wrtr.Flush()
				if err != nil {
					fmt.Fprintf(os.Stderr, "%s\n", err.Error())
					return ""
				}

				err = zpr.Close()
				if err != nil {
					fmt.Fprintf(os.Stderr, "%s\n", err.Error())
					return ""
				}

			} else {

				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			}

		} else {

			// copy uncompressed record to file
			fl.WriteString(str)
			if !strings.HasSuffix(str, "\n") {
				fl.WriteString("\n")
			}
		}

		// fl.Sync()

		err = fl.Close()
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			return ""
		}

		// progress monitor prints dot every 1000 (.xml) or 50000 (.e2x) records
		countSuccess()

		return res
	}

	// xmlStasher reads from channel and calls stashRecord
	xmlStasher := func(wg *sync.WaitGroup, inp <-chan Extract, out chan<- string) {

		defer wg.Done()

		for ext := range inp {

			ext.Ident = FindIdentifier(ext.Text, parent, find)

			hsh := stashRecord(ext.Text, ext.Ident, ext.Index)

			res := ext.Ident
			if hash {
				res += "\t" + hsh
			}
			res += "\n"

			runtime.Gosched()

			out <- res
		}
	}

	var wg sync.WaitGroup

	// launch multiple stasher goroutines
	for i := 0; i < NumServe; i++ {
		wg.Add(1)
		go xmlStasher(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all stashers are done
	go func() {
		wg.Wait()
		close(out)
		// print newline after rows of dots (progress monitor)
		fmt.Fprintf(os.Stderr, "\n")
	}()

	return out
}

func CreateFetchers(stash, sfx string, zipp bool, inp <-chan Extract) <-chan Extract {

	if inp == nil {
		return nil
	}

	out := make(chan Extract, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create fetcher channel\n")
		os.Exit(1)
	}

	if zipp {
		sfx += ".gz"
	}

	fetchRecord := func(file string, buf bytes.Buffer) string {

		var arry [132]rune
		trie := MakeArchiveTrie(file, arry)

		if file == "" || trie == "" {
			return ""
		}

		fpath := path.Join(stash, trie, file+sfx)
		if fpath == "" {
			return ""
		}

		iszip := zipp

		inFile, err := os.Open(fpath)

		// if failed to find ".xml" or ".e2x" file, try ".xml.gz" or ".e2x.gz" without requiring -gzip
		if err != nil && os.IsNotExist(err) && !zipp {
			iszip = true
			fpath := path.Join(stash, trie, file+sfx+".gz")
			if fpath == "" {
				return ""
			}
			inFile, err = os.Open(fpath)
		}
		if err != nil {
			msg := err.Error()
			if !strings.HasSuffix(msg, "no such file or directory") && !strings.HasSuffix(msg, "cannot find the path specified.") {
				fmt.Fprintf(os.Stderr, "%s\n", msg)
			}
			return ""
		}

		defer inFile.Close()

		brd := bufio.NewReader(inFile)

		if iszip {

			zpr, err := gzip.NewReader(brd)

			defer zpr.Close()

			if err == nil {
				// copy and decompress cached file contents
				buf.ReadFrom(zpr)
			}

		} else {

			// copy cached file contents
			buf.ReadFrom(brd)
		}

		str := buf.String()

		return str
	}

	// xmlFetcher reads XML from file
	xmlFetcher := func(wg *sync.WaitGroup, inp <-chan Extract, out chan<- Extract) {

		// report when more records to process
		defer wg.Done()

		var buf bytes.Buffer

		for ext := range inp {

			buf.Reset()

			str := fetchRecord(ext.Text, buf)

			runtime.Gosched()

			out <- Extract{ext.Index, "", str, nil}
		}
	}

	var wg sync.WaitGroup

	// launch multiple fetcher goroutines
	for i := 0; i < NumServe; i++ {
		wg.Add(1)
		go xmlFetcher(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all fetchers are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

func CreateStreamers(stash string, inp <-chan Extract) <-chan Extract {

	if inp == nil {
		return nil
	}

	out := make(chan Extract, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create streamer channel\n")
		os.Exit(1)
	}

	sfx := ".xml.gz"

	getRecord := func(file string, buf bytes.Buffer) []byte {

		var arry [132]rune
		trie := MakeArchiveTrie(file, arry)

		if file == "" || trie == "" {
			return nil
		}

		fpath := path.Join(stash, trie, file+sfx)
		if fpath == "" {
			return nil
		}

		inFile, err := os.Open(fpath)

		if err != nil {
			msg := err.Error()
			if !strings.HasSuffix(msg, "no such file or directory") && !strings.HasSuffix(msg, "cannot find the path specified.") {
				fmt.Fprintf(os.Stderr, "%s\n", msg)
			}
			return nil
		}

		brd := bufio.NewReader(inFile)

		// copy cached file contents
		buf.ReadFrom(brd)

		data := buf.Bytes()

		inFile.Close()

		return data
	}

	// xmlStreamer reads compressed XML from file
	xmlStreamer := func(wg *sync.WaitGroup, inp <-chan Extract, out chan<- Extract) {

		// report when more records to process
		defer wg.Done()

		var buf bytes.Buffer

		for ext := range inp {

			buf.Reset()

			data := getRecord(ext.Text, buf)

			runtime.Gosched()

			out <- Extract{ext.Index, "", "", data}
		}
	}

	var wg sync.WaitGroup

	// launch multiple streamer goroutines
	for i := 0; i < NumServe; i++ {
		wg.Add(1)
		go xmlStreamer(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all streamers are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

func CreateDispensers(inp <-chan Extract) <-chan []string {

	if inp == nil {
		return nil
	}

	out := make(chan []string, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create dispenser channel\n")
		os.Exit(1)
	}

	var ilock sync.Mutex

	// map for inverted index
	inverted := make(map[string][]string)

	// add single posting
	addPost := func(fld, term, pos, uid string) {

		// protect map with mutex
		ilock.Lock()

		data, ok := inverted[term]
		if !ok {
			data = make([]string, 0, 4)
			// first entry on new slice is term
			data = append(data, term)
		}
		data = append(data, fld)
		data = append(data, uid)
		data = append(data, pos)
		// always need to update inverted, since data may be reallocated
		inverted[term] = data

		// unlock at end to avoid defer overhead
		ilock.Unlock()
	}

	// xmlDispenser prepares UID, term, and position strings for inversion
	xmlDispenser := func(wg *sync.WaitGroup, inp <-chan Extract, out chan<- []string) {

		defer wg.Done()

		currUID := ""

		doDispense := func(tag, attr, content string) {

			if tag == "IdxUid" {
				currUID = content
			} else {

				// expand Greek letters, anglicize characters in other alphabets
				if IsNotASCII(content) {

					if HasGreek(content) {
						content = SpellGreek(content)
						content = CompressRunsOfSpaces(content)
					}
					content = unidecode.Unidecode(content)

					content = DoAccentTransform(content)
					content = UnicodeToASCII(content)

					content = strings.TrimSpace(content)
				}

				// remove punctuation from term
				content = strings.Map(func(c rune) rune {
					if !unicode.IsLetter(c) && !unicode.IsDigit(c) && c != ' ' && c != '-' && c != '_' {
						return -1
					}
					return c
				}, content)

				content = strings.Replace(content, "_", " ", -1)
				content = strings.Replace(content, "-", " ", -1)

				content = CompressRunsOfSpaces(content)
				content = strings.TrimSpace(content)

				if content != "" && currUID != "" {
					addPost(tag, content, attr, currUID)
				}
			}
		}

		// read partitioned XML from producer channel
		for ext := range inp {

			StreamValues(ext.Text[:], "IdxDocument", doDispense)
		}
	}

	var wg sync.WaitGroup

	// launch multiple dispenser goroutines
	for i := 0; i < NumServe; i++ {
		wg.Add(1)
		go xmlDispenser(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all dispensers are done
	go func() {
		wg.Wait()

		// send results to inverters
		for _, data := range inverted {
			out <- data

			runtime.Gosched()
		}

		close(out)
	}()

	return out
}

func CreateInverters(inp <-chan []string) <-chan Extract {

	if inp == nil {
		return nil
	}

	out := make(chan Extract, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create inverter channel\n")
		os.Exit(1)
	}

	// xmlInverter sorts and prints one posting list
	xmlInverter := func(wg *sync.WaitGroup, inp <-chan []string, out chan<- Extract) {

		defer wg.Done()

		var buffer strings.Builder

		printPosting := func(key string, data []string) string {

			fields := make(map[string]map[string]string)

			for len(data) > 1 {
				fld := data[0]
				uid := data[1]
				att := data[2]
				positions, ok := fields[fld]
				if !ok {
					positions = make(map[string]string)
					fields[fld] = positions
				}
				// store position attribute string by uid
				positions[uid] = att
				// skip to next position
				data = data[3:]
			}

			buffer.Reset()

			buffer.WriteString("  <InvDocument>\n")
			buffer.WriteString("    <InvKey>")
			buffer.WriteString(key)
			buffer.WriteString("</InvKey>\n")
			buffer.WriteString("    <InvIDs>\n")

			// sort fields in alphabetical order
			var keys []string
			for ky := range fields {
				keys = append(keys, ky)
			}
			sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })

			for _, fld := range keys {

				positions := fields[fld]

				var arry []string

				for item := range positions {
					arry = append(arry, item)
				}

				if len(arry) > 1 {
					sort.Slice(arry, func(i, j int) bool {
						// numeric sort on strings checks lengths first
						lni := len(arry[i])
						lnj := len(arry[j])
						// shorter string is numerically less, assuming no leading zeros
						if lni < lnj {
							return true
						}
						if lni > lnj {
							return false
						}
						// same length, can now do string comparison on contents
						return arry[i] < arry[j]
					})
				}

				// print list of UIDs, skipping duplicates
				prev := ""
				for _, uid := range arry {
					if uid == prev {
						continue
					}

					buffer.WriteString("      <")
					buffer.WriteString(fld)
					atr := positions[uid]
					if atr != "" {
						buffer.WriteString(" ")
						buffer.WriteString(atr)
					}
					buffer.WriteString(">")
					buffer.WriteString(uid)
					buffer.WriteString("</")
					buffer.WriteString(fld)
					buffer.WriteString(">\n")

					prev = uid
				}
			}

			buffer.WriteString("    </InvIDs>\n")
			buffer.WriteString("  </InvDocument>\n")

			str := buffer.String()

			return str
		}

		for inv := range inp {

			key := inv[0]
			data := inv[1:]

			str := printPosting(key, data)

			out <- Extract{0, key, str, nil}

			runtime.Gosched()
		}
	}

	var wg sync.WaitGroup

	// launch multiple inverter goroutines
	for i := 0; i < NumServe; i++ {
		wg.Add(1)
		go xmlInverter(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all inverters are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

func CreateResolver(inp <-chan Extract) <-chan string {

	if inp == nil {
		return nil
	}

	out := make(chan string, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create resolver channel\n")
		os.Exit(1)
	}

	// xmlResolver prints inverted postings alphabetized by identifier prefix
	xmlResolver := func(inp <-chan Extract, out chan<- string) {

		// close channel when all records have been processed
		defer close(out)

		// map for inverted index
		inverted := make(map[string]string)

		// drain channel, populate map for alphabetizing
		for curr := range inp {

			inverted[curr.Ident] = curr.Text
		}

		var ordered []string

		for item := range inverted {
			ordered = append(ordered, item)
		}

		if len(ordered) > 1 {
			sort.Slice(ordered, func(i, j int) bool { return ordered[i] < ordered[j] })
		}

		// iterate through alphabetized results
		for _, curr := range ordered {

			txt := inverted[curr]

			// send result to output
			out <- txt

			runtime.Gosched()
		}
	}

	// launch single resolver goroutine
	go xmlResolver(inp, out)

	return out
}

type Plex struct {
	Which int
	Ident string
	Text  string
	Index int
	Sibs  []string
}

type PlexHeap []Plex

// methods that satisfy heap.Interface
func (h PlexHeap) Len() int {
	return len(h)
}
func (h PlexHeap) Less(i, j int) bool {
	return h[i].Ident < h[j].Ident
}
func (h PlexHeap) Swap(i, j int) {
	h[i], h[j] = h[j], h[i]
}
func (h *PlexHeap) Push(x interface{}) {
	*h = append(*h, x.(Plex))
}
func (h *PlexHeap) Pop() interface{} {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[0 : n-1]
	return x
}

func CreatePresenters(args []string) []<-chan Plex {

	if args == nil {
		return nil
	}

	numFiles := len(args)
	if numFiles < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Not enough inverted files to merge\n")
		os.Exit(1)
	}

	chns := make([]<-chan Plex, numFiles)
	if chns == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create presenter channel array\n")
		os.Exit(1)
	}

	// xmlPresenter sends partitioned XML strings through channel
	xmlPresenter := func(fileNum int, fileName string, out chan<- Plex) {

		// close channel when all records have been processed
		defer close(out)

		f, err := os.Open(fileName)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open input file '%s'\n", fileName)
			os.Exit(1)
		}

		// close input file when all records have been processed
		defer f.Close()

		var in io.Reader

		in = f

		// if suffix is ".gz", use decompressor
		iszip := false
		if strings.HasSuffix(fileName, ".gz") {
			iszip = true
		}

		if iszip {
			brd := bufio.NewReader(f)
			if brd == nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create buffered reader on '%s'\n", fileName)
				os.Exit(1)
			}
			zpr, err := gzip.NewReader(brd)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create decompressor on '%s'\n", fileName)
				os.Exit(1)
			}

			// close decompressor when all records have been processed
			defer zpr.Close()

			// use decompressor for reading file
			in = zpr
		}

		rdr := CreateReader(in)

		if rdr == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create XML Block Reader\n")
			os.Exit(1)
		}

		find := ParseIndex("InvKey")

		// partition all input by pattern and send XML substring through channel
		PartitionPattern("InvDocument", "", rdr,
			func(str string) {
				id := FindIdentifier(str[:], "InvDocument", find)

				out <- Plex{fileNum, id, str, 0, nil}
			})
	}

	// launch multiple presenter goroutines
	for i, str := range args {

		chn := make(chan Plex, ChanDepth)
		if chn == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create presenter channel\n")
			os.Exit(1)
		}

		go xmlPresenter(i, str, chn)

		chns[i] = chn
	}

	// no need for separate anonymous goroutine to wait until all presenters are done

	return chns
}

func CreateManifold(inp []<-chan Plex) <-chan Plex {

	if inp == nil {
		return nil
	}

	out := make(chan Plex, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create manifold channel\n")
		os.Exit(1)
	}

	// xmlManifold restores alphabetical order of merged postings
	xmlManifold := func(inp []<-chan Plex, out chan<- Plex) {

		// close channel when all records have been processed
		defer close(out)

		// initialize empty heap
		hp := &PlexHeap{}
		heap.Init(hp)

		// read first object from all input channels in turn
		for _, chn := range inp {
			plx, ok := <-chn
			if ok {
				heap.Push(hp, plx)
			}
		}

		// array to collect strings with same identifier
		var arry []string

		prevIdent := ""
		rec := 0

		// reading from heap returns objects in alphabetical order
		for hp.Len() > 0 {

			// remove lowest item from heap, use interface type assertion
			curr := heap.Pop(hp).(Plex)

			// compare adjacent record identifiers
			if prevIdent == curr.Ident {

				// save next inverted object string in slice
				arry = append(arry, curr.Text)

			} else {

				if len(arry) > 0 {

					rec++
					// send set from previous identifier to output channel
					out <- Plex{0, prevIdent, "", rec, arry}

					// empty the slice
					arry = nil
				}

				// remember new identifier
				prevIdent = curr.Ident

				// save first inverted object with this identifier
				arry = append(arry, curr.Text)
			}

			// read next object from channel that just supplied lowest item
			chn := inp[curr.Which]
			plx, ok := <-chn
			if ok {
				heap.Push(hp, plx)
			}
		}

		if len(arry) > 0 {

			rec++
			// send last record
			out <- Plex{0, prevIdent, "", rec, arry}

			arry = nil
		}
	}

	// launch single manifold goroutine
	go xmlManifold(inp, out)

	return out
}

func CreateFusers(inp <-chan Extract) <-chan Plex {

	if inp == nil {
		return nil
	}

	out := make(chan Plex, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create fuser channel\n")
		os.Exit(1)
	}

	var flock sync.Mutex

	// map for combining inverted indices
	inverted := make(map[string][]string)

	addInverts := func(id, str string) {

		// protect map with mutex
		flock.Lock()

		data, ok := inverted[id]
		if !ok {
			data = make([]string, 0, 1)
		}

		data = append(data, str)
		// always need to update inverted, since data may be reallocated
		inverted[id] = data

		// unlock at end to avoid defer overhead
		flock.Unlock()
	}

	xmlFuser := func(wg *sync.WaitGroup, inp <-chan Extract, out chan<- Plex) {

		defer wg.Done()

		find := ParseIndex("InvKey")

		// read partitioned XML from producer channel
		for ext := range inp {

			str := ext.Text[:]
			id := FindIdentifier(str, "InvDocument", find)
			addInverts(id, str)
		}
	}

	var wg sync.WaitGroup

	// launch multiple fuser goroutines
	for i := 0; i < NumServe; i++ {
		wg.Add(1)
		go xmlFuser(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all fusers are done
	go func() {
		wg.Wait()

		// sort id keys in alphabetical order
		var keys []string
		for ky := range inverted {
			keys = append(keys, ky)
		}
		sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })

		rec := 0

		for _, id := range keys {

			arry := inverted[id]

			rec++
			// send array of records with same identifier to output channel
			out <- Plex{0, id, "", rec, arry}

			// empty the slice
			arry = nil
		}

		close(out)
	}()

	return out
}

func CreateMergers(inp <-chan Plex) <-chan Extract {

	if inp == nil {
		return nil
	}

	out := make(chan Extract, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create merger channel\n")
		os.Exit(1)
	}

	// xmlMerger fuses adjacent InvDocument records with the same identifier
	xmlMerger := func(wg *sync.WaitGroup, inp <-chan Plex, out chan<- Extract) {

		defer wg.Done()

		var buffer strings.Builder

		fusePostings := func(key string, data []string) string {

			fields := make(map[string]map[string]string)

			addIdents := func(fld, pos, uid string) {

				// no need for mutex here
				positions, ok := fields[fld]
				if !ok {
					positions = make(map[string]string)
					fields[fld] = positions
				}

				positions[uid] = pos
			}

			addUID := func(tag, attr, content string) {

				if tag != "InvKey" {

					addIdents(tag, attr, content)
				}
			}

			for _, str := range data {
				StreamValues(str[:], "InvDocument", addUID)
			}

			buffer.Reset()

			buffer.WriteString("  <InvDocument>\n")
			buffer.WriteString("    <InvKey>")
			buffer.WriteString(key)
			buffer.WriteString("</InvKey>\n")
			buffer.WriteString("    <InvIDs>\n")

			// sort fields in alphabetical order
			var keys []string
			for ky := range fields {
				keys = append(keys, ky)
			}
			sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })

			for _, fld := range keys {

				positions := fields[fld]

				var arry []string

				for item := range positions {
					arry = append(arry, item)
				}

				if len(arry) > 1 {
					sort.Slice(arry, func(i, j int) bool {
						// numeric sort on strings checks lengths first
						lni := len(arry[i])
						lnj := len(arry[j])
						// shorter string is numerically less, assuming no leading zeros
						if lni < lnj {
							return true
						}
						if lni > lnj {
							return false
						}
						// same length, can now do string comparison on contents
						return arry[i] < arry[j]
					})
				}

				// print list of UIDs, skipping duplicates
				last := ""
				for _, uid := range arry {
					// detect duplicate UIDs, now in same list after conversion of one term entry from foreign alphabet
					if uid == last {
						continue
					}
					buffer.WriteString("      <")
					buffer.WriteString(fld)
					atr := positions[uid]
					if atr != "" {
						buffer.WriteString(" ")
						buffer.WriteString(atr)
					}
					buffer.WriteString(">")
					buffer.WriteString(uid)
					buffer.WriteString("</")
					buffer.WriteString(fld)
					buffer.WriteString(">\n")

					last = uid
				}
			}

			buffer.WriteString("    </InvIDs>\n")
			buffer.WriteString("  </InvDocument>\n")

			txt := buffer.String()

			return txt
		}

		for plx := range inp {

			rec := plx.Index
			key := plx.Ident
			data := plx.Sibs

			if len(data) < 1 {
				continue
			}

			str := fusePostings(key, data)

			out <- Extract{rec, key, str, nil}

			runtime.Gosched()
		}
	}

	var wg sync.WaitGroup

	// launch multiple merger goroutines
	for i := 0; i < NumServe; i++ {
		wg.Add(1)
		go xmlMerger(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all mergers are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

func CreateSplitter(merg string, zipp bool, inp <-chan Extract) <-chan string {

	if inp == nil {
		return nil
	}

	out := make(chan string, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create splitter channel\n")
		os.Exit(1)
	}

	openSaver := func(merg, key string, zipp bool) (*os.File, *bufio.Writer, *gzip.Writer) {

		var (
			fl   *os.File
			wrtr *bufio.Writer
			zpr  *gzip.Writer
			err  error
		)

		sfx := ".mrg"
		if zipp {
			sfx += ".gz"
		}

		fpath := path.Join(merg, key+sfx)
		if fpath == "" {
			return nil, nil, nil
		}

		// overwrites and truncates existing file
		fl, err = os.Create(fpath)
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			return nil, nil, nil
		}

		var out io.Writer

		out = fl

		if zipp {

			zpr, err = gzip.NewWriterLevel(fl, gzip.BestSpeed)
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
				return nil, nil, nil
			}

			out = zpr
		}

		// create buffered writer layer
		wrtr = bufio.NewWriter(out)
		if wrtr == nil {
			fmt.Fprintf(os.Stderr, "Unable to create bufio.NewWriter\n")
			return nil, nil, nil
		}

		return fl, wrtr, zpr
	}

	closeSaver := func(fl *os.File, wrtr *bufio.Writer, zpr *gzip.Writer) {

		wrtr.Flush()
		if zpr != nil {
			zpr.Close()
		}
		// fl.Sync()
	}

	// xmlSplitter distributes adjacent records with the same identifier prefix
	xmlSplitter := func(inp <-chan Extract, out chan<- string) {

		// close channel when all records have been processed
		defer close(out)

		var (
			fl   *os.File
			wrtr *bufio.Writer
			zpr  *gzip.Writer
		)

		currTag := ""
		prevTag := ""

		// remember previous record
		prev := Extract{}

		for curr := range inp {

			// use first few characters of identifier
			currTag = IdentifierKey(curr.Ident)
			if currTag == "" {
				continue
			}

			// then truncate to 2 or 3 character prefix
			if len(currTag) > 2 {
				key := currTag[:2]
				num, ok := trieLen[key]
				if ok && num > 2 {
					currTag = currTag[:3]
				} else {
					currTag = currTag[:2]
				}
			}

			if fl == nil {
				// open initial file
				fl, wrtr, zpr = openSaver(merg, currTag, zipp)
				if wrtr == nil {
					continue
				}

				// send first opening tag and indent
				wrtr.WriteString("<InvDocumentSet>\n  ")
			}

			// compare keys from adjacent term lists
			if prev.Text != "" && prevTag != currTag {

				// after IdentifierKey converts space to underscore,
				// okay that x_ and x0 will be out of alphabetical order

				// send closing tag
				wrtr.WriteString("</InvDocumentSet>\n")

				closeSaver(fl, wrtr, zpr)

				out <- currTag

				// force garbage collection
				debug.FreeOSMemory()

				runtime.Gosched()

				// open next file
				fl, wrtr, zpr = openSaver(merg, currTag, zipp)
				if wrtr == nil {
					continue
				}

				// send opening tag and indent
				wrtr.WriteString("<InvDocumentSet>\n  ")
			}

			// send one InvDocument
			str := strings.TrimSpace(curr.Text)

			wrtr.WriteString(str)
			if !strings.HasSuffix(str, "\n") {
				wrtr.WriteString("\n")
			}

			// now remember this record
			prev = curr

			prevTag = currTag
		}

		if prev.Text != "" {

			// send last closing tag
			wrtr.WriteString("</InvDocumentSet>\n")

			closeSaver(fl, wrtr, zpr)

			out <- currTag

			// force garbage collection
			debug.FreeOSMemory()

			runtime.Gosched()
		}
	}

	// launch single splitter goroutine
	go xmlSplitter(inp, out)

	return out
}

func CreatePromoters(args []string, prom, field string) <-chan string {

	if args == nil {
		return nil
	}

	out := make(chan string, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create promoter channel\n")
		os.Exit(1)
	}

	// xmlPromoter splits inverted index groups into subdirectories for term lists and postings files
	xmlPromoter := func(wg *sync.WaitGroup, fileName string, out chan<- string) {

		defer wg.Done()

		f, err := os.Open(fileName)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open input file '%s'\n", fileName)
			os.Exit(1)
		}

		// close input file when all records have been processed
		defer f.Close()

		var in io.Reader

		in = f

		// if suffix is ".gz", use decompressor
		iszip := false
		if strings.HasSuffix(fileName, ".gz") {
			iszip = true
		}

		if iszip {
			brd := bufio.NewReader(f)
			if brd == nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create buffered reader on '%s'\n", fileName)
				os.Exit(1)
			}
			zpr, err := gzip.NewReader(brd)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create decompressor on '%s'\n", fileName)
				os.Exit(1)
			}

			// close decompressor when all records have been processed
			defer zpr.Close()

			// use decompressor for reading file
			in = zpr
		}

		rdr := CreateReader(in)

		if rdr == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create XML Block Reader\n")
			os.Exit(1)
		}

		getOnePosting := func(text string) (string, []int32, []string) {

			var data []int32
			var atts []string

			term := ""

			doPromote := func(tag, attr, content string) {

				if tag == "InvKey" {

					// term used for postings file name
					term = content

					term = strings.ToLower(term)

				} else if tag == field {

					// convert UID string to integer
					if content == "" {
						fmt.Fprintf(os.Stderr, "\nERROR: Empty UID for term '%s'\n", term)
						return
					}
					value, err := strconv.ParseInt(content, 10, 32)
					if err != nil {
						fmt.Fprintf(os.Stderr, "%s\n", err.Error())
						return
					}
					data = append(data, int32(value))

					if strings.HasPrefix(attr, "pos=\"") {
						attr = attr[5:]
						lgth := len(attr)
						if lgth > 1 && attr[lgth-1] == '"' {
							// "
							attr = attr[:lgth-1]
						}
						atts = append(atts, attr)
					}
				}
			}

			// explore data fields
			StreamValues(text[:], "InvDocument", doPromote)

			if term == "" || len(data) < 1 {
				return "", nil, nil
			}

			return term, data, atts
		}

		var (
			termPos int32
			postPos int32
			ofstPos int32

			indxList bytes.Buffer
			termList bytes.Buffer
			postList bytes.Buffer
			uqidList bytes.Buffer
			ofstList bytes.Buffer
		)

		retlength := len("\n")

		addOnePosting := func(term string, data []int32, atts []string) {

			tlength := len(term)
			dlength := len(data)
			alength := len(atts)

			// write to term list buffer
			termList.WriteString(term[:])
			termList.WriteString("\n")

			// write to postings buffer
			binary.Write(&postList, binary.LittleEndian, data)

			// write to master index buffer
			binary.Write(&indxList, binary.LittleEndian, termPos)
			binary.Write(&indxList, binary.LittleEndian, postPos)

			postPos += int32(dlength * 4)
			termPos += int32(tlength + retlength)

			// return if no position attributes
			if alength < 1 {
				return
			}
			if dlength != alength {
				fmt.Fprintf(os.Stderr, "dlength %d, alength %d\n", dlength, alength)
				return
			}

			// write term offset list for each UID
			for _, attr := range atts {

				binary.Write(&uqidList, binary.LittleEndian, ofstPos)

				atrs := strings.Split(attr, ",")
				atln := len(atrs)
				for _, att := range atrs {
					if att == "" {
						continue
					}
					value, err := strconv.ParseInt(att, 10, 32)
					if err != nil {
						fmt.Fprintf(os.Stderr, "%s\n", err.Error())
						return
					}
					binary.Write(&ofstList, binary.LittleEndian, int16(value))
				}

				ofstPos += int32(atln * 2)
			}
		}

		topOffMaster := func() {

			// phantom term and postings positions eliminates special case calculation at end
			binary.Write(&indxList, binary.LittleEndian, termPos)
			binary.Write(&indxList, binary.LittleEndian, postPos)
			binary.Write(&uqidList, binary.LittleEndian, ofstPos)
		}

		writeFile := func(dpath, fname string, bfr bytes.Buffer) {

			fpath := path.Join(dpath, fname)
			if fpath == "" {
				return
			}

			// overwrites and truncates existing file
			fl, err := os.Create(fpath)
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
				return
			}

			data := bfr.Bytes()

			wrtr := bufio.NewWriter(fl)

			_, err = wrtr.Write(data)
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			}

			wrtr.Flush()

			// fl.Sync()

			fl.Close()
		}

		writeFiveFiles := func(key string) {

			var arry [516]rune
			dpath, key := PostingPath(prom, field, key, arry)
			if dpath == "" {
				return
			}

			// make subdirectories, if necessary
			err = os.MkdirAll(dpath, os.ModePerm)
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
				return
			}

			writeFile(dpath, key+"."+field+".trm", termList)

			writeFile(dpath, key+"."+field+".pst", postList)

			writeFile(dpath, key+"."+field+".mst", indxList)

			// do not write position index and offset data files for fields with no position attributes recorded
			if uqidList.Len() > 0 && ofstList.Len() > 0 {

				writeFile(dpath, key+"."+field+".uqi", uqidList)

				writeFile(dpath, key+"."+field+".ofs", ofstList)
			}
		}

		currTag := ""
		prevTag := ""

		ok := false

		PartitionPattern("InvDocument", "", rdr,
			func(str string) {

				term, data, atts := getOnePosting(str)

				if term == "" || data == nil {
					return
				}

				ok = true

				// use first few characters of identifier
				currTag = IdentifierKey(term)

				if prevTag != currTag {

					// after IdentifierKey converts space to underscore,
					// okay that xxx_ and xxx0 will be out of alphabetical order

					// directory prefix changed from last posting
					if prevTag != "" {

						topOffMaster()
						writeFiveFiles(prevTag)
						out <- prevTag
					}

					// reset buffers and position counters
					termPos = 0
					postPos = 0
					ofstPos = 0

					indxList.Reset()
					termList.Reset()
					postList.Reset()
					uqidList.Reset()
					ofstList.Reset()
				}

				addOnePosting(term, data, atts)

				prevTag = currTag
			})

		if ok {

			// write last set of files
			topOffMaster()
			writeFiveFiles(prevTag)
			out <- prevTag
		}
	}

	var wg sync.WaitGroup

	// launch multiple promoter goroutines
	for _, str := range args {
		wg.Add(1)
		go xmlPromoter(&wg, str, out)
	}

	// launch separate anonymous goroutine to wait until all promoters are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

func CreateMatchers(phrs string, exclude bool, inp <-chan Extract) <-chan Extract {

	if inp == nil {
		return nil
	}

	if phrs == "" {
		// if -phrase (-require, -exclude) argument is not present, simply return input channel
		return inp
	}

	phrs = PrepareQuery(phrs)

	phrs = ProcessStopWords(phrs)

	clauses := PartitionQuery(phrs)

	clauses = SetFieldQualifiers(clauses, false)

	if clauses == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to parse phrase\n")
		os.Exit(1)
	}

	out := make(chan Extract, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create phrase matcher channel\n")
		os.Exit(1)
	}

	// split at punctuation, but leave < and > in to delimit content strings
	cleanupRecord := func(str string) string {

		if IsNotASCII(str) {
			str = DoAccentTransform(str)
			if HasUnicodeMarkup(str) {
				str = RepairUnicodeMarkup(str, SPACE)
			}
		}

		if HasBadSpace(str) {
			str = CleanupBadSpaces(str)
		}

		if HasAmpOrNotASCII(str) {
			str = html.UnescapeString(str)
		}

		if IsNotASCII(str) {
			if HasGreek(str) {
				str = SpellGreek(str)
			}
		}

		if HasHyphenOrApostrophe(str) {
			str = FixSpecialCases(str)
		}

		var buffer strings.Builder

		for _, ch := range str {
			if ch > 127 {
				buffer.WriteRune(' ')
			} else if unicode.IsLetter(ch) || unicode.IsDigit(ch) {
				buffer.WriteRune(ch)
			} else if ch == '<' || ch == '>' {
				buffer.WriteRune(' ')
				buffer.WriteRune(ch)
				buffer.WriteRune(' ')
			} else {
				buffer.WriteRune(' ')
			}
		}

		res := buffer.String()

		res = strings.TrimSpace(res)
		res = CompressRunsOfSpaces(res)
		res = strings.ToLower(res)

		var chain []string

		terms := strings.Fields(res)

		// replace unwanted and stop words with plus sign
		for _, item := range terms {

			// allow tilde proximity indicator
			if item == "~" {
				chain = append(chain, item)
				continue
			}

			// skip if stop word, breaking word pair chain
			if DeStop && IsStopWord(item) {
				chain = append(chain, "+")
				continue
			}

			// apply stemming algorithm
			if DoStem {
				isWildCard := strings.HasSuffix(item, "*")
				if isWildCard {
					// temporarily remove trailing asterisk
					item = strings.TrimSuffix(item, "*")
				}
				item = porter2.Stem(item)
				item = strings.TrimSpace(item)
				if isWildCard {
					// do wildcard search in stemmed term list
					item += "*"
				}
			}

			// record single term
			chain = append(chain, item)
		}

		// rejoin into processed sentence
		tmp := strings.Join(chain, " ")

		tmp = CompressRunsOfSpaces(tmp)
		tmp = strings.TrimSpace(tmp)

		return tmp
	}

	proximitySearch := func(srch, str string) bool {

		// split into two words separated by run of tildes
		words := strings.Split(str, " ")
		// proximity variables
		first := ""
		secnd := ""
		dist := 0
		for _, item := range words {
			if strings.Contains(item, "~") {
				dist = len(item)
			} else if first == "" {
				first = item
			} else if secnd == "" {
				secnd = item
			} else {
				fmt.Fprintf(os.Stderr, "\nERROR: More than two terms in proximity search\n")
				os.Exit(1)
			}
		}
		if first == "" || secnd == "" || dist < 1 {
			fmt.Fprintf(os.Stderr, "\nERROR: Fields missing in proximity search\n")
			os.Exit(1)
		}

		terms := strings.Fields(srch)

		for j, fst := range terms {
			if fst != first {
				continue
			}
			rest := terms[j+1:]
			for k, sec := range rest {
				if sec == secnd {
					return true
				}
				if k >= dist {
					break
				}
			}
		}

		return false
	}

	// check each phrase against record
	testPhrase := func(srch string, tokens []string) bool {

		eval := func(str string) bool {

			if strings.Contains(str, "~") {
				return proximitySearch(srch, str)
			}

			return strings.Contains(srch, str)
		}

		nextToken := func() string {

			if len(tokens) < 1 {
				return ""
			}

			tkn := tokens[0]
			tokens = tokens[1:]

			return tkn
		}

		// recursive definitions
		var excl func() (bool, string)
		var expr func() (bool, string)
		var fact func() (bool, string)
		var term func() (bool, string)

		fact = func() (bool, string) {

			var res bool

			tkn := nextToken()
			if tkn == "(" {
				res, tkn = expr()
				if tkn == ")" {
					tkn = nextToken()
				}
			} else {
				res = eval(tkn)
				tkn = nextToken()
			}

			return res, tkn
		}

		excl = func() (bool, string) {

			var val bool

			res, tkn := fact()
			for tkn == "!" {
				val, tkn = fact()
				if val {
					res = false
				}
			}

			return res, tkn
		}

		term = func() (bool, string) {

			var val bool

			res, tkn := excl()
			for tkn == "&" {
				val, tkn = excl()
				if !val {
					res = false
				}
			}

			return res, tkn
		}

		expr = func() (bool, string) {

			var val bool

			res, tkn := term()
			for tkn == "|" {
				val, tkn = term()
				if val {
					res = true
				}
			}

			return res, tkn
		}

		// enter recursive descent parser
		found, _ := expr()

		return found
	}

	// xmlMatcher reads partitioned XML from channel and removes records that do not contain the phrase(s)
	xmlMatcher := func(wg *sync.WaitGroup, inp <-chan Extract, out chan<- Extract) {

		// report when this phrase matcher has no more records to process
		defer wg.Done()

		// read partitioned XML from producer channel
		for ext := range inp {

			idx := ext.Index
			text := ext.Text

			if text == "" {
				// should never see empty input data
				out <- Extract{idx, "", text, nil}
				continue
			}

			srch := cleanupRecord(text)

			ok := testPhrase(srch, clauses)

			if exclude != ok {
				// send text of record if phrase match succeeded with -require, or failed with -exclude
				out <- Extract{idx, "", text, nil}
				continue
			}

			// otherwise send empty text so unshuffler does not have to deal with record index gaps
			out <- Extract{idx, "", "", nil}
		}
	}

	var wg sync.WaitGroup

	// launch multiple phrase matcher goroutines
	for i := 0; i < NumServe; i++ {
		wg.Add(1)
		go xmlMatcher(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all phrase matchers are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

func CreateExternalIndexer(args []string, zipp bool, in io.Reader) int {

	recordCount := 0

	transform := make(map[string]string)

	readTransformTable := func(tf string) {

		inFile, err := os.Open(tf)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Unable to open transformation file %s - %s\n", tf, err.Error())
			os.Exit(1)
		}

		scant := bufio.NewScanner(inFile)

		// populate transformation map
		for scant.Scan() {

			line := scant.Text()
			frst, scnd := SplitInTwoAt(line, "\t", LEFT)

			transform[frst] = scnd
		}

		inFile.Close()
	}

	// BIOCONCEPTS INDEXER

	// create intermediate table for {chemical|disease|gene}2pubtatorcentral.gz (undocumented)
	if args[0] == "-bioconcepts" {

		if len(args) < 2 {
			fmt.Fprintf(os.Stderr, "\nERROR: Insufficient arguments for -bioconcepts\n")
			os.Exit(1)
		}

		// read transformation file
		tf := args[1]
		readTransformTable(tf)

		var buffer strings.Builder
		count := 0
		okay := false

		wrtr := bufio.NewWriter(os.Stdout)

		scanr := bufio.NewScanner(in)

		currpmid := ""

		// read lines of PMIDs and extracted concepts
		for scanr.Scan() {

			line := scanr.Text()

			cols := strings.Split(line, "\t")
			if len(cols) != 5 {
				continue
			}

			pmid := cols[0]
			if currpmid != pmid {
				// end current block
				currpmid = pmid

				if pmid == "" {
					continue
				}

				recordCount++
				count++

				if count >= 1000 {
					count = 0
					txt := buffer.String()
					if txt != "" {
						// print current buffer
						wrtr.WriteString(txt[:])
					}
					buffer.Reset()
				}

				okay = true
			}

			typ := cols[1]
			val := cols[2]
			switch typ {
			case "Gene":
				genes := strings.Split(val, ";")
				for _, gene := range genes {
					if gene == "None" {
						continue
					}
					buffer.WriteString(pmid)
					buffer.WriteString("\t")
					buffer.WriteString("GENE")
					buffer.WriteString("\t")
					buffer.WriteString(gene)
					buffer.WriteString("\n")
					gn, ok := transform[gene]
					if !ok || gn == "" {
						continue
					}
					buffer.WriteString(pmid)
					buffer.WriteString("\t")
					buffer.WriteString("GENE")
					buffer.WriteString("\t")
					buffer.WriteString(gn)
					buffer.WriteString("\n")
				}
			case "Disease":
				if strings.HasPrefix(val, "MESH:") {
					diszs := strings.Split(val[5:], "|")
					for _, disz := range diszs {
						buffer.WriteString(pmid)
						buffer.WriteString("\t")
						buffer.WriteString("DISZ")
						buffer.WriteString("\t")
						buffer.WriteString(disz)
						buffer.WriteString("\n")
						dn, ok := transform[disz]
						if !ok || dn == "" {
							continue
						}
						buffer.WriteString(pmid)
						buffer.WriteString("\t")
						buffer.WriteString("DISZ")
						buffer.WriteString("\t")
						buffer.WriteString(dn)
						buffer.WriteString("\n")
					}
				} else if strings.HasPrefix(val, "OMIM:") {
					omims := strings.Split(val[5:], "|")
					for _, omim := range omims {
						// was OMIM, now fused with DISZ
						buffer.WriteString(pmid)
						buffer.WriteString("\t")
						buffer.WriteString("DISZ")
						buffer.WriteString("\t")
						buffer.WriteString(omim)
						buffer.WriteString("\n")
					}
				}
			case "Chemical":
				if strings.HasPrefix(val, "MESH:") {
					chems := strings.Split(val[5:], "|")
					for _, chem := range chems {
						buffer.WriteString(pmid)
						buffer.WriteString("\t")
						buffer.WriteString("CHEM")
						buffer.WriteString("\t")
						buffer.WriteString(chem)
						buffer.WriteString("\n")
						ch, ok := transform[chem]
						if !ok || ch == "" {
							continue
						}
						buffer.WriteString(pmid)
						buffer.WriteString("\t")
						buffer.WriteString("CHEM")
						buffer.WriteString("\t")
						buffer.WriteString(ch)
						buffer.WriteString("\n")
					}
				} else if strings.HasPrefix(val, "CHEBI:") {
					// was CEBI, now fused with CHEM
					buffer.WriteString(pmid)
					buffer.WriteString("\t")
					buffer.WriteString("CHEM")
					buffer.WriteString("\t")
					buffer.WriteString(val[6:])
					buffer.WriteString("\n")
				}
			case "Species":
			case "Mutation":
			case "CellLine":
			default:
			}
		}

		if okay {
			txt := buffer.String()
			if txt != "" {
				// print current buffer
				wrtr.WriteString(txt[:])
			}
		}
		buffer.Reset()

		wrtr.Flush()

		return recordCount
	}

	// GENERIF INDEXER

	// create intermediate table for generifs_basic.gz (undocumented)
	if args[0] == "-generif" || args[0] == "-generifs" {

		if len(args) < 2 {
			fmt.Fprintf(os.Stderr, "\nERROR: Insufficient arguments for -generif\n")
			os.Exit(1)
		}

		// read transformation file
		tf := args[1]
		readTransformTable(tf)

		var buffer strings.Builder
		count := 0
		okay := false

		wrtr := bufio.NewWriter(os.Stdout)

		scanr := bufio.NewScanner(in)

		currpmid := ""

		// skip first line with column heading names
		for scanr.Scan() {

			line := scanr.Text()
			cols := strings.Split(line, "\t")
			if len(cols) != 5 {
				fmt.Fprintf(os.Stderr, "Unexpected number of columns (%d) in generifs_basic.gz\n", len(cols))
				os.Exit(1)
			}
			if len(cols) != 5 || cols[0] != "#Tax ID" {
				fmt.Fprintf(os.Stderr, "Unrecognized contents in generifs_basic.gz\n")
				os.Exit(1)
			}
			break
		}

		// read lines of PMIDs and gene references
		for scanr.Scan() {

			line := scanr.Text()

			cols := strings.Split(line, "\t")
			if len(cols) != 5 {
				continue
			}

			val := cols[2]
			pmids := strings.Split(val, ",")
			for _, pmid := range pmids {
				if currpmid != pmid {
					// end current block
					currpmid = pmid

					if pmid == "" {
						continue
					}

					recordCount++
					count++

					if count >= 1000 {
						count = 0
						txt := buffer.String()
						if txt != "" {
							// print current buffer
							wrtr.WriteString(txt[:])
						}
						buffer.Reset()
					}

					okay = true
				}

				gene := cols[1]
				// was GRIF, now fused with GENE
				buffer.WriteString(pmid)
				buffer.WriteString("\t")
				buffer.WriteString("GENE")
				buffer.WriteString("\t")
				buffer.WriteString(gene)
				buffer.WriteString("\n")
				gn, ok := transform[gene]
				if !ok || gn == "" {
					continue
				}
				buffer.WriteString(pmid)
				buffer.WriteString("\t")
				buffer.WriteString("GENE")
				buffer.WriteString("\t")
				buffer.WriteString(gn)
				buffer.WriteString("\n")
			}
		}

		if okay {
			txt := buffer.String()
			if txt != "" {
				// print current buffer
				wrtr.WriteString(txt[:])
			}
		}
		buffer.Reset()

		wrtr.Flush()

		return recordCount
	}

	// THEME INDEXER

	// create intermediate table for chemical-gene-disease themes (undocumented)
	if args[0] == "-theme" || args[0] == "-themes" {

		if len(args) < 4 {
			fmt.Fprintf(os.Stderr, "\nERROR: Insufficient arguments for -theme\n")
			os.Exit(1)
		}

		one := args[1]
		two := args[2]
		tag := args[3]

		// for disambiguating B, E, E+, and J themes, in CHDI, CHGE, GEDI, and GEGE data sets
		sfx := ""
		if len(tag) > 0 {
			switch tag[0] {
			case 'C':
				sfx = "c"
			case 'G':
				sfx = "g"
			}
		}

		fl, err := os.Open(one)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open input file '%s'\n", one)
			os.Exit(1)
		}

		scanr := bufio.NewScanner(fl)

		var columns []string
		numCols := 0

		// read first line with column heading names
		if scanr.Scan() {

			line := scanr.Text()
			line = strings.Replace(line, "+", "plus", -1)
			line = strings.Replace(line, "-", "minus", -1)
			columns = strings.Split(line, "\t")
			numCols = len(columns)

			if numCols < 3 {
				fmt.Fprintf(os.Stderr, "Unexpected number of columns (%d) in part-i file\n", numCols)
				os.Exit(1)
			}
			if columns[0] != "path" {
				fmt.Fprintf(os.Stderr, "Unrecognized contents in part-i file\n")
				os.Exit(1)
			}
			// increment by 2 to ignore flagship indicator fields
			for i := 1; i < numCols; i += 2 {
				// disambiguate B, E, E+, and J themes that appear in two data sets
				theme := columns[i]
				switch theme {
				case "B":
					columns[i] = "B,B" + sfx
				case "E":
					columns[i] = "E,E" + sfx
				case "Eplus":
					columns[i] = "Eplus,E" + sfx + "plus"
				case "Eminus":
					columns[i] = "Eminus,E" + sfx + "minus"
				case "J":
					columns[i] = "J,J" + sfx
				}
			}
		}

		var scores []int

		for i := 0; i < numCols; i++ {
			scores = append(scores, 0)
		}

		mapper := make(map[string]string)

		// read lines of dependency paths, scores for each theme
		for scanr.Scan() {

			line := scanr.Text()

			cols := strings.Split(line, "\t")
			if len(cols) != numCols {
				fmt.Fprintf(os.Stderr, "Mismatched columns in '%s'\n", line)
				continue
			}

			for i := 0; i < numCols; i++ {
				scores[i] = 0
			}

			sum := 0
			for i := 1; i < numCols; i += 2 {
				str := cols[i]
				str, _ = SplitInTwoAt(str, ".", LEFT)
				val, err := strconv.Atoi(str)
				if err != nil {
					fmt.Fprintf(os.Stderr, "Unrecognized number '%s'\n", str)
					continue
				}
				scores[i] = val
				sum += val
			}
			if sum == 0 {
				continue
			}

			path := cols[0]
			path = strings.ToLower(path)
			themes := ""
			comma := ""
			for i := 1; i < numCols; i += 2 {
				// find scores over cutoff
				if scores[i]*3 > sum {
					themes += comma
					themes += columns[i]
					comma = ","
				}
			}
			if themes == "" {
				continue
			}
			// populate theme lookup table
			mapper[path] = themes
		}

		fl.Close()

		fl, err = os.Open(two)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open input file '%s'\n", two)
			os.Exit(1)
		}

		var buffer strings.Builder
		count := 0
		okay := false

		wrtr := bufio.NewWriter(os.Stdout)

		scanr = bufio.NewScanner(fl)

		// read lines of PMIDs and dependency paths
		for scanr.Scan() {

			line := scanr.Text()

			cols := strings.Split(line, "\t")
			if len(cols) != 14 {
				fmt.Fprintf(os.Stderr, "Mismatched columns in '%s'\n", line)
				continue
			}

			pmid := cols[0]
			path := cols[12]
			path = strings.ToLower(path)
			themes, ok := mapper[path]
			if !ok {
				continue
			}

			thms := strings.Split(themes, ",")
			for _, theme := range thms {
				if theme == "" {
					continue
				}
				buffer.WriteString(pmid)
				buffer.WriteString("\t")
				buffer.WriteString("THME")
				buffer.WriteString("\t")
				buffer.WriteString(theme)
				buffer.WriteString("\n")
			}

			recordCount++
			count++

			if count >= 1000 {
				count = 0
				txt := buffer.String()
				if txt != "" {
					// print current buffer
					wrtr.WriteString(txt[:])
				}
				buffer.Reset()
			}

			okay = true
		}

		if okay {
			txt := buffer.String()
			if txt != "" {
				// print current buffer
				wrtr.WriteString(txt[:])
			}
		}
		buffer.Reset()

		wrtr.Flush()

		fl.Close()

		return recordCount
	}

	// DEPENDENCY PATH INDEXER

	// create intermediate table for chemical-gene-disease dependency paths (undocumented)
	if args[0] == "-dpath" || args[0] == "-dpaths" {

		var buffer strings.Builder
		count := 0
		okay := false

		replr := strings.NewReplacer(
			">", "_gtrthan_",
			"<", "_lssthan_",
			"/", "_slash_",
			"%", "_prcnt_",
			":", "_colln_",
			"+", "_pluss_",
			"!", "_exclam_",
			"?", "_qmark_",
			"'", "_squot_",
			"(", "_lparen_",
			")", "_rparen_",
		)
		if replr == nil {
			fmt.Fprintf(os.Stderr, "Unable to create replacer\n")
			os.Exit(1)
		}

		wrtr := bufio.NewWriter(os.Stdout)

		scanr := bufio.NewScanner(in)

		// read lines of PMIDs and dependency paths
		for scanr.Scan() {

			line := scanr.Text()

			cols := strings.Split(line, "\t")
			if len(cols) != 14 {
				fmt.Fprintf(os.Stderr, "Mismatched columns in '%s'\n", line)
				continue
			}

			pmid := cols[0]
			path := cols[12]
			path = strings.ToLower(path)

			// rescue known characters
			tmp := CompressRunsOfSpaces(path)
			tmp = strings.TrimSpace(tmp)

			tmp = " " + tmp + " "

			tmp = replr.Replace(tmp)

			tmp = CompressRunsOfSpaces(tmp)
			tmp = strings.TrimSpace(tmp)

			// final cleanup
			tmp = strings.Replace(tmp, "|", "_", -1)
			tmp = strings.Replace(tmp, "__", "_", -1)

			pths := strings.Split(tmp, " ")
			for _, pth := range pths {
				if pth == "" {
					continue
				}
				buffer.WriteString(pmid)
				buffer.WriteString("\t")
				buffer.WriteString("PATH")
				buffer.WriteString("\t")
				buffer.WriteString(pth)
				buffer.WriteString("\n")
			}

			recordCount++
			count++

			if count >= 1000 {
				count = 0
				txt := buffer.String()
				if txt != "" {
					// print current buffer
					wrtr.WriteString(txt[:])
				}
				buffer.Reset()
			}

			okay = true
		}

		if okay {
			txt := buffer.String()
			if txt != "" {
				// print current buffer
				wrtr.WriteString(txt[:])
			}
		}
		buffer.Reset()

		wrtr.Flush()

		return recordCount
	}

	// THESIS INDEXER

	// create -e2index file for bioconcepts, geneRIFs, and themes and their dependency paths (undocumented)
	if args[0] == "-thesis" {

		// e.g., -thesis 250000 "$target" "biocchem"
		if len(args) < 4 {
			fmt.Fprintf(os.Stderr, "\nERROR: Insufficient arguments for -thesis\n")
			os.Exit(1)
		}

		chunk, err := strconv.Atoi(args[1])
		if err != nil {
			fmt.Fprintf(os.Stderr, "Unrecognized count - '%s'\n", err.Error())
			os.Exit(1)
		}
		target := strings.TrimSuffix(args[2], "/")
		prefix := args[3]

		suffix := "e2x"
		sfx := suffix
		if zipp {
			sfx += ".gz"
		}

		fnum := 0

		fld := ""

		scanr := bufio.NewScanner(os.Stdin)

		processChunk := func() bool {

			// map for combined index
			indexed := make(map[string][]string)

			writeChunk := func() {

				var (
					fl   *os.File
					wrtr *bufio.Writer
					zpr  *gzip.Writer
					err  error
				)

				fnum++
				fpath := fmt.Sprintf("%s/%s%03d.%s", target, prefix, fnum, sfx)
				fl, err = os.Create(fpath)
				if err != nil {
					fmt.Fprintf(os.Stderr, "%s\n", err.Error())
					return
				}
				defer fl.Close()

				pth := fmt.Sprintf("%s%03d.%s", prefix, fnum, suffix)
				os.Stderr.WriteString(pth + "\n")

				var out io.Writer

				out = fl

				if zipp {

					zpr, err = gzip.NewWriterLevel(fl, gzip.BestSpeed)
					if err != nil {
						fmt.Fprintf(os.Stderr, "%s\n", err.Error())
						return
					}

					out = zpr
				}

				wrtr = bufio.NewWriter(out)
				if wrtr == nil {
					fmt.Fprintf(os.Stderr, "Unable to create bufio.NewWriter\n")
					return
				}

				var buffer strings.Builder
				count := 0

				buffer.WriteString("<IdxDocumentSet>\n")

				// sort fields in alphabetical order
				var keys []string
				for ky := range indexed {
					keys = append(keys, ky)
				}

				if len(keys) > 1 {
					sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })
				}

				for _, idx := range keys {

					item, ok := indexed[idx]
					if !ok {
						continue
					}

					uid := item[0]
					data := item[1:]

					if uid == "" || len(data) < 1 {
						continue
					}

					if len(data) > 1 {
						sort.Slice(data, func(i, j int) bool { return data[i] < data[j] })
					}

					buffer.WriteString("  <IdxDocument>\n")
					buffer.WriteString("    <IdxUid>")
					buffer.WriteString(uid)
					buffer.WriteString("</IdxUid>\n")
					buffer.WriteString("    <IdxSearchFields>\n")

					prev := ""
					for len(data) > 0 {
						val := data[0]
						data = data[1:]

						if val == prev {
							continue
						}

						buffer.WriteString("      <")
						buffer.WriteString(fld)
						buffer.WriteString(">")
						buffer.WriteString(val)
						buffer.WriteString("</")
						buffer.WriteString(fld)
						buffer.WriteString(">\n")

						prev = val
					}

					buffer.WriteString("    </IdxSearchFields>\n")
					buffer.WriteString("  </IdxDocument>\n")

					recordCount++
					count++

					if count >= 1000 {
						count = 0
						txt := buffer.String()
						if txt != "" {
							// print current buffer
							wrtr.WriteString(txt[:])
						}
						buffer.Reset()
					}
				}

				buffer.WriteString("</IdxDocumentSet>\n")

				txt := buffer.String()
				if txt != "" {
					// print current buffer
					wrtr.WriteString(txt[:])
				}
				buffer.Reset()

				wrtr.Flush()

				if zpr != nil {
					zpr.Close()
				}
			}

			lineCount := 0
			okay := false

			// read lines of dependency paths, scores for each theme
			for scanr.Scan() {

				line := scanr.Text()

				cols := strings.Split(line, "\t")
				if len(cols) != 3 {
					fmt.Fprintf(os.Stderr, "Mismatched columns in '%s'\n", line)
					continue
				}

				uid := cols[0]
				fd := cols[1]
				val := cols[2]
				if uid == "" || fd == "" || val == "" {
					continue
				}
				if fld == "" {
					fld = fd
				}
				if fld != fd {
					fmt.Fprintf(os.Stderr, "Field '%s' expected, '%s' found\n", fld, fd)
					continue
				}

				val = strings.ToLower(val)
				// convert angle brackets in chemical names
				val = html.EscapeString(val)

				data, ok := indexed[uid]
				if !ok {
					data = make([]string, 0, 2)
					// first entry on new slice is uid
					data = append(data, uid)
				}
				data = append(data, val)
				// always need to update indexed, since data may be reallocated
				indexed[uid] = data

				okay = true

				lineCount++
				if lineCount > chunk {
					break
				}
			}

			if okay {
				writeChunk()
				return true
			}

			return false
		}

		for processChunk() {
			// loop until scanner runs out of lines
		}

		return recordCount
	}

	return 0
}

func CreateExternalArchive(stash string, args []string) <-chan string {

	createPresenters := func(args []string) []<-chan Plex {

		if args == nil {
			return nil
		}

		numFiles := len(args)
		if numFiles < 1 {
			fmt.Fprintf(os.Stderr, "\nERROR: Not enough indexed files to merge\n")
			os.Exit(1)
		}

		chns := make([]<-chan Plex, numFiles)
		if chns == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create presenter channel array\n")
			os.Exit(1)
		}

		// xmlPresenter sends partitioned XML strings through channel
		xmlPresenter := func(fileNum int, fileName string, out chan<- Plex) {

			// close channel when all records have been processed
			defer close(out)

			f, err := os.Open(fileName)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to open input file '%s'\n", fileName)
				os.Exit(1)
			}

			// close input file when all records have been processed
			defer f.Close()

			var in io.Reader

			in = f

			// if suffix is ".gz", use decompressor
			iszip := false
			if strings.HasSuffix(fileName, ".gz") {
				iszip = true
			}

			if iszip {
				brd := bufio.NewReader(f)
				if brd == nil {
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to create buffered reader on '%s'\n", fileName)
					os.Exit(1)
				}
				zpr, err := gzip.NewReader(brd)
				if err != nil {
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to create decompressor on '%s'\n", fileName)
					os.Exit(1)
				}

				// close decompressor when all records have been processed
				defer zpr.Close()

				// use decompressor for reading file
				in = zpr
			}

			rdr := CreateReader(in)

			if rdr == nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create XML Block Reader\n")
				os.Exit(1)
			}

			find := ParseIndex("IdxUid")

			// partition all input by pattern and send XML substring through channel
			PartitionPattern("IdxDocument", "", rdr,
				func(str string) {
					id := FindIdentifier(str[:], "IdxDocument", find)

					out <- Plex{fileNum, id, str, 0, nil}
				})
		}

		// launch multiple presenter goroutines
		for i, str := range args {

			chn := make(chan Plex, ChanDepth)
			if chn == nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create presenter channel\n")
				os.Exit(1)
			}

			go xmlPresenter(i, str, chn)

			chns[i] = chn
		}

		// no need for separate anonymous goroutine to wait until all presenters are done

		return chns
	}

	createManifold := func(inp []<-chan Plex) <-chan Plex {

		if inp == nil {
			return nil
		}

		out := make(chan Plex, ChanDepth)
		if out == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create manifold channel\n")
			os.Exit(1)
		}

		// xmlManifold restores alphabetical order of merged postings
		xmlManifold := func(inp []<-chan Plex, out chan<- Plex) {

			// close channel when all records have been processed
			defer close(out)

			// initialize empty heap
			hp := &PlexHeap{}
			heap.Init(hp)

			// read first object from all input channels in turn
			for _, chn := range inp {
				plx, ok := <-chn
				if ok {
					heap.Push(hp, plx)
				}
			}

			// array to collect strings with same identifier
			var arry []string

			prevIdent := ""
			rec := 0

			// reading from heap returns objects in alphabetical order
			for hp.Len() > 0 {

				// remove lowest item from heap, use interface type assertion
				curr := heap.Pop(hp).(Plex)

				// compare adjacent record identifiers
				if prevIdent == curr.Ident {

					// save next indexed object string in slice
					arry = append(arry, curr.Text)

				} else {

					if len(arry) > 0 {

						rec++
						// send set from previous identifier to output channel
						out <- Plex{0, prevIdent, "", rec, arry}

						// empty the slice
						arry = nil
					}

					// remember new identifier
					prevIdent = curr.Ident

					// save first indexed object with this identifier
					arry = append(arry, curr.Text)
				}

				// read next object from channel that just supplied lowest item
				chn := inp[curr.Which]
				plx, ok := <-chn
				if ok {
					heap.Push(hp, plx)
				}
			}

			if len(arry) > 0 {

				rec++
				// send last record
				out <- Plex{0, prevIdent, "", rec, arry}

				arry = nil
			}
		}

		// launch single manifold goroutine
		go xmlManifold(inp, out)

		return out
	}

	createMergers := func(inp <-chan Plex) <-chan Extract {

		if inp == nil {
			return nil
		}

		out := make(chan Extract, ChanDepth)
		if out == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create merger channel\n")
			os.Exit(1)
		}

		// xmlMerger fuses adjacent IdxDocument records with the same identifier
		xmlMerger := func(wg *sync.WaitGroup, inp <-chan Plex, out chan<- Extract) {

			defer wg.Done()

			var buffer strings.Builder

			fusePostings := func(key string, data []string) string {

				fields := make(map[string]map[string]string)

				addIdents := func(fld, pos, uid string) {

					// no need for mutex here
					positions, ok := fields[fld]
					if !ok {
						positions = make(map[string]string)
						fields[fld] = positions
					}

					positions[uid] = pos
				}

				addUID := func(tag, attr, content string) {

					if tag != "IdxUid" {

						addIdents(tag, attr, content)
					}
				}

				for _, str := range data {
					StreamValues(str[:], "IdxDocument", addUID)
				}

				buffer.Reset()

				buffer.WriteString("<IdxDocument>\n")
				buffer.WriteString("<IdxUid>")
				buffer.WriteString(key)
				buffer.WriteString("</IdxUid>\n")
				buffer.WriteString("<InvIDs>\n")

				// sort fields in alphabetical order
				var keys []string
				for ky := range fields {
					keys = append(keys, ky)
				}
				sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })

				for _, fld := range keys {

					positions := fields[fld]

					var arry []string

					for item := range positions {
						arry = append(arry, item)
					}

					if len(arry) > 1 {
						sort.Slice(arry, func(i, j int) bool {
							// numeric sort on strings checks lengths first
							lni := len(arry[i])
							lnj := len(arry[j])
							// shorter string is numerically less, assuming no leading zeros
							if lni < lnj {
								return true
							}
							if lni > lnj {
								return false
							}
							// same length, can now do string comparison on contents
							return arry[i] < arry[j]
						})
					}

					// print list of UIDs, skipping duplicates
					last := ""
					for _, uid := range arry {
						// detect duplicate UIDs, now in same list after conversion of one term entry from foreign alphabet
						if uid == last {
							continue
						}
						buffer.WriteString("<")
						buffer.WriteString(fld)
						atr := positions[uid]
						if atr != "" {
							buffer.WriteString(" ")
							buffer.WriteString(atr)
						}
						buffer.WriteString(">")
						buffer.WriteString(uid)
						buffer.WriteString("</")
						buffer.WriteString(fld)
						buffer.WriteString(">\n")

						last = uid
					}
				}

				buffer.WriteString("</InvIDs>\n")
				buffer.WriteString("</IdxDocument>\n")

				txt := buffer.String()

				return txt
			}

			for plx := range inp {

				rec := plx.Index
				key := plx.Ident
				data := plx.Sibs

				if len(data) < 1 {
					continue
				}

				str := fusePostings(key, data)

				out <- Extract{rec, key, str, nil}

				runtime.Gosched()
			}
		}

		var wg sync.WaitGroup

		// launch multiple merger goroutines
		for i := 0; i < NumServe; i++ {
			wg.Add(1)
			go xmlMerger(&wg, inp, out)
		}

		// launch separate anonymous goroutine to wait until all mergers are done
		go func() {
			wg.Wait()
			close(out)
		}()

		return out
	}

	chns := createPresenters(args)
	mfld := createManifold(chns)
	mrgr := createMergers(mfld)
	stsq := CreateStashers(stash, "IdxDocument", "IdxDocument/IdxUid", ".e2x", false, true, 50000, mrgr)

	if chns == nil || mfld == nil || mrgr == nil || stsq == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create extra index stasher\n")
		os.Exit(1)
	}

	return stsq
}

// MAIN FUNCTION

func main() {

	// skip past executable name
	args := os.Args[1:]

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: No command-line arguments supplied to rchive\n")
		os.Exit(1)
	}

	// CONCURRENCY, CLEANUP, AND DEBUGGING FLAGS

	// do these first because -defcpu and -maxcpu can be sent from wrapper before other arguments

	ncpu := runtime.NumCPU()
	if ncpu < 1 {
		ncpu = 1
	}

	// wrapper can limit maximum number of processors to use (undocumented)
	maxProcs := ncpu
	defProcs := 0

	// concurrent performance tuning parameters, can be overridden by -proc and -cons
	numProcs := 0
	serverRatio := 4

	// number of servers usually calculated by -cons server ratio, but can be overridden by -serv
	NumServe = 0

	// number of channels usually equals number of servers, but can be overridden by -chan
	ChanDepth = 0

	// miscellaneous tuning parameters
	HeapSize = 16
	FarmSize = 64

	// garbage collector control can be set by environment variable or default value with -gogc 0
	goGc := 200
	gcdefault := true

	// -flag sets -strict or -mixed cleanup flags from argument
	flgs := ""

	DeStop = true

	// read data from file instead of stdin
	fileName := ""

	// debugging
	dbug := false
	stts := false
	timr := false

	// profiling
	prfl := false

	// element to use as local data index
	indx := ""

	// file of index values for removing duplicates
	unqe := ""

	// path for local data indexed as trie
	stsh := ""
	ftch := ""
	strm := ""

	// path for local extra link data
	smmn := ""

	// flag for inverted index
	nvrt := false

	// flag for combining sets of inverted files
	join := false

	// flag for combining sets of inverted files
	fuse := false

	// destination directory for merging and splitting inverted files
	merg := ""

	// base destination directory for promoting inverted index to retrieval indices
	prom := ""

	// field for promoting inverted index files
	fild := ""

	// base for queries
	base := ""

	// query by phrase, normalized terms (with truncation wildcarding)
	phrs := ""
	rlxd := false
	xact := false
	mock := false
	btch := false

	// print term list with counts
	trms := ""
	plrl := false
	psns := false

	// use gzip compression on local data files
	zipp := false

	// print UIDs and hash values
	hshv := false

	// convert UIDs to archive trie
	trei := false

	// compare input record against stash
	cmpr := false
	cmprType := ""
	ignr := ""

	// flag missing identifiers
	msng := false

	// flag records with damaged embedded HTML tags
	dmgd := false
	dmgdType := ""

	// kludge to use non-threaded fetching for windows
	windows := false

	// get numeric value
	getNumericArg := func(name string, zer, min, max int) int {

		if len(args) < 2 {
			fmt.Fprintf(os.Stderr, "\nERROR: %s is missing\n", name)
			os.Exit(1)
		}
		value, err := strconv.Atoi(args[1])
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: %s (%s) is not an integer\n", name, args[1])
			os.Exit(1)
		}
		// skip past first of two arguments
		args = args[1:]

		// special case for argument value of 0
		if value < 1 {
			return zer
		}
		// limit value to between specified minimum and maximum
		if value < min {
			return min
		}
		if value > max {
			return max
		}
		return value
	}

	inSwitch := true

	// get concurrency, cleanup, and debugging flags in any order
	for {

		inSwitch = true

		switch args[0] {

		// concurrency override arguments can be passed in by local wrapper script (undocumented)
		case "-maxcpu":
			maxProcs = getNumericArg("Maximum number of processors", 1, 1, ncpu)
		case "-defcpu":
			defProcs = getNumericArg("Default number of processors", ncpu, 1, ncpu)
		// performance tuning flags
		case "-proc":
			numProcs = getNumericArg("Number of processors", ncpu, 1, ncpu)
		case "-cons":
			serverRatio = getNumericArg("Parser to processor ratio", 4, 1, 32)
		case "-serv":
			NumServe = getNumericArg("Concurrent parser count", 0, 1, 128)
		case "-chan":
			ChanDepth = getNumericArg("Communication channel depth", 0, ncpu, 128)
		case "-heap":
			HeapSize = getNumericArg("Unshuffler heap size", 8, 8, 64)
		case "-farm":
			FarmSize = getNumericArg("Node buffer length", 4, 4, 2048)
		case "-gogc":
			goGc = getNumericArg("Garbage collection percentage", 0, 50, 1000)
			gcdefault = false

		// read data from file
		case "-input":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Input file name is missing\n")
				os.Exit(1)
			}
			fileName = args[1]
			// skip past first of two arguments
			args = args[1:]

		// file with selected indexes for removing duplicates
		case "-unique":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Unique identifier file is missing\n")
				os.Exit(1)
			}
			unqe = args[1]
			// skip past first of two arguments
			args = args[1:]

		// local directory path for indexing
		case "-archive", "-stash":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Archive path is missing\n")
				os.Exit(1)
			}
			stsh = args[1]
			if stsh != "" && !strings.HasSuffix(stsh, "/") {
				stsh += "/"
			}
			// skip past first of two arguments
			args = args[1:]
		// local directory path for retrieval
		case "-fetch":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Fetch path is missing\n")
				os.Exit(1)
			}
			ftch = args[1]
			if ftch != "" && !strings.HasSuffix(ftch, "/") {
				ftch += "/"
			}
			// skip past first of two arguments
			args = args[1:]
		// local directory path for retrieval of compressed XML
		case "-stream":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Stream path is missing\n")
				os.Exit(1)
			}
			strm = args[1]
			if strm != "" && !strings.HasSuffix(strm, "/") {
				strm += "/"
			}
			// skip past first of two arguments
			args = args[1:]

		// local directory path for extra link retrieval
		case "-summon":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Summon path is missing\n")
				os.Exit(1)
			}
			smmn = args[1]
			if smmn != "" && !strings.HasSuffix(smmn, "/") {
				smmn += "/"
			}
			// skip past first of two arguments
			args = args[1:]

		// data element for indexing
		case "-index":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Index element is missing\n")
				os.Exit(1)
			}
			indx = args[1]
			// skip past first of two arguments
			args = args[1:]

		// build inverted index
		case "-invert":
			nvrt = true

		// combine sets of inverted index files
		case "-join":
			join = true

		case "-fuse":
			fuse = true

		// merge inverted index files, distribute by prefix
		case "-merge":
			if len(args) < 3 {
				fmt.Fprintf(os.Stderr, "\nERROR: Merge field is missing\n")
				os.Exit(1)
			}
			merg = args[1]
			// skip past first of two arguments
			args = args[1:]

		// promote inverted index to term-specific postings files
		case "-promote":
			if len(args) < 3 {
				fmt.Fprintf(os.Stderr, "\nERROR: Promote path is missing\n")
				os.Exit(1)
			}
			prom = args[1]
			fild = args[2]
			// skip past first and second arguments
			args = args[2:]

		case "-path":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Postings path is missing\n")
				os.Exit(1)
			}
			base = args[1]
			// skip past first of two arguments
			args = args[1:]

		case "-exact":
			xact = true
			fallthrough
		case "-search":
			rlxd = true
			fallthrough
		case "-query":
			if xact && rlxd {
				rlxd = false
			}
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Query argument is missing\n")
				os.Exit(1)
			}
			phrs = args[1]
			// skip past first of two arguments
			args = args[1:]

		case "-batch":
			btch = true

		case "-mockx":
			xact = true
			fallthrough
		case "-mocks":
			rlxd = true
			fallthrough
		case "-mock":
			if xact && rlxd {
				rlxd = false
			}
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Query argument is missing\n")
				os.Exit(1)
			}
			phrs = args[1]
			mock = true
			// skip past first of two arguments
			args = args[1:]

		// -countp tests the files containing positions of terms per UID (undocumented)
		case "-countp":
			psns = true
			fallthrough
		case "-counts":
			plrl = true
			fallthrough
		case "-countr":
			rlxd = true
			fallthrough
		case "-count":
			if plrl && rlxd {
				rlxd = false
			}
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Count argument is missing\n")
				os.Exit(1)
			}
			trms = args[1]
			// skip past first of two arguments
			args = args[1:]

		case "-gzip":
			zipp = true
		case "-hash":
			hshv = true
		case "-trie":
			trei = true
		// check for missing records
		case "-missing":
			msng = true

		// use non-threaded fetch function for windows (undocumented)
		case "-windows":
			windows = true

		// data cleanup flags
		case "-compress", "-compressed":
			DoCompress = true
		case "-spaces", "-cleanup":
			DoCleanup = true
		case "-strict":
			DoStrict = true
		case "-mixed":
			DoMixed = true
		case "-accent":
			DeAccent = true
		case "-ascii":
			DoASCII = true

		// previously visible processing flags (undocumented)
		case "-stems", "-stem":
			DoStem = true
		case "-stops", "-stop":
			DeStop = false

		case "-unicode", "-repair":
			DoUnicode = true
		case "-script":
			DoScript = true
		case "-mathml":
			DoMathML = true

		case "-flag", "-flags":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Flags argument is missing\n")
				os.Exit(1)
			}
			flgs = args[1]
			// skip past first of two arguments
			args = args[1:]

		// debugging flags
		case "-damaged", "-damage", "-broken":
			dmgd = true
			if len(args) > 1 {
				next := args[1]
				// if next argument is not another flag
				if next != "" && next[0] != '-' {
					// get optional extraction class (SELF, SINGLE, DOUBLE, AMPER, or ALL)
					dmgdType = next
					// skip past first of two arguments
					args = args[1:]
				}
			}
		case "-prepare":
			cmpr = true
			if len(args) > 1 {
				next := args[1]
				// if next argument is not another flag
				if next != "" && next[0] != '-' {
					// get optional data source specifier
					cmprType = next
					// skip past first of two arguments
					args = args[1:]
				}
			}
		case "-ignore":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: -ignore value is missing\n")
				os.Exit(1)
			}
			ignr = args[1]
			// skip past first of two arguments
			args = args[1:]

		// debugging flags
		case "-debug":
			dbug = true
		case "-stats", "-stat":
			stts = true
		case "-timer":
			timr = true
		case "-profile":
			prfl = true

		default:
			// if not any of the controls, set flag to break out of for loop
			inSwitch = false
		}

		if !inSwitch {
			break
		}

		// skip past argument
		args = args[1:]

		if len(args) < 1 {
			break
		}
	}

	// -flag allows script to set -strict or -mixed (or -stems, or -stops) from argument
	switch flgs {
	case "strict":
		DoStrict = true
	case "mixed":
		DoMixed = true
	case "stems", "stem":
		DoStem = true
	case "stops", "stop":
		DeStop = false
	case "none", "default":
	default:
		if flgs != "" {
			fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized -flag value '%s'\n", flgs)
			os.Exit(1)
		}
	}

	CountLines = DoMixed
	AllowEmbed = DoStrict || DoMixed
	ContentMods = AllowEmbed || DoCompress || DoUnicode || DoScript || DoMathML || DeAccent || DoASCII

	// reality checks on number of processors to use
	// performance degrades if capacity is above maximum number of partitions per second (context switching?)
	if numProcs == 0 {
		if defProcs > 0 {
			numProcs = defProcs
		} else {
			// best performance measurement with current code is obtained when 6 to 8 processors are assigned,
			// varying slightly among queries on PubmedArticle, gene DocumentSummary, and INSDSeq sequence records
			numProcs = 8
			if cpuid.CPU.ThreadsPerCore > 1 {
				cores := ncpu / cpuid.CPU.ThreadsPerCore
				if cores > 4 && cores < 8 {
					numProcs = cores
				}
			}
		}
	}
	if numProcs > ncpu {
		numProcs = ncpu
	}
	if numProcs > maxProcs {
		numProcs = maxProcs
	}

	// allow simultaneous threads for multiplexed goroutines
	runtime.GOMAXPROCS(numProcs)

	// adjust garbage collection target percentage
	if goGc >= 50 && goGc <= 1000 {
		debug.SetGCPercent(goGc)
	}

	// explicit -serv argument overrides -cons ratio
	if NumServe > 0 {
		serverRatio = NumServe / numProcs
		// if numServers / numProcs is not a whole number, do not print serverRatio in -stats
		if NumServe != numProcs*serverRatio {
			serverRatio = 0
		}
	} else {
		NumServe = numProcs * serverRatio
	}
	// server limits
	if NumServe > 128 {
		NumServe = 128
	} else if NumServe < 1 {
		NumServe = numProcs
	}

	// explicit -chan argument overrides default to number of servers
	if ChanDepth == 0 {
		ChanDepth = NumServe
	}

	// -stats prints number of CPUs and performance tuning values if no other arguments (undocumented)
	if stts && len(args) < 1 {

		fmt.Fprintf(os.Stderr, "Core %d\n", ncpu/cpuid.CPU.ThreadsPerCore)
		fmt.Fprintf(os.Stderr, "Thrd %d\n", ncpu)
		fmt.Fprintf(os.Stderr, "Sock %d\n", ncpu/cpuid.CPU.LogicalCores)
		fmt.Fprintf(os.Stderr, "Mmry %d\n", memory.TotalMemory()/(1024*1024*1024))

		fmt.Fprintf(os.Stderr, "Proc %d\n", numProcs)
		if serverRatio > 0 {
			fmt.Fprintf(os.Stderr, "Cons %d\n", serverRatio)
		}
		fmt.Fprintf(os.Stderr, "Serv %d\n", NumServe)
		fmt.Fprintf(os.Stderr, "Chan %d\n", ChanDepth)
		fmt.Fprintf(os.Stderr, "Heap %d\n", HeapSize)
		fmt.Fprintf(os.Stderr, "Farm %d\n", FarmSize)
		fmt.Fprintf(os.Stderr, "Gogc %d\n", goGc)

		fi, err := os.Stdin.Stat()
		if err == nil {
			mode := fi.Mode().String()
			fmt.Fprintf(os.Stderr, "Mode %s\n", mode)
		}

		fmt.Fprintf(os.Stderr, "\n")

		return
	}

	// if copying from local files accessed by identifier, add dummy argument to bypass length tests
	if stsh != "" && indx == "" {
		args = append(args, "-dummy")
	} else if ftch != "" || strm != "" || smmn != "" {
		args = append(args, "-dummy")
	} else if base != "" {
		args = append(args, "-dummy")
	} else if trei || dmgd || cmpr {
		args = append(args, "-dummy")
	}

	// expand -archive ~/ to home directory path
	if stsh != "" {

		if stsh[:2] == "~/" {
			cur, err := user.Current()
			if err == nil {
				hom := cur.HomeDir
				stsh = strings.Replace(stsh, "~/", hom+"/", 1)
			}
		}
	}

	// expand -fetch ~/ to home directory path
	if ftch != "" {

		if ftch[:2] == "~/" {
			cur, err := user.Current()
			if err == nil {
				hom := cur.HomeDir
				ftch = strings.Replace(ftch, "~/", hom+"/", 1)
			}
		}
	}

	// expand -stream ~/ to home directory path
	if strm != "" {

		if strm[:2] == "~/" {
			cur, err := user.Current()
			if err == nil {
				hom := cur.HomeDir
				strm = strings.Replace(strm, "~/", hom+"/", 1)
			}
		}
	}

	// expand -promote ~/ to home directory path
	if prom != "" {

		if prom[:2] == "~/" {
			cur, err := user.Current()
			if err == nil {
				hom := cur.HomeDir
				prom = strings.Replace(prom, "~/", hom+"/", 1)
			}
		}
	}

	// expand -summon ~/ to home directory path
	if smmn != "" {

		if smmn[:2] == "~/" {
			cur, err := user.Current()
			if err == nil {
				hom := cur.HomeDir
				smmn = strings.Replace(smmn, "~/", hom+"/", 1)
			}
		}
	}

	// DOCUMENTATION COMMANDS

	if len(args) > 0 {

		inSwitch = true

		switch args[0] {
		case "-version":
			fmt.Printf("%s\n", rchiveVersion)
		case "-help":
			fmt.Printf("rchive %s\n%s\n", rchiveVersion, rchiveHelp)
		case "-extras", "-extra", "-advanced":
			fmt.Printf("rchive %s\n%s\n", rchiveVersion, rchiveExtras)
		case "-internal", "-internals":
			fmt.Printf("rchive %s\n%s\n", rchiveVersion, rchiveInternal)
		default:
			// if not any of the documentation commands, keep going
			inSwitch = false
		}

		if inSwitch {
			return
		}
	}

	// FILE NAME CAN BE SUPPLIED WITH -input COMMAND

	in := os.Stdin

	// check for data being piped into stdin
	isPipe := false
	fi, staterr := os.Stdin.Stat()
	if staterr == nil {
		isPipe = bool((fi.Mode() & os.ModeNamedPipe) != 0)
	}

	usingFile := false

	if fileName != "" {

		inFile, err := os.Open(fileName)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open input file '%s'\n", fileName)
			os.Exit(1)
		}

		defer inFile.Close()

		// use indicated file instead of stdin
		in = inFile
		usingFile = true

		if isPipe && runtime.GOOS != "windows" {
			mode := fi.Mode().String()
			fmt.Fprintf(os.Stderr, "\nERROR: Input data from both stdin and file '%s', mode is '%s'\n", fileName, mode)
			os.Exit(1)
		}
	}

	// check for -input command after extraction arguments
	for _, str := range args {
		if str == "-input" {
			fmt.Fprintf(os.Stderr, "\nERROR: Misplaced -input command\n")
			os.Exit(1)
		}
	}

	// START PROFILING IF REQUESTED

	if prfl {

		f, err := os.Create("cpu.pprof")
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create profile output file\n")
			os.Exit(1)
		}

		pprof.StartCPUProfile(f)

		defer pprof.StopCPUProfile()
	}

	// INITIALIZE PROCESS TIMER AND RECORD COUNT

	startTime := time.Now()
	recordCount := 0
	byteCount := 0

	// print processing rate and program duration
	printDuration := func(name string) {

		stopTime := time.Now()
		duration := stopTime.Sub(startTime)
		seconds := float64(duration.Nanoseconds()) / 1e9

		if recordCount >= 1000000 {
			throughput := float64(recordCount/100000) / 10.0
			fmt.Fprintf(os.Stderr, "\nRchive processed %.1f million %s in %.3f seconds", throughput, name, seconds)
		} else {
			fmt.Fprintf(os.Stderr, "\nRchive processed %d %s in %.3f seconds", recordCount, name, seconds)
		}

		if seconds >= 0.001 && recordCount > 0 {
			rate := int(float64(recordCount) / seconds)
			if rate >= 1000000 {
				fmt.Fprintf(os.Stderr, " (%d million %s/second", rate/1000000, name)
			} else {
				fmt.Fprintf(os.Stderr, " (%d %s/second", rate, name)
			}
			if byteCount > 0 {
				rate := int(float64(byteCount) / seconds)
				if rate >= 1000000 {
					fmt.Fprintf(os.Stderr, ", %d megabytes/second", rate/1000000)
				} else if rate >= 1000 {
					fmt.Fprintf(os.Stderr, ", %d kilobytes/second", rate/1000)
				} else {
					fmt.Fprintf(os.Stderr, ", %d bytes/second", rate)
				}
			}
			fmt.Fprintf(os.Stderr, ")")
		}

		fmt.Fprintf(os.Stderr, "\n\n")
	}

	// EXTERNAL INDEXERS AND LINK ARCHIVER

	if len(args) > 0 {
		switch args[0] {
		case "-bioconcepts", "-generif", "-generifs":
			recordCount = CreateExternalIndexer(args, zipp, in)

			debug.FreeOSMemory()

			if timr {
				printDuration("records")
			}

			return
		case "-theme", "-themes", "-dpath", "-dpaths", "-thesis":
			recordCount = CreateExternalIndexer(args, zipp, in)

			debug.FreeOSMemory()

			if timr {
				printDuration("lines")
			}

			return
		default:
		}
	}

	if len(args) > 1 {
		switch args[0] {
		case "-distribute":
			args = args[1:]

			// first argument is archive path
			path := args[0]
			if path == "" {
				fmt.Fprintf(os.Stderr, "\nERROR: Need path in order to create extra index stasher\n")
				os.Exit(1)
			}
			if path[:2] == "~/" {
				cur, err := user.Current()
				if err == nil {
					hom := cur.HomeDir
					path = strings.Replace(path, "~/", hom+"/", 1)
				}
			}

			// remaining arguments are *.e2x files
			// e.g., rchive -timer -distribute archive_directory *.e2x
			args = args[1:]
			stsq := CreateExternalArchive(path, args)

			if stsq == nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create extra index stasher\n")
				os.Exit(1)
			}

			// drain output channel
			for str := range stsq {

				if hshv {
					// print table of UIDs and hash values
					os.Stdout.WriteString(str)
				}

				recordCount++
				runtime.Gosched()
			}

			debug.FreeOSMemory()

			if timr {
				printDuration("records")
			}

			return
		default:
		}
	}

	// JOIN SUBSETS OF INVERTED INDEX FILES

	// -join combines subsets of inverted files for subsequent -merge operation
	if join {

		// environment variable can override garbage collector (undocumented)
		gcEnv := os.Getenv("EDIRECT_JOIN_GOGC")
		if gcEnv != "" {
			val, err := strconv.Atoi(gcEnv)
			if err == nil {
				if val >= 50 && val <= 1000 {
					debug.SetGCPercent(val)
				} else {
					debug.SetGCPercent(100)
				}
			}
		} else if gcdefault {
			// default to 100 for join and merge
			debug.SetGCPercent(100)
		}

		// environment variable can override number of servers (undocumented)
		svEnv := os.Getenv("EDIRECT_JOIN_SERV")
		if svEnv != "" {
			val, err := strconv.Atoi(svEnv)
			if err == nil {
				if val >= 1 && val <= 128 {
					NumServe = val
				} else {
					NumServe = 1
				}
			}
		}

		chns := CreatePresenters(args)
		mfld := CreateManifold(chns)
		mrgr := CreateMergers(mfld)
		unsq := CreateUnshuffler(mrgr)

		if chns == nil || mfld == nil || mrgr == nil || unsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create inverted index joiner\n")
			os.Exit(1)
		}

		if dbug {

			// drain results, but suppress normal output
			for range unsq {
				recordCount++
				runtime.Gosched()
			}

			// force garbage collection, return memory to operating system
			debug.FreeOSMemory()

			// print processing parameters as XML object
			stopTime := time.Now()
			duration := stopTime.Sub(startTime)
			seconds := float64(duration.Nanoseconds()) / 1e9

			// Threads is a more easily explained concept than GOMAXPROCS
			fmt.Printf("<Xtract>\n")
			fmt.Printf("  <Threads>%d</Threads>\n", numProcs)
			fmt.Printf("  <Parsers>%d</Parsers>\n", NumServe)
			fmt.Printf("  <Time>%.3f</Time>\n", seconds)
			if seconds >= 0.001 && recordCount > 0 {
				rate := int(float64(recordCount) / seconds)
				fmt.Printf("  <Rate>%d</Rate>\n", rate)
			}
			fmt.Printf("</Xtract>\n")

			return
		}

		var out io.Writer

		out = os.Stdout

		if zipp {

			zpr, err := gzip.NewWriterLevel(out, gzip.BestSpeed)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create compressor\n")
				os.Exit(1)
			}

			// close decompressor when all records have been processed
			defer zpr.Close()

			// use compressor for writing file
			out = zpr
		}

		// create buffered writer layer
		wrtr := bufio.NewWriter(out)

		wrtr.WriteString("<InvDocumentSet>\n")

		// drain channel of alphabetized results
		for curr := range unsq {

			str := curr.Text

			if str == "" {
				continue
			}

			// send result to output
			wrtr.WriteString(str)

			recordCount++
			runtime.Gosched()
		}

		wrtr.WriteString("</InvDocumentSet>\n\n")

		wrtr.Flush()
	}

	// MERGE INVERTED INDEX FILES AND GROUP BY TERM

	// -merge combines inverted files, distributes by prefix
	if merg != "" {

		// environment variable can override garbage collector (undocumented)
		gcEnv := os.Getenv("EDIRECT_MERGE_GOGC")
		if gcEnv != "" {
			val, err := strconv.Atoi(gcEnv)
			if err == nil {
				if val >= 50 && val <= 1000 {
					debug.SetGCPercent(val)
				} else {
					debug.SetGCPercent(100)
				}
			}
		} else if gcdefault {
			// default to 100 for join and merge
			debug.SetGCPercent(100)
		}

		// environment variable can override number of servers (undocumented)
		svEnv := os.Getenv("EDIRECT_MERGE_SERV")
		if svEnv != "" {
			val, err := strconv.Atoi(svEnv)
			if err == nil {
				if val >= 1 && val <= 128 {
					NumServe = val
				} else {
					NumServe = 1
				}
			}
		}

		chns := CreatePresenters(args)
		mfld := CreateManifold(chns)
		mrgr := CreateMergers(mfld)
		unsq := CreateUnshuffler(mrgr)
		sptr := CreateSplitter(merg, zipp, unsq)

		if chns == nil || mfld == nil || mrgr == nil || unsq == nil || sptr == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create inverted index merger\n")
			os.Exit(1)
		}

		if dbug {

			// drain results, but suppress normal output
			for range sptr {
				recordCount++
				runtime.Gosched()
			}

			// force garbage collection, return memory to operating system
			debug.FreeOSMemory()

			// print processing parameters as XML object
			stopTime := time.Now()
			duration := stopTime.Sub(startTime)
			seconds := float64(duration.Nanoseconds()) / 1e9

			// Threads is a more easily explained concept than GOMAXPROCS
			fmt.Printf("<Xtract>\n")
			fmt.Printf("  <Threads>%d</Threads>\n", numProcs)
			fmt.Printf("  <Parsers>%d</Parsers>\n", NumServe)
			fmt.Printf("  <Time>%.3f</Time>\n", seconds)
			if seconds >= 0.001 && recordCount > 0 {
				rate := int(float64(recordCount) / seconds)
				fmt.Printf("  <Rate>%d</Rate>\n", rate)
			}
			fmt.Printf("</Xtract>\n")

			return
		}

		// drain channel, print two-character index name
		for str := range sptr {

			fmt.Fprintf(os.Stdout, "%s\n", str)

			recordCount++
			runtime.Gosched()
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("terms")
		}

		return
	}

	// PROMOTE MERGED INVERTED INDEX TO TERM LIST AND POSTINGS FILES

	if prom != "" && fild != "" {

		prmq := CreatePromoters(args, prom, fild)

		if prmq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create new postings file generator\n")
			os.Exit(1)
		}

		// drain channel, print 2-4 character file prefix
		for str := range prmq {

			fmt.Fprintf(os.Stdout, "%s\n", str)

			recordCount++
			runtime.Gosched()
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("terms")
		}

		return
	}

	// QUERY POSTINGS FILES

	if phrs != "" || trms != "" || btch {
		if base == "" {
			// obtain path from environment variable within rchive as a convenience
			base = os.Getenv("EDIRECT_PUBMED_MASTER")
			if base != "" {
				if !strings.HasSuffix(base, "/") {
					base += "/"
				}
				base += "Postings"
			}
		}
	}

	if base != "" && btch {

		// read query lines for exact match
		scanr := bufio.NewScanner(in)

		for scanr.Scan() {
			txt := scanr.Text()

			// deStop should match value used in building the indices
			recordCount += ProcessSearch(base, txt, true, false)
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	if base != "" && phrs != "" {

		// deStop should match value used in building the indices
		if mock {
			recordCount = ProcessMock(base, phrs, xact, rlxd)
		} else {
			recordCount = ProcessSearch(base, phrs, xact, rlxd)
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	if base != "" && trms != "" {

		// deStop should match value used in building the indices
		recordCount = ProcessCount(base, trms, plrl, psns, rlxd)

		debug.FreeOSMemory()

		if timr {
			printDuration("terms")
		}

		return
	}

	// CONFIRM INPUT DATA AVAILABILITY AFTER RUNNING COMMAND GENERATORS

	if fileName == "" && runtime.GOOS != "windows" {

		fromStdin := bool((fi.Mode() & os.ModeCharDevice) == 0)
		if !isPipe || !fromStdin {
			mode := fi.Mode().String()
			fmt.Fprintf(os.Stderr, "\nERROR: No data supplied to rchive from stdin or file, mode is '%s'\n", mode)
			os.Exit(1)
		}
	}

	if !usingFile && !isPipe {

		fmt.Fprintf(os.Stderr, "\nERROR: No XML input data supplied to rchive\n")
		os.Exit(1)
	}

	// SPECIFY STRINGS TO GO BEFORE AND AFTER ENTIRE OUTPUT OR EACH RECORD

	head := ""
	tail := ""

	hd := ""
	tl := ""

	for {

		if len(args) < 1 {
			break
		}

		inSwitch = true

		switch args[0] {
		case "-head":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -head command\n")
				os.Exit(1)
			}
			head = ConvertSlash(args[1])
		case "-tail":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -tail command\n")
				os.Exit(1)
			}
			tail = ConvertSlash(args[1])
		case "-hd":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -hd command\n")
				os.Exit(1)
			}
			hd = ConvertSlash(args[1])
		case "-tl":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -tl command\n")
				os.Exit(1)
			}
			tl = ConvertSlash(args[1])
		default:
			// if not any of the controls, set flag to break out of for loop
			inSwitch = false
		}

		if !inSwitch {
			break
		}

		// skip past arguments
		args = args[2:]
	}

	// PRODUCE ARCHIVE SUBPATH FROM IDENTIFIER

	// -trie converts identifier to directory subpath plus file name (undocumented)
	if trei {

		scanr := bufio.NewScanner(in)

		sfx := ".xml"
		if zipp {
			sfx += ".gz"
		}

		// read lines of identifiers
		for scanr.Scan() {

			file := scanr.Text()

			var arry [132]rune
			trie := MakeArchiveTrie(file, arry)
			if trie == "" || file == "" {
				continue
			}

			fpath := path.Join(trie, file+sfx)
			if fpath == "" {
				continue
			}

			os.Stdout.WriteString(fpath)
			os.Stdout.WriteString("\n")
		}

		return
	}

	// CHECK FOR MISSING RECORDS IN LOCAL DIRECTORY INDEXED BY TRIE ON IDENTIFIER

	// -archive plus -missing checks for missing records
	if stsh != "" && msng {

		scanr := bufio.NewScanner(in)

		sfx := ".xml"
		if zipp {
			sfx += ".gz"
		}

		// read lines of identifiers
		for scanr.Scan() {

			file := scanr.Text()

			pos := strings.Index(file, ".")
			if pos >= 0 {
				// remove version suffix
				file = file[:pos]
			}

			var arry [132]rune
			trie := MakeArchiveTrie(file, arry)

			if file == "" || trie == "" {
				continue
			}

			fpath := path.Join(stsh, trie, file+sfx)
			if fpath == "" {
				continue
			}

			_, err := os.Stat(fpath)

			// if failed to find ".xml" file, try ".xml.gz" without requiring -gzip
			if err != nil && os.IsNotExist(err) && !zipp {
				fpath := path.Join(stsh, trie, file+".xml.gz")
				if fpath == "" {
					continue
				}
				_, err = os.Stat(fpath)
			}
			if err != nil && os.IsNotExist(err) {
				// record is missing from local file cache
				os.Stdout.WriteString(file)
				os.Stdout.WriteString("\n")
			}
		}

		return
	}

	// RETRIEVE XML COMPONENT RECORDS FROM LOCAL DIRECTORY INDEXED BY TRIE ON IDENTIFIER

	// alternative windows version limits memory by not using goroutines
	if ftch != "" && indx == "" && runtime.GOOS == "windows" && windows {

		scanr := bufio.NewScanner(in)
		if scanr == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create UID scanner\n")
			os.Exit(1)
		}

		sfx := ".xml"
		if zipp {
			sfx += ".gz"
		}

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		var buf bytes.Buffer

		for scanr.Scan() {

			// read next identifier
			file := scanr.Text()

			pos := strings.Index(file, ".")
			if pos >= 0 {
				// remove version suffix
				file = file[:pos]
			}

			var arry [132]rune
			trie := MakeArchiveTrie(file, arry)

			if file == "" || trie == "" {
				continue
			}

			fpath := path.Join(ftch, trie, file+sfx)
			if fpath == "" {
				continue
			}

			iszip := zipp

			inFile, err := os.Open(fpath)

			// if failed to find ".xml" file, try ".xml.gz" without requiring -gzip
			if err != nil && os.IsNotExist(err) && !zipp {
				iszip = true
				fpath := path.Join(ftch, trie, file+".xml.gz")
				if fpath == "" {
					continue
				}
				inFile, err = os.Open(fpath)
			}
			if err != nil {
				continue
			}

			buf.Reset()

			brd := bufio.NewReader(inFile)

			if iszip {

				zpr, err := gzip.NewReader(brd)

				if err == nil {
					// copy and decompress cached file contents
					buf.ReadFrom(zpr)
				}

				zpr.Close()

			} else {

				// copy cached file contents
				buf.ReadFrom(brd)
			}

			inFile.Close()

			str := buf.String()

			if str == "" {
				continue
			}

			recordCount++

			if hd != "" {
				os.Stdout.WriteString(hd)
				os.Stdout.WriteString("\n")
			}

			if hshv {
				// calculate hash code for verification table
				hsh := crc32.NewIEEE()
				hsh.Write([]byte(str))
				val := hsh.Sum32()
				res := strconv.FormatUint(uint64(val), 10)
				txt := file + "\t" + res + "\n"
				os.Stdout.WriteString(txt)
			} else {
				// send result to output
				os.Stdout.WriteString(str)
				if !strings.HasSuffix(str, "\n") {
					os.Stdout.WriteString("\n")
				}
			}

			if tl != "" {
				os.Stdout.WriteString(tl)
				os.Stdout.WriteString("\n")
			}

			debug.FreeOSMemory()
		}

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// -fetch without -index retrieves XML files in trie-based directory structure
	if ftch != "" && indx == "" {

		uidq := CreateUIDReader(in)
		strq := CreateFetchers(ftch, ".xml", zipp, uidq)
		unsq := CreateUnshuffler(strq)

		if uidq == nil || strq == nil || unsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create archive reader\n")
			os.Exit(1)
		}

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		// drain output channel
		for curr := range unsq {

			str := curr.Text

			if str == "" {
				continue
			}

			if hd != "" {
				os.Stdout.WriteString(hd)
				os.Stdout.WriteString("\n")
			}

			if hshv {
				// calculate hash code for verification table
				hsh := crc32.NewIEEE()
				hsh.Write([]byte(str))
				val := hsh.Sum32()
				res := strconv.FormatUint(uint64(val), 10)
				txt := curr.Ident + "\t" + res + "\n"
				os.Stdout.WriteString(txt)
			} else {
				// send result to output
				os.Stdout.WriteString(str)
				if !strings.HasSuffix(str, "\n") {
					os.Stdout.WriteString("\n")
				}
			}

			if tl != "" {
				os.Stdout.WriteString(tl)
				os.Stdout.WriteString("\n")
			}

			recordCount++
			runtime.Gosched()
		}

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// -stream without -index retrieves compressed XML files in trie-based directory structure
	if strm != "" && indx == "" {

		uidq := CreateUIDReader(in)
		strq := CreateStreamers(strm, uidq)
		unsq := CreateUnshuffler(strq)

		if uidq == nil || strq == nil || unsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create archive reader\n")
			os.Exit(1)
		}

		// drain output channel
		for curr := range unsq {

			data := curr.Data

			if data == nil {
				continue
			}

			recordCount++
			runtime.Gosched()

			_, err := os.Stdout.Write(data)
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			}
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// -summon retrieves link files in trie-based directory structure
	if smmn != "" && indx == "" {

		uidq := CreateUIDReader(in)
		strq := CreateFetchers(smmn, ".e2x", zipp, uidq)
		unsq := CreateUnshuffler(strq)

		if uidq == nil || strq == nil || unsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create link reader\n")
			os.Exit(1)
		}

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		// drain output channel
		for curr := range unsq {

			str := curr.Text

			if str == "" {
				continue
			}

			if hd != "" {
				os.Stdout.WriteString(hd)
				os.Stdout.WriteString("\n")
			}

			if hshv {
				// calculate hash code for verification table
				hsh := crc32.NewIEEE()
				hsh.Write([]byte(str))
				val := hsh.Sum32()
				res := strconv.FormatUint(uint64(val), 10)
				txt := curr.Ident + "\t" + res + "\n"
				os.Stdout.WriteString(txt)
			} else {
				// send result to output
				os.Stdout.WriteString(str)
				if !strings.HasSuffix(str, "\n") {
					os.Stdout.WriteString("\n")
				}
			}

			if tl != "" {
				os.Stdout.WriteString(tl)
				os.Stdout.WriteString("\n")
			}

			recordCount++
			runtime.Gosched()
		}

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// CREATE XML BLOCK READER FROM STDIN OR FILE

	rdr := CreateReader(in)
	if rdr == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create XML Block Reader\n")
		os.Exit(1)
	}

	// ENTREZ INDEX INVERSION

	// -invert reads IdxDocumentSet XML and creates an inverted index
	if nvrt {

		// environment variable can override garbage collector (undocumented)
		gcEnv := os.Getenv("EDIRECT_INVERT_GOGC")
		if gcEnv != "" {
			val, err := strconv.Atoi(gcEnv)
			if err == nil {
				if val >= 50 && val <= 1000 {
					debug.SetGCPercent(val)
				} else {
					debug.SetGCPercent(100)
				}
			}
		}

		// environment variable can override number of servers (undocumented)
		svEnv := os.Getenv("EDIRECT_INVERT_SERV")
		if svEnv != "" {
			val, err := strconv.Atoi(svEnv)
			if err == nil {
				if val >= 1 && val <= 128 {
					NumServe = val
				} else {
					NumServe = 1
				}
			}
		}

		colq := CreateProducer("IdxDocument", "", rdr)
		dspq := CreateDispensers(colq)
		invq := CreateInverters(dspq)
		rslq := CreateResolver(invq)

		if colq == nil || dspq == nil || invq == nil || rslq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create inverter\n")
			os.Exit(1)
		}

		if dbug {

			// drain results, but suppress normal output
			for range rslq {
				recordCount++
				runtime.Gosched()
			}

			// force garbage collection, return memory to operating system
			debug.FreeOSMemory()

			// print processing parameters as XML object
			stopTime := time.Now()
			duration := stopTime.Sub(startTime)
			seconds := float64(duration.Nanoseconds()) / 1e9

			// Threads is a more easily explained concept than GOMAXPROCS
			fmt.Printf("<Xtract>\n")
			fmt.Printf("  <Threads>%d</Threads>\n", numProcs)
			fmt.Printf("  <Parsers>%d</Parsers>\n", NumServe)
			fmt.Printf("  <Time>%.3f</Time>\n", seconds)
			if seconds >= 0.001 && recordCount > 0 {
				rate := int(float64(recordCount) / seconds)
				fmt.Printf("  <Rate>%d</Rate>\n", rate)
			}
			fmt.Printf("</Xtract>\n")

			return
		}

		var out io.Writer

		out = os.Stdout

		if zipp {

			zpr, err := gzip.NewWriterLevel(out, gzip.BestSpeed)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create compressor\n")
				os.Exit(1)
			}

			// close decompressor when all records have been processed
			defer zpr.Close()

			// use compressor for writing file
			out = zpr
		}

		// create buffered writer layer
		wrtr := bufio.NewWriter(out)

		wrtr.WriteString("<InvDocumentSet>\n")

		// drain channel of alphabetized results
		for str := range rslq {

			// send result to output
			wrtr.WriteString(str)

			recordCount++
			runtime.Gosched()
		}

		wrtr.WriteString("</InvDocumentSet>\n\n")

		wrtr.Flush()

		debug.FreeOSMemory()

		if timr {
			printDuration("terms")
		}

		return
	}

	// FUSE SUBSETS OF INVERTED INDEX FILES

	// -fuse combines subsets of inverted files for subsequent -merge operation
	if fuse {

		// environment variable can override garbage collector (undocumented)
		gcEnv := os.Getenv("EDIRECT_FUSE_GOGC")
		if gcEnv != "" {
			val, err := strconv.Atoi(gcEnv)
			if err == nil {
				if val >= 50 && val <= 1000 {
					debug.SetGCPercent(val)
				} else {
					debug.SetGCPercent(100)
				}
			}
		} else if gcdefault {
			// default to 100 for fuse and merge
			debug.SetGCPercent(100)
		}

		// environment variable can override number of servers (undocumented)
		svEnv := os.Getenv("EDIRECT_FUSE_SERV")
		if svEnv != "" {
			val, err := strconv.Atoi(svEnv)
			if err == nil {
				if val >= 1 && val <= 128 {
					NumServe = val
				} else {
					NumServe = 1
				}
			}
		}

		chns := CreateProducer("InvDocument", "", rdr)
		fusr := CreateFusers(chns)
		mrgr := CreateMergers(fusr)
		unsq := CreateUnshuffler(mrgr)

		if chns == nil || fusr == nil || mrgr == nil || unsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create inverted index fuser\n")
			os.Exit(1)
		}

		var out io.Writer

		out = os.Stdout

		if zipp {

			zpr, err := gzip.NewWriterLevel(out, gzip.BestSpeed)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create compressor\n")
				os.Exit(1)
			}

			// close decompressor when all records have been processed
			defer zpr.Close()

			// use compressor for writing file
			out = zpr
		}

		// create buffered writer layer
		wrtr := bufio.NewWriter(out)

		wrtr.WriteString("<InvDocumentSet>\n")

		// drain channel of alphabetized results
		for curr := range unsq {

			str := curr.Text

			if str == "" {
				continue
			}

			// send result to output
			wrtr.WriteString(str)

			recordCount++
			runtime.Gosched()
		}

		wrtr.WriteString("</InvDocumentSet>\n\n")

		wrtr.Flush()

		debug.FreeOSMemory()

		if timr {
			printDuration("terms")
		}

		return
	}

	// ENSURE PRESENCE OF PATTERN ARGUMENT

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to rchive\n")
		os.Exit(1)
	}

	// allow -record as synonym of -pattern (undocumented)
	if args[0] == "-record" || args[0] == "-Record" {
		args[0] = "-pattern"
	}

	// make sure top-level -pattern command is next
	if args[0] != "-pattern" && args[0] != "-Pattern" {
		fmt.Fprintf(os.Stderr, "\nERROR: No -pattern in command-line arguments\n")
		os.Exit(1)
	}
	if len(args) < 2 {
		fmt.Fprintf(os.Stderr, "\nERROR: Item missing after -pattern command\n")
		os.Exit(1)
	}

	topPat := args[1]
	if topPat == "" {
		fmt.Fprintf(os.Stderr, "\nERROR: Item missing after -pattern command\n")
		os.Exit(1)
	}
	if strings.HasPrefix(topPat, "-") {
		fmt.Fprintf(os.Stderr, "\nERROR: Misplaced %s command\n", topPat)
		os.Exit(1)
	}

	// look for -pattern Parent/* construct for heterogeneous data, e.g., -pattern PubmedArticleSet/*
	topPattern, star := SplitInTwoAt(topPat, "/", LEFT)
	if topPattern == "" {
		return
	}

	parent := ""
	if star == "*" {
		parent = topPattern
	} else if star != "" {
		fmt.Fprintf(os.Stderr, "\nERROR: -pattern Parent/Child construct is not supported\n")
		os.Exit(1)
	}

	// FILTER XML RECORDS BY PRESENCE OF ONE OR MORE PHRASES

	// -pattern plus -phrase (-require, -exclude) filters by phrase in XML (undocumented)
	if len(args) > 2 && (args[2] == "-phrase" || args[2] == "-require" || args[2] == "-exclude") {

		exclude := false
		if args[2] == "-exclude" {
			exclude = true
		}

		if len(args) < 4 {
			fmt.Fprintf(os.Stderr, "\nERROR: Missing argument after %s\n", args[2])
			os.Exit(1)
		} else if len(args) > 4 {
			fmt.Fprintf(os.Stderr, "\nERROR: No arguments allowed after %s value\n", args[2])
			os.Exit(1)
		}

		// phrase to find anywhere in XML
		phrs := args[3]

		// convert old "+" phrase separator to new "AND" convention for entrez-phrase-search backward compatibility
		phrs = strings.Replace(phrs, " + ", " AND ", -1)
		// remove wildcard asterisk characters
		phrs = strings.Replace(phrs, "*", " ", -1)
		phrs = CompressRunsOfSpaces(phrs)
		phrs = strings.TrimSpace(phrs)

		if phrs == "" {
			fmt.Fprintf(os.Stderr, "\nERROR: Missing argument after %s\n", args[2])
			os.Exit(1)
		}

		xmlq := CreateProducer(topPattern, star, rdr)
		mchq := CreateMatchers(phrs, exclude, xmlq)
		unsq := CreateUnshuffler(mchq)

		if xmlq == nil || mchq == nil || unsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create phrase matcher\n")
			os.Exit(1)
		}

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		// drain output channel
		for curr := range unsq {

			str := curr.Text

			if str == "" {
				continue
			}

			if hd != "" {
				os.Stdout.WriteString(hd)
				os.Stdout.WriteString("\n")
			}

			// send result to output
			os.Stdout.WriteString(str)
			if !strings.HasSuffix(str, "\n") {
				os.Stdout.WriteString("\n")
			}

			if tl != "" {
				os.Stdout.WriteString(tl)
				os.Stdout.WriteString("\n")
			}

			recordCount++
			runtime.Gosched()
		}

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// REPORT RECORDS THAT CONTAIN DAMAGED EMBEDDED HTML TAGS

	// -damaged plus -index plus -pattern reports records with multiply-encoded HTML tags
	if dmgd && indx != "" {

		find := ParseIndex(indx)

		PartitionPattern(topPattern, star, rdr,
			func(str string) {
				recordCount++

				id := FindIdentifier(str[:], parent, find)
				if id == "" {
					return
				}

				// remove default version suffix
				if strings.HasSuffix(id, ".1") {
					idlen := len(id)
					id = id[:idlen-2]
				}

				ReportEncodedMarkup(dmgdType, id, str)
			})

		if timr {
			printDuration("records")
		}

		return
	}

	// COMPARE XML UPDATES TO LOCAL DIRECTORY, RETAIN NEW OR SUBSTANTIVELY CHANGED RECORDS

	// -prepare plus -archive plus -index plus -pattern compares XML files against stash
	if stsh != "" && indx != "" && cmpr {

		doReport := false
		if cmprType == "" || cmprType == "report" {
			doReport = true
		} else if cmprType != "release" {
			fmt.Fprintf(os.Stderr, "\nERROR: -prepare argument must be release or report\n")
			os.Exit(1)
		}

		find := ParseIndex(indx)

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		PartitionPattern(topPattern, star, rdr,
			func(str string) {
				recordCount++

				id := FindIdentifier(str[:], parent, find)
				if id == "" {
					return
				}

				pos := strings.Index(id, ".")
				if pos >= 0 {
					// remove version suffix
					id = id[:pos]
				}

				var arry [132]rune
				trie := MakeArchiveTrie(id, arry)

				if id == "" || trie == "" {
					return
				}

				fpath := path.Join(stsh, trie, id+".xml")
				if fpath == "" {
					return
				}

				// print new or updated XML record
				printRecord := func(stn string, isNew bool) {

					if stn == "" {
						return
					}

					if doReport {
						if isNew {
							os.Stdout.WriteString("NW ")
							os.Stdout.WriteString(id)
							os.Stdout.WriteString("\n")
						} else {
							os.Stdout.WriteString("UP ")
							os.Stdout.WriteString(id)
							os.Stdout.WriteString("\n")
						}
						return
					}

					if hd != "" {
						os.Stdout.WriteString(hd)
						os.Stdout.WriteString("\n")
					}

					os.Stdout.WriteString(stn)
					os.Stdout.WriteString("\n")

					if tl != "" {
						os.Stdout.WriteString(tl)
						os.Stdout.WriteString("\n")
					}
				}

				_, err := os.Stat(fpath)
				if err != nil && os.IsNotExist(err) {
					// new record
					printRecord(str, true)
					return
				}
				if err != nil {
					return
				}

				buf, err := ioutil.ReadFile(fpath)
				if err != nil {
					return
				}

				txt := string(buf[:])
				txt = strings.TrimSuffix(txt, "\n")

				// check for optional -ignore argument
				if ignr != "" {

					// ignore differences inside specified object
					ltag := "<" + ignr + ">"
					sleft, _ := SplitInTwoAt(str, ltag, LEFT)
					tleft, _ := SplitInTwoAt(txt, ltag, LEFT)

					rtag := "</" + ignr + ">"
					_, srght := SplitInTwoAt(str, rtag, RIGHT)
					_, trght := SplitInTwoAt(txt, rtag, RIGHT)

					if sleft == tleft && srght == trght {
						if doReport {
							os.Stdout.WriteString("NO ")
							os.Stdout.WriteString(id)
							os.Stdout.WriteString("\n")
						}
						return
					}

				} else {

					// compare entirety of objects
					if str == txt {
						if doReport {
							os.Stdout.WriteString("NO ")
							os.Stdout.WriteString(id)
							os.Stdout.WriteString("\n")
						}
						return
					}
				}

				// substantively modified record
				printRecord(str, false)
			})

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		if timr {
			printDuration("records")
		}

		return
	}

	// SAVE XML COMPONENT RECORDS TO LOCAL DIRECTORY INDEXED BY TRIE ON IDENTIFIER

	// -archive plus -index plus -pattern saves XML files in trie-based directory structure
	if stsh != "" && indx != "" {

		xmlq := CreateProducer(topPattern, star, rdr)
		stsq := CreateStashers(stsh, parent, indx, ".xml", hshv, zipp, 1000, xmlq)

		if xmlq == nil || stsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create stash generator\n")
			os.Exit(1)
		}

		// drain output channel
		for str := range stsq {

			if hshv {
				// print table of UIDs and hash values
				os.Stdout.WriteString(str)
			}

			recordCount++
			runtime.Gosched()
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// READ FILE OF IDENTIFIERS AND EXTRACT SELECTED RECORDS FROM XML INPUT FILE

	// -index plus -unique [plus -head/-tail/-hd/-tl] plus -pattern with no other extraction arguments
	// takes an XML input file and a file of its UIDs and keeps only the last version of each record
	if indx != "" && unqe != "" && len(args) == 2 {

		// read file of identifiers to use for filtering
		fl, err := os.Open(unqe)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open identifier file '%s'\n", unqe)
			os.Exit(1)
		}

		// create map that counts instances of each UID
		order := make(map[string]int)

		scanr := bufio.NewScanner(fl)

		// read lines of identifiers
		for scanr.Scan() {

			id := scanr.Text()

			// map records count for given identifier
			val := order[id]
			val++
			order[id] = val
		}

		fl.Close()

		find := ParseIndex(indx)

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		PartitionPattern(topPattern, star, rdr,
			func(str string) {
				recordCount++

				id := FindIdentifier(str[:], parent, find)
				if id == "" {
					return
				}

				val, ok := order[id]
				if !ok {
					// not in identifier list, skip
					return
				}
				// decrement count in map
				val--
				order[id] = val
				if val > 0 {
					// only write last record with a given identifier
					return
				}

				if hd != "" {
					os.Stdout.WriteString(hd)
					os.Stdout.WriteString("\n")
				}

				// write selected record
				os.Stdout.WriteString(str[:])
				os.Stdout.WriteString("\n")

				if tl != "" {
					os.Stdout.WriteString(tl)
					os.Stdout.WriteString("\n")
				}
			})

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		if timr {
			printDuration("records")
		}

		return
	}

	// GENERATE RECORD INDEX ON XML INPUT FILE

	// -index plus -pattern prints record identifier and XML size
	if indx != "" {

		lbl := ""
		// check for optional filename label after -pattern argument (undocumented)
		if len(args) > 3 && args[2] == "-lbl" {
			lbl = args[3]

			lbl = strings.TrimSpace(lbl)
			if strings.HasPrefix(lbl, "pubmed") {
				lbl = lbl[7:]
			}
			if strings.HasSuffix(lbl, ".xml.gz") {
				xlen := len(lbl)
				lbl = lbl[:xlen-7]
			}
			lbl = strings.TrimSpace(lbl)
		}

		// legend := "ID\tREC\tSIZE"

		find := ParseIndex(indx)

		PartitionPattern(topPattern, star, rdr,
			func(str string) {
				recordCount++

				id := FindIdentifier(str[:], parent, find)
				if id == "" {
					return
				}
				if lbl != "" {
					fmt.Printf("%s\t%d\t%s\n", id, len(str), lbl)
				} else {
					fmt.Printf("%s\t%d\n", id, len(str))
				}
			})

		if timr {
			printDuration("records")
		}

		return
	}

	// SORT XML RECORDS BY IDENTIFIER

	// -pattern record_name -sort parent/element@attribute^version, strictly alphabetic sort order (undocumented)
	if len(args) == 4 && args[2] == "-sort" {

		indx := args[3]

		// create map that records each UID
		order := make(map[string][]string)

		find := ParseIndex(indx)

		PartitionPattern(topPattern, star, rdr,
			func(str string) {
				recordCount++

				id := FindIdentifier(str[:], parent, find)
				if id == "" {
					return
				}

				data, ok := order[id]
				if !ok {
					data = make([]string, 0, 1)
				}
				data = append(data, str)
				// always need to update order, since data may be reallocated
				order[id] = data
			})

		var keys []string
		for ky := range order {
			keys = append(keys, ky)
		}
		// sort fields in alphabetical order, unlike xtract version, which sorts numbers by numeric value
		sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		for _, id := range keys {

			strs := order[id]
			for _, str := range strs {
				os.Stdout.WriteString(str)
				os.Stdout.WriteString("\n")
			}
		}

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// REPORT UNRECOGNIZED COMMAND

	fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized rchive command\n")
	os.Exit(1)
}
