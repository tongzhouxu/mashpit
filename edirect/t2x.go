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
// File Name:  t2x.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

/*
  Compile application by running:

  go build t2x.go
*/

package main

import (
	"bufio"
	"fmt"
	"html"
	"io"
	"os"
	"runtime"
	"runtime/debug"
	"strconv"
	"strings"
	"time"
)

const t2xHelp = `
Data Files

  -input     Read JSON from file instead of stdin
  -output    Write XML to file instead of stdout

Parent Object Names

  -set       Replace set wrapper
  -rec       Replace record wrapper

Ignore Heading Line

  -skip      Lines to skip

Performance Measurement

  -timer     Report processing duration and rate

Example

  nquire -ftp ftp.ncbi.nlm.nih.gov gene/DATA gene_info.gz |
  gunzip -c | grep -v NEWENTRY | cut -f 2,3 |
  t2x -set Set -rec Rec -skip 1 Code Name

`

// TableConverter parses tab-delimited files into XML
func TableConverter(inp io.Reader, out io.Writer, args []string) int {

	if inp == nil {
		return 0
	}

	head := ""
	tail := ""

	hd := ""
	tl := ""

	skip := 0
	lower := false
	upper := false
	indent := true

	var fields []string
	numFlds := 0

	for len(args) > 0 {
		str := args[0]
		switch str {
		case "-help":
			fmt.Printf("t2x\n%s\n", t2xHelp)
			return 0
		case "-set":
			args = args[1:]
			if len(args) < 1 {
				fmt.Fprintf(os.Stderr, "\nERROR: No argument after -set\n")
				os.Exit(1)
			}
			set := args[0]
			if set != "" && set != "-" {
				head = "<" + set + ">"
				tail = "</" + set + ">"
			}
			args = args[1:]
		case "-rec":
			args = args[1:]
			if len(args) < 1 {
				fmt.Fprintf(os.Stderr, "\nERROR: No argument after -rec\n")
				os.Exit(1)
			}
			rec := args[0]
			if rec != "" && rec != "-" {
				hd = "<" + rec + ">"
				tl = "</" + rec + ">"
			}
			args = args[1:]
		case "-skip":
			args = args[1:]
			if len(args) < 1 {
				fmt.Fprintf(os.Stderr, "\nERROR: No argument after -skip\n")
				os.Exit(1)
			}
			tmp := args[0]
			val, err := strconv.Atoi(tmp)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: -skip argument (%s) is not an integer\n", tmp)
				os.Exit(1)
			}
			skip = val
			args = args[1:]
		case "-lower":
			lower = true
			args = args[1:]
		case "-upper":
			upper = true
			args = args[1:]
		case "-indent":
			indent = true
			args = args[1:]
		case "-flush":
			indent = false
			args = args[1:]
		default:
			// remaining arguments are names for columns
			fields = append(fields, str)
			numFlds++
			args = args[1:]
		}
	}

	if numFlds < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient arguments for t2x\n")
		os.Exit(1)
	}

	var buffer strings.Builder
	count := 0
	recordCount := 0
	okay := false
	row := 0

	wrtr := bufio.NewWriter(out)

	scanr := bufio.NewScanner(inp)

	if head != "" {
		buffer.WriteString(head)
		buffer.WriteString("\n")
	}

	for scanr.Scan() {

		line := scanr.Text()

		row++

		if skip > 0 {
			skip--
			continue
		}

		cols := strings.Split(line, "\t")

		if len(cols) != numFlds {
			fmt.Fprintf(os.Stderr, "Mismatched columns in row %d - '%s'\n", row, line)
			continue
		}

		if hd != "" {
			if indent {
				buffer.WriteString("  ")
			}
			buffer.WriteString(hd)
			buffer.WriteString("\n")
		}

		for i, fld := range fields {
			val := cols[i]
			if lower {
				val = strings.ToLower(val)
			}
			if upper {
				val = strings.ToUpper(val)
			}
			val = html.EscapeString(val)
			val = strings.TrimSpace(val)
			if indent {
				buffer.WriteString("    ")
			}
			buffer.WriteString("<")
			buffer.WriteString(fld)
			buffer.WriteString(">")
			buffer.WriteString(val)
			buffer.WriteString("</")
			buffer.WriteString(fld)
			buffer.WriteString(">")
			buffer.WriteString("\n")
		}

		if tl != "" {
			if indent {
				buffer.WriteString("  ")
			}
			buffer.WriteString(tl)
			buffer.WriteString("\n")
		}

		count++
		recordCount++

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

	if tail != "" {
		buffer.WriteString(tail)
		buffer.WriteString("\n")
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

func main() {

	// skip past executable name
	args := os.Args[1:]

	goOn := true

	timr := false

	infile := ""
	outfile := ""

	for len(args) > 0 && goOn {
		str := args[0]
		switch str {
		case "-help":
			fmt.Printf("t2x\n%s\n", t2xHelp)
			return
		case "-i", "-input":
			// read data from file instead of stdin
			args = args[1:]
			if len(args) < 1 {
				fmt.Fprintf(os.Stderr, "Input file name is missing\n")
				os.Exit(1)
			}
			infile = args[0]
			if infile == "-" {
				infile = ""
			}
			args = args[1:]
		case "-o", "-output":
			// write data to file instead of stdout
			args = args[1:]
			if len(args) < 1 {
				fmt.Fprintf(os.Stderr, "Output file name is missing\n")
				os.Exit(1)
			}
			outfile = args[0]
			if outfile == "-" {
				outfile = ""
			}
			args = args[1:]
		case "-timer":
			timr = true
			args = args[1:]
		default:
			goOn = false
		}
	}

	in := os.Stdin

	isPipe := false
	fi, err := os.Stdin.Stat()
	if err == nil {
		// check for data being piped into stdin
		isPipe = bool((fi.Mode() & os.ModeNamedPipe) != 0)
	}

	usingFile := false
	if infile != "" {

		fl, err := os.Open(infile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			os.Exit(1)
		}

		defer fl.Close()

		// use indicated file instead of stdin
		in = fl
		usingFile = true

		if isPipe && runtime.GOOS != "windows" {
			mode := fi.Mode().String()
			fmt.Fprintf(os.Stderr, "Input data from both stdin and file '%s', mode is '%s'\n", infile, mode)
			os.Exit(1)
		}
	}

	if !usingFile && !isPipe {
		fmt.Fprintf(os.Stderr, "No input data supplied\n")
		os.Exit(1)
	}

	op := os.Stdout

	if outfile != "" {

		fl, err := os.Create(outfile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			os.Exit(1)
		}

		defer fl.Close()

		// use indicated file instead of stdout
		op = fl
	}

	// initialize process timer
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
			fmt.Fprintf(os.Stderr, "\nXtract processed %.1f million %s in %.3f seconds", throughput, name, seconds)
		} else {
			fmt.Fprintf(os.Stderr, "\nXtract processed %d %s in %.3f seconds", recordCount, name, seconds)
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

	recordCount = TableConverter(in, op, args)

	// explicitly freeing memory before exit is useful for finding leaks when profiling
	debug.FreeOSMemory()

	if timr {
		printDuration("lines")
	}
}
