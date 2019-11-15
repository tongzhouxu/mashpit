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
// File Name:  j2x.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

/*
  Download external Go libraries by running:

  cd "$GOPATH"
  go get -u github.com/gedex/inflector

  Then compile application by running:

  go build j2x.go
*/

package main

import (
	"encoding/json"
	"fmt"
	"github.com/gedex/inflector"
	"html"
	"io"
	"os"
	"runtime"
	"runtime/debug"
	"strconv"
	"strings"
	"time"
)

const j2xHelp = `
Data Files

  -input     Read JSON from file instead of stdin
  -output    Write XML to file instead of stdout

Parent Object Names

  -set       Replace set wrapper
  -rec       Replace record wrapper

Nested Array Naming

  -nest      Nested array naming policy
               [flat|recurse|plural|depth]

Performance Measurement

  -timer     Report processing duration and rate

Example

  nquire -get -url "http://mygene.info/v3" gene 2652 |
  j2x -set - -rec GeneRec -nest plural

`

// global variables, initialized (recursively) to "zero value" of type
var (
	ByteCount int
	ChanDepth int
	InBlank   [256]bool
	InElement [256]bool
)

// init function(s) run after creation of variables, before main function
func init() {

	// set communication channel buffer size
	ChanDepth = 16

	// range iterates over all elements of slice
	for i := range InBlank {
		// (would already have been zeroed at creation in this case)
		InBlank[i] = false
	}
	InBlank[' '] = true
	InBlank['\t'] = true
	InBlank['\n'] = true
	InBlank['\r'] = true
	InBlank['\f'] = true

	for i := range InElement {
		// (would already have been zeroed at creation in this case)
		InElement[i] = false
	}
	for ch := 'A'; ch <= 'Z'; ch++ {
		InElement[ch] = true
	}
	for ch := 'a'; ch <= 'z'; ch++ {
		InElement[ch] = true
	}
	for ch := '0'; ch <= '9'; ch++ {
		InElement[ch] = true
	}
	InElement['_'] = true
	InElement['-'] = true
	InElement['.'] = true
	InElement[':'] = true
}

func HasAdjacentSpacesOrNewline(str string) bool {

	whiteSpace := false

	for _, ch := range str {
		if ch == '\n' {
			return true
		}
		if ch == ' ' {
			if whiteSpace {
				return true
			}
			whiteSpace = true
		} else {
			whiteSpace = false
		}
	}

	return false
}

func CompressRunsOfSpaces(str string) string {

	whiteSpace := false
	var buffer strings.Builder

	for _, ch := range str {
		if ch < 127 && InBlank[ch] {
			if !whiteSpace {
				buffer.WriteRune(' ')
			}
			whiteSpace = true
		} else {
			buffer.WriteRune(ch)
			whiteSpace = false
		}
	}

	return buffer.String()
}

// JSONTokenizer sends sequential JSON tokens down a channel
func JSONTokenizer(inp io.Reader) <-chan string {

	if inp == nil {
		return nil
	}

	out := make(chan string, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "Unable to create JSON tokenizer channel\n")
		os.Exit(1)
	}

	tokenizeJSON := func(inp io.Reader, out chan<- string) {

		// close channel when all tokens have been sent
		defer close(out)

		// use token decoder from encoding/json package
		dec := json.NewDecoder(inp)
		if dec == nil {
			fmt.Fprintf(os.Stderr, "Unable to create JSON Decoder\n")
			os.Exit(1)
		}
		dec.UseNumber()

		for {
			t, err := dec.Token()
			if err == io.EOF {
				return
			}
			if err != nil {
				fmt.Fprintf(os.Stderr, "Unable to read JSON token '%s'\n", err)
				os.Exit(1)
			}

			// type switch performs sequential type assertions until match is found
			switch v := t.(type) {
			case json.Delim:
				// opening or closing braces (for objects) or brackets (for arrays)
				out <- string(v)
			case string:
				str := v
				if HasAdjacentSpacesOrNewline(str) {
					str = CompressRunsOfSpaces(str)
				}
				out <- str
			case json.Number:
				out <- v.String()
			case float64:
				out <- strconv.FormatFloat(v, 'f', -1, 64)
			case bool:
				if v {
					out <- "true"
				} else {
					out <- "false"
				}
			case nil:
				out <- "null"
			default:
				out <- t.(string)
			}
		}
	}

	// launch single tokenizer goroutine
	go tokenizeJSON(inp, out)

	return out
}

// JSONConverter parses JSON token stream into XML object stream
func JSONConverter(inp <-chan string, set, rec, nest string) <-chan string {

	if inp == nil {
		return nil
	}

	out := make(chan string, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "Unable to create JSON converter channel\n")
		os.Exit(1)
	}

	// opt is used for anonymous top-level objects, anon for anonymous top-level arrays
	opt := "opt"
	anon := "anon"
	if rec != "" {
		// override record delimiter
		opt = rec
		anon = rec
	}

	flat := false
	plural := false
	depth := false

	switch nest {
	case "flat":
		flat = true
	case "plural", "name":
		plural = true
	case "depth", "deep", "level":
		depth = true
	}

	// convertJSON sends XML records down a channel
	convertJSON := func(inp <-chan string, out chan<- string) {

		// close channel when all tokens have been processed
		defer close(out)

		// ensure that XML tags are legal (initial digit allowed by xtract for biological data in JSON)
		fixTag := func(tag string) string {

			if tag == "" {
				return tag
			}

			okay := true
			for _, ch := range tag {
				if !InElement[ch] {
					okay = false
				}
			}
			if okay {
				return tag
			}

			var temp strings.Builder

			// replace illegal characters with underscore
			for _, ch := range tag {
				if InElement[ch] {
					temp.WriteRune(ch)
				} else {
					temp.WriteRune('_')
				}
			}

			return temp.String()
		}

		// closure silently places local variable pointer onto inner function call stack
		var buffer strings.Builder

		// array to speed up indentation
		indentSpaces := []string{
			"",
			"  ",
			"    ",
			"      ",
			"        ",
			"          ",
			"            ",
			"              ",
			"                ",
			"                  ",
		}

		indent := 0
		if set != "" {
			indent = 1
		}

		// indent a specified number of spaces
		doIndent := func(indt int) {
			i := indt
			for i > 9 {
				buffer.WriteString("                    ")
				i -= 10
			}
			if i < 0 {
				return
			}
			buffer.WriteString(indentSpaces[i])
		}

		// recursive function definitions
		var parseObject func(tag string)
		var parseArray func(tag, pfx string, lvl int)

		// recursive descent parser uses mutual recursion
		parseValue := func(tag, pfx, tkn string, lvl int) {

			switch tkn {
			case "{":
				parseObject(tag)
				// no break needed, would use fallthrough to explicitly cause program control to flow to the next case
			case "[":
				if flat {
					parseArray(tag, pfx, lvl+1)
				} else if lvl > 0 {
					// nested JSON arrays create recursive XML objects
					doIndent(indent)
					indent++
					tg := tag
					if plural {
						tg = inflector.Pluralize(tag)
					}
					buffer.WriteString("<")
					buffer.WriteString(tg)
					buffer.WriteString(">\n")
					if depth {
						parseArray(pfx+"_"+strconv.Itoa(lvl), tag, lvl+1)
					} else {
						parseArray(tag, pfx, lvl+1)
					}
					indent--
					doIndent(indent)
					buffer.WriteString("</")
					buffer.WriteString(tg)
					buffer.WriteString(">\n")
				} else {
					parseArray(tag, pfx, lvl+1)
				}
			case "}", "]":
				// should not get here, decoder tracks nesting of braces and brackets
			case "":
				// empty value string generates self-closing object
				doIndent(indent)
				buffer.WriteString("<")
				buffer.WriteString(tag)
				buffer.WriteString("/>\n")
			default:
				// write object and contents to string builder
				doIndent(indent)
				tkn = strings.TrimSpace(tkn)
				tkn = html.EscapeString(tkn)
				buffer.WriteString("<")
				buffer.WriteString(tag)
				buffer.WriteString(">")
				buffer.WriteString(tkn)
				buffer.WriteString("</")
				buffer.WriteString(tag)
				buffer.WriteString(">\n")
			}
		}

		parseObject = func(tag string) {

			doIndent(indent)
			indent++
			buffer.WriteString("<")
			buffer.WriteString(tag)
			buffer.WriteString(">\n")

			for {
				// shadowing tag variable inside for loop does not step on value of tag argument in outer scope
				tag, ok := <-inp
				if !ok {
					break
				}

				ByteCount += len(tag)

				if tag == "}" || tag == "]" {
					break
				}

				tag = fixTag(tag)

				tkn, ok := <-inp
				if !ok {
					break
				}

				ByteCount += len(tkn)

				if tkn == "}" || tkn == "]" {
					break
				}

				parseValue(tag, tag, tkn, 0)
			}

			indent--
			doIndent(indent)
			buffer.WriteString("</")
			buffer.WriteString(tag)
			buffer.WriteString(">\n")
		}

		parseArray = func(tag, pfx string, lvl int) {

			for {
				tkn, ok := <-inp
				if !ok {
					break
				}

				ByteCount += len(tkn)

				if tkn == "}" || tkn == "]" {
					break
				}

				parseValue(tag, pfx, tkn, lvl)
			}
		}

		// process stream of catenated top-level JSON objects or arrays
		for {
			tkn, ok := <-inp
			if !ok {
				break
			}

			ByteCount += len(tkn)

			if tkn == "{" {
				parseObject(opt)
			} else if tkn == "[" {
				parseArray(anon, anon, 0)
			} else {
				break
			}

			txt := buffer.String()

			// send result through output channel
			out <- txt

			buffer.Reset()

			runtime.Gosched()
		}
	}

	// launch single converter goroutine
	// to wait until all consumers are done
	go convertJSON(inp, out)

	return out
}

func main() {

	// skip past executable name
	args := os.Args[1:]

	set := "root"
	rec := ""
	nest := ""

	timr := false

	infile := ""
	outfile := ""

	nextArg := func() (string, bool) {

		if len(args) < 1 {
			return "", false
		}

		// remove next token from slice
		nxt := args[0]
		args = args[1:]

		return nxt, true
	}

	// look for optional arguments
	for {
		arg, ok := nextArg()
		if !ok {
			break
		}

		switch arg {
		case "-help":
			fmt.Printf("j2x\n%s\n", j2xHelp)
			return
		case "-i", "-input":
			// read data from file instead of stdin
			infile, ok = nextArg()
			if !ok {
				fmt.Fprintf(os.Stderr, "Input file name is missing\n")
				os.Exit(1)
			}
			if infile == "-" {
				infile = ""
			}
		case "-o", "-output":
			// write data to file instead of stdout
			outfile, ok = nextArg()
			if !ok {
				fmt.Fprintf(os.Stderr, "Output file name is missing\n")
				os.Exit(1)
			}
			if outfile == "-" {
				outfile = ""
			}
		case "-set":
			// override set wrapper
			set, ok = nextArg()
			if ok && set == "-" {
				set = ""
			}
		case "-rec", "-record":
			// override record wrapper
			rec, ok = nextArg()
			if ok && rec == "-" {
				rec = ""
			}
		case "-nest":
			// specify nested array naming policy
			nest, ok = nextArg()
			if !ok {
				fmt.Fprintf(os.Stderr, "Nested array naming policy is missing\n")
				os.Exit(1)
			}
			if ok && nest == "-" {
				nest = "flat"
			}
			switch nest {
			case "flat", "plural", "name", "recurse", "recursive", "same", "depth", "deep", "level":
			default:
				fmt.Fprintf(os.Stderr, "Unrecognized nested array naming policy\n")
				os.Exit(1)
			}
		case "-timer":
			timr = true
		default:
			// alternative form uses positional arguments to override set and rec
			set = arg
			if set == "-" {
				set = ""
			}
			rec, ok = nextArg()
			if ok && rec == "-" {
				rec = ""
			}
			break
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

	// use output channel of tokenizer as input channel of converter
	jtkn := JSONTokenizer(in)
	jcnv := JSONConverter(jtkn, set, rec, nest)

	if jtkn == nil || jcnv == nil {
		fmt.Fprintf(os.Stderr, "Unable to create JSON to XML converter\n")
		os.Exit(1)
	}

	if set != "" {
		fmt.Fprintf(op, "<%s>\n", set)
	}

	// drain output of last channel in service chain
	// runtime assigns concurrent goroutines to execute on separate CPUs for maximum speed
	for str := range jcnv {

		if str == "" {
			continue
		}

		// send result to stdout
		op.WriteString(str)
		if !strings.HasSuffix(str, "\n") {
			op.WriteString("\n")
		}

		recordCount++

		runtime.Gosched()
	}

	if set != "" {
		fmt.Fprintf(op, "</%s>\n", set)
	}

	byteCount = ByteCount

	// explicitly freeing memory before exit is useful for finding leaks when profiling
	debug.FreeOSMemory()

	if timr {
		printDuration("records")
	}
}
