---
title: 'Mashpit: sketching out genomic epidemiology'
tags:
  - Python
  - Min-hash 
  - Mash
  - Sourmash
  - Outbreak
authors:
  - name: Tongzhou Xu
    affiliation: 1
    orcid: 0000-0002-3829-8688
  - name: Henk C. den Bakker
    orcid: 0000-0002-4086-1580
    affiliation: 1
  - name: Xiangyu Deng
    orcid: 0000-0002-7251-2529
    affiliation: 1
  - name: Lee S. Katz
    affiliation: "1, 2"
    orcid: 0000-0002-2533-9161
affiliations:
 - name: Center for Food Safety, University of Georgia, Griffin, GA, USA
   index: 1
 - name: Enteric Diseases Laboratory Branch (EDLB), Centers for Disease Control and Prevention, Atlanta, GA, USA
   index: 2
date: 13 August 2017
bibliography: paper.bib

---

# Summary

We are in the era of genomic epidemiology.
For many transmissible diseases, large percentages of the pathogenic agent are being whole genome sequenced.
Then, their sequences are being compared in various ways.
One example is SARS-CoV-2, where the virus is being sequenced and then shared through either Global Initiative on Sharing Avian Influenza Data (GISAID) or National Institute for Biotechnology Information (NCBI).
Next, these sequences are downloaded and compared in various platforms such as NextStrain [@hadfield2018nextstrain].
From these comparisons, scientists can help track down transmission of the virus or describe its population structure.
Another example is _Salmonella_, where it is sequenced through the PulseNet molecular surveillance network and then uploaded to NCBI [@nadon2017pulsenet].
From there, these sequences can be analyzed in a variety of ways including multilocus sequence typing (MLST).
These comparisons can help track down a food vehicle and lead to actionable results
such as a food recall or inspection.

Amazingly for these organisms, at the time of this writing in January 2022, there are more than 3 million SARS-CoV-2 genomes and more than 500 thousand _Salmonella_ genomes.
These numbers are expected to increase dramatically and therefore faster methods are needed.

There have been some major advances to increase the speed of these comparisons.
We observed that a new algorithm for genomics called Min-Hash became available
to independent researchers.
A very commonly used software for Min-Hash is called Mash [@ondov2016mash].
Querying with Mash can be about 3 orders of magnitude faster than other common methods like
Basic Local Alignment Search Tool (BLAST)
and can have a smaller disk footprint.
Therefore it can be run on more common scientific workstations.
However, it does not yield metadata (e.g., date of isolation)
which is necessary to gain meaning (e.g., if it is an ongoing outbreak or if it happened 20 years ago).
On the other hand, NCBI has also created the Pathogen Detection (PD) pipeline and site.
It combines information from at least three databases: SRA, GenBank, and BioSample.
About once a day, it compares all genomes of a given taxon, separates all genomes into individual clusters using MLST, and then creates a phylogeny for each cluster using single nucleotide polymorphism (SNP) analysis.
This method is quite comprehensive but it relies on each sample being public, and it cannot be executed locally.

We present Mashpit, a new rapid genomic epidemiology platform to query against these large groups of genomes on a local computer.

# Statement of need 

Querying a sample against these magnitudes of genomes is becoming less sustainable,
especially for smaller laboratories.
Currently, GISAID and NCBI are staying ahead of the curve by producing a global tree of each organism
every day.
This is done through hurculean efforts, cutting edge algorithms, and powerful computers.
However, smaller laboratories usually have a scientific workstation or similar,
much different than a cluster computing system.

We note that for some organisms like SARS-CoV-2 and _Salmonella_, queries
can be of a sensitive nature.
For example, a state laboratory finding a new variant could be very sensitive
until it has been conclusively confirmed.
If it is a negative result, then the public would have been needlessly alerted
and therefore the query is sensitive.
Another example is that if a regulatory agency uncovers _Salmonella_ in a food
processing plant, then it is possible there would need to be regulatory action.
However, this result would need to be confirmed conclusively such that, again,
the public would not be needlessly alerted and a business needlessly affected.

To address any needs for speed and sensitivity, we created Mashpit.

# Mashpit design

Mashpit is comprised of three major parts: A min-hash database, its associated metadata, and the min-hash querying.

The database is created with an interface to Mash, called Sourmash (REF).
Each genome is imported by sketching it and adding it to a Sourmash signature database, which integrates the sketches into Python.
Each genome can also have an entry in the associated metadata. These data include date of isolation, geography, host age range, and other information that could be useful in an epidemiological investigation.
If the import is performed by downloading the genome from NCBI, then its associated BioSample data are also imported into the associated metadata.
We have calculated that it takes about 1-2 seconds to download each genome from NCBI, 0.5 seconds to sketch it and add it to the signature database.
For 399162 samples, it takes 2-3 minutes to import all NCBI metadata.
Finally, it takes 3-4 hours total to select a representative for each NCBI cluster.

Because it takes more than a few hours to create a large comprehensive database,
we have created a versioned Mashpit database available from our GitHub site.

With the database and its metadata complete, a user could perform a query.
The query is an assembly fasta file, which is then sketched and compared against the signature database.
The query then returns a tab delimited spreadsheet, sorted by Mash distance.
All associated metadata are included in the spreadsheet.

The speed of the query is determined by the database size (Figure X). 
_Tongzhou please finish this paragraph after finishing the figure_

# Discussion

We present Mashpit, a rapid genomic epidemiology platform.
Due to the underlying algorithm Min-Hash, it is excedingly fast.
It also has such a small hard drive and computational footprint that it can basically be used on common scientific workstations.
However, we note that the Mash distance does not correlate well to well-established distances such as MLST.
Therefore we recommend that this platform is used as a first-pass to filter unrelated samples before using a more established protocol such as MLST.
In conclusion, we believe that Mashpit is an essential genomic epidemiology tool.

# Figures

* Time to integrate one Salmonella genome (maybe)
* size of database per sample (maybe)
* speed of query per database size

# Acknowledgements

Financial support for the development of Mashpit was provided by the Center for Food Safety at the University of Georgia, USA.
The findings and conclusions in this report are those of the authors and do not necessarily represent the official position of the Centers for Disease Control and Prevention.

# References
