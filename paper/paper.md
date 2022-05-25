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
date: 18 May 2022
bibliography: paper.bib

---

# Summary

We are in the era of genomic epidemiology.
The surveillance of many transmissible diseases is increasingly being conducted through whole genome sequencing of pathogenic agents. 
One notable example is _Salmonella_, a major foodborne pathogen routinely sequenced by surveillance programs such as PulseNet. 
Large volumes of _Salmonella_ genomes from these programs are deposited in database systems including NCBI [@nadon2017pulsenet]. 
These publicly available genomes can be analyzed in a variety of ways such as serotyping [@zhang2019seqsero2],
multilocus sequence typing (MLST) [@zhou2020enterobase], and single nucleotide polymorphism (SNP) typing [@katz2017comparative].
These analyses provide important laboratory evidence for outbreak surveillance and investigation.

At the time of this writing in March 2022, there are more than 400 thousand _Salmonella_ genomes and more than half a million other pathogen genomes at NCBI Pathogen Detection (https://www.ncbi.nlm.nih.gov/pathogens).
These numbers are expected to increase dramatically and therefore faster methods are needed.

There have been some major advances to scale up bioinformatic analyses to large volumes of pathogenic genomes.
One approach is to provide centralized resources that integrate data and analytical tools.
For example, Pathogen Detection combines information from three databases: SRA, GenBank, and BioSample.
About once a day, it compares all genomes of a given taxon, separates all genomes into individual clusters using MLST, and then creates a phylogeny for each cluster using SNP analysis.
This method is quite comprehensive, but it relies on each sample being public, and it cannot be executed locally.

Another approach is to provide new tools for decentralized and customized manipulation of genomics resources.
We observed that an algorithm for genomics called Min-Hash is well positioned for this purpose.
A commonly used software for Min-Hash is called Mash [@ondov2016mash].
Querying with Mash can be about 3 orders of magnitude faster than other common methods like Basic Local Alignment Search Tool (BLAST) and can have a smaller disk footprint [@camacho2009blast].
Therefore it can be run on more common scientific workstations.

We present Mashpit, a new rapid genomic epidemiology platform to query against these large groups of genomes on a local computer.

# Statement of need 

Querying a sample against these magnitudes of genomes is becoming less sustainable, especially for smaller laboratories.
Currently, GISAID and NCBI are staying ahead of the curve by producing a global tree of each organism every day.
This requires herculean efforts, cutting edge algorithms, and powerful computers.
However, smaller laboratories usually have a scientific workstation or similar equipment, much different than a cluster computing system.

We note that for some organisms like Salmonella, queries can be of a sensitive nature.
For example, harboring isolates in food production environments that are related to outbreak isolates is often perceived as a potential liability by food establishments, therefore thwarting the efforts to use and share the genomes of these organisms.

To address any needs for speed and sensitivity, we created Mashpit.
Mashpit queries genomes locally using Mash, thereby achieving speedy results while keeping any sensitive queries offline.

# Mashpit design

Mashpit is comprised of three major parts: A min-hash database, its associated metadata, and the min-hash querying.

The database is created with an interface to Mash, called Sourmash [@Brown2016].
Each genome is imported by sketching it and adding it to a Sourmash signature database.
Each genome can also have an entry in the associated metadata.
These data include date of isolation, geography, host age range, and other information that could be useful in an epidemiological investigation.
If the import is performed by downloading the genome from NCBI, then its associated BioSample data are also imported into the associated metadata.
We have calculated that it takes about 1-2 seconds to download each genome from NCBI, 0.5 seconds to sketch it and add it to the signature database.
For 399162 samples, it takes 2-3 minutes to import all NCBI metadata.
Finally, it takes 3-4 hours total to select a representative for each NCBI cluster.

Because it takes more than a few hours to create a large comprehensive database, we have created a versioned Salmonella Mashpit database available from our GitHub site.

With the database and its metadata complete, a user could perform a query.
The query is an assembly fasta file, which is then sketched and compared against the signature database.
The query then returns a tab delimited spreadsheet, sorted by Mash distance.
All associated metadata are included in the spreadsheet.

The speed of the query is determined by the database size (\autoref{fig:queryTime}).

# Discussion

We present Mashpit, a rapid genomic epidemiology platform.
Due to the underlying algorithm Min-Hash, it is exceedingly fast.
It also has such a small hard drive and computational footprint that it can basically be used on common scientific workstations.
However, we note that the Mash distance does not correlate well to well-established distances such as MLST.
Therefore we recommend that this platform is used as a first-pass to filter unrelated samples before using a more established protocol such as MLST.
In conclusion, we believe that Mashpit is an essential genomic epidemiology tool.

# Figures

![The speed of an individual query strongly correlates to the size of the database with a linear relationship. For every 10^4 samples, the time of each query increases about 6.8 seconds.\label{fig:queryTime}](query_time.png)

# Acknowledgements

Financial support for the development of Mashpit was provided by the Center for Food Safety at the University of Georgia, USA.
The findings and conclusions in this report are those of the authors and do not necessarily represent the official position of the Centers for Disease Control and Prevention.

# References

