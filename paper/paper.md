---
title: 'Mashpit: sketching out genomic epidemiology'
tags:
  - Python
  - MinHash 
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
 - name: Center for Food Safety, University of Georgia, Griffin, GA, United States of America
   index: 1
 - name: Enteric Diseases Laboratory Branch (EDLB), Centers for Disease Control and Prevention, Atlanta, GA, United States of America
   index: 2
date: 7 Aug 2024
bibliography: paper.bib

---

# Summary

We are in the era of genomic epidemiology.
The surveillance of many transmissible diseases is increasingly being conducted through whole genome sequencing of pathogenic agents. 
One notable example is _Salmonella_, a major foodborne pathogen routinely sequenced by surveillance programs such as PulseNet [@swaminathan2001pulsenet]. 
Large volumes of _Salmonella_ genomes from these programs are deposited in database systems including NCBI [@nadon2017pulsenet]. 
These publicly available genomes can be analyzed in a variety of ways such as serotyping [@zhang2019seqsero2],
multilocus sequence typing (MLST) [@zhou2020enterobase], and single nucleotide polymorphism (SNP) typing [@katz2017comparative].
These analyses provide important laboratory evidence for outbreak surveillance and investigation.

As of August 2024, there are more than 600 thousand _Salmonella_ genomes and more than half a million other pathogen genomes at NCBI Pathogen Detection (https://www.ncbi.nlm.nih.gov/pathogens).
These numbers are expected to increase dramatically and therefore faster analysis methods are needed.

There have been some major advances to scale up bioinformatic analyses to large volumes of pathogenic genomes.
One approach is to provide centralized resources that integrate data and analytical tools.
For example, NCBI Pathogen Detection combines information from three databases: SRA, GenBank, and BioSample.
About once a day, it compares all genomes of a given taxon, separates all genomes into individual clusters using MLST, and then creates a phylogeny for each cluster using SNP analysis.
This method is quite comprehensive, but it relies on each sample being public, and it cannot be executed locally.

Another approach is to provide new tools for decentralized and customized manipulation of genomics resources.
We observed that an algorithm for genomics called MinHash is well positioned for this purpose.
A commonly used software for MinHash is called Mash [@ondov2016mash].
Querying with Mash can be about 4 orders of magnitude faster than other common methods like Basic Local Alignment Search Tool (BLAST) and can have a smaller disk footprint [@camacho2009blast;@topaz2018bmscan].
Therefore it can be run on more common scientific workstations.

We present Mashpit, a new rapid genomic epidemiology platform to query against these large groups of genomes on a local computer.

# Statement of need 

Querying a sample against these magnitudes of genomes is becoming less sustainable, especially for smaller laboratories.
Currently, GISAID and NCBI are staying ahead of the curve by producing a global tree of each organism every day [@shu2017gisaid].
This requires herculean efforts, cutting-edge algorithms, and powerful computers.
However, smaller laboratories usually have a scientific workstation or similar equipment, much different than a cluster computing system.

We note that for some organisms like _Salmonella_, queries can be of a sensitive nature.
For example, harboring isolates in food production environments that are related to outbreak isolates is often perceived as a potential liability by food establishments, therefore thwarting the efforts to use and share the genomes of these organisms.

To address any needs for speed and sensitivity, we created Mashpit.
Mashpit queries genomes locally using Mash, thereby achieving speedy results while keeping any sensitive queries offline.

# Mashpit design

Mashpit consists of three major parts: A MinHash database, its associated metadata, and the MinHash querying.

The database is created with an interface to Mash, called Sourmash [@Brown2016].
Each genome is imported by sketching it and adding it to a Sourmash signature database.
Each genome can also have an entry in the associated metadata.
These data include date of isolation, geography, host age range, and other information that could be useful in an epidemiological investigation.
Mashpit can build a species database from NCBI Pathogen Detection, termed a Mashpit taxon database, or a custom database from a list of biosample accessions. For each SNP cluster of one species on NCBI Pathogen Detection, the set of all genomes is defined as:
$$G=\{g_1,g_2,…,g_n\}$$
where n is the number of genomes in the cluster.
The centroid genome $g_c$ is calculated as:
$$g_c=\underset{g_i∈G}{argmin}\sum_{j=1}^{n} d(g_i,g_j)$$
where $d(g_i,g_j)$ is the distance between two genomes.
The centroid genome represents the most central genome in each SNP cluster, reducing redundancy while retaining representative information for queries. By default, Mashpit will download the latest SNP cluster for specified species and uses a kmer size of 31 and kmer number of 1000 for sketching the genomes. 

With the database and its metadata complete, a user could perform a query.
The query is an assembly fasta file, which is then sketched and compared against the signature database.
The query then returns a tab delimited spreadsheet, sorted by Mash distance, and a phylogenetic tree based on the Mash distance.
All associated metadata are included in the spreadsheet.

Mashpit also provides a webserver interface for users to query the database. 
The webserver is built using Flask and can be run locally or deployed on a server. 
The webserver provides a user-friendly interface for users to upload their query genomes and view the results.

# Performance
To evaluate the performance of Mashpit, we tested it on a server that runs Ubuntu 20.04.2 with an Intel Xeon CPU E5-2697 v4 2.30GHz and 256GB RAM.
We used NCBI Pathogen Detection SNP clusters that were versioned before January 2024. We then randomly selected 1000 newly added genomes for each species added to NCBI Pathogen Detection after January 2024. 
We measured the elapsed time for querying four major foodborne pathogens: _Salmonella_, _Listeria_, _E. coli_, and _Campylobacter_ (\autoref{fig:time_query}).
We also compared the query results with the true SNP cluster of the query genomes. 
We calculated the proportion of true SNP clusters appearing among the top hits at various thresholds (\autoref{fig:accuracy}). 
The 'threshold' indicates whether the correct SNP cluster is among the top 'threshold number' of query hits. 
For instance, a threshold of 25 indicates that the correct cluster is among the top 25 hits. 
Our findings indicate that _Salmonella_ achieved a 70% success rate for true clusters appearing within the top 25 hits, compared to approximately 90% for _Campylobacter_. This variability reflects differences in how species are represented in the database and the limitations of MinHash-based methods for resolving closely related clusters.

For _Salmonella_, which is the most frequently sequenced organism in NCBI Pathogen Detection, many closely related SNP clusters exist due to its extensive representation. Mash, being a MinHash-based method, operates at a resolution that is not always sufficient to distinguish fine-scale differences between these closely related clusters. 
As a result, users analyzing _Salmonella_ should interpret Mashpit results as preliminary and consider following up with higher-resolution methods for definitive SNP cluster assignments.

# Discussion

Mashpit provides a fast and lightweight platform for genomic epidemiology. 
Its MinHash-based approach enables rapid querying of large datasets on standard scientific workstations, addressing key challenges for laboratories with limited computational resources or privacy concerns.

However, we note that the Mash distance does not correlate well with established distances such as Average Nucleotide Identity (ANI) for closely related genomes [@jain2018high]. Therefore it has resolution limits when differentiating closely related clusters, particularly for species like _Salmonella_ that are highly represented in databases such as NCBI Pathogen Detection.

Therefore we recommend that this platform is used as a first-pass to filter unrelated samples before using a more established protocol such as MLST.
In conclusion, we believe that Mashpit is an essential genomic epidemiology tool.

# Figures

![Average query time for four Mashpit taxon databases. \label{fig:time_query}](figure_1.png){ width=100% }

![Probability of the true SNP cluster being included among the highest-ranking hits at varying thresholds. \label{fig:accuracy}](figure_2.png){ width=60% }

# Acknowledgements

Financial support for the development of Mashpit was provided by the Center for Food Safety at the University of Georgia, United States of America.
The findings and conclusions in this report are those of the authors and do not necessarily represent the official position of the Centers for Disease Control and Prevention.

# References

