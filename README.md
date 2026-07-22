# Mashpit

[![DOI](https://joss.theoj.org/papers/10.21105/joss.07306/status.svg)](https://doi.org/10.21105/joss.07306)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/mashpit/README.html)
![unittest](https://github.com/tongzhouxu/mashpit/actions/workflows/python-app.yml/badge.svg)
[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![PyPI release](https://img.shields.io/pypi/v/mashpit)](https://pypi.python.org/pypi/mashpit/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

**Mashpit** is a lightweight, local genomic epidemiology platform for rapidly comparing a query genome against representative genomes from NCBI Pathogen Detection SNP clusters.

Mashpit is designed for offline or resource-limited analysis, including deployment on small computers such as a Raspberry Pi. It combines compact genomic sketches, epidemiological metadata, and tree-guided representative selection to support fast preliminary genome screening without requiring a centralized computing service.

> **Mashpit 1.0 is under development.** The database-building workflow described below reflects the new tree-based design and may differ from the latest packaged release.

## Key features

- **Rapid local queries:** Compare an assembled genome against a compact database of representative pathogen genomes.
- **Tree-based representative selection:** Select representatives directly from NCBI Pathogen Detection SNP trees rather than relying on a single centroid for every cluster.
- **Adaptive database size:** Use more representatives for large or genetically diverse SNP clusters and fewer representatives for compact clusters.
- **Optional radius control:** Limit the maximum tree distance between cluster members and their selected representatives.
- **Robust assembly retrieval:** Exclude records without assembly accessions, retry failed genome downloads, and remove persistently unavailable representatives.
- **Integrated metadata:** Return epidemiological information such as isolation date, location, host, BioSample accession, assembly accession, and SNP-cluster identifiers.
- **Local visualization:** Produce a query result table and a tree showing the relationship between the query and its closest database representatives.
- **Offline operation:** After a database has been built, routine queries can be performed without uploading genomic data.

## Mashpit 1.0 database design

An NCBI Pathogen Detection taxon database contains many SNP clusters. A single representative is not sufficient for every cluster because cluster size and within-cluster diversity vary substantially.

Mashpit 1.0 therefore uses an adaptive tree-based strategy:

1. Download the selected NCBI Pathogen Detection metadata and SNP trees.
2. Remove records that do not have a valid assembly accession.
3. Determine the number of representatives required for each SNP cluster from its tree structure.
4. Select representative leaves that cover the cluster, optionally subject to a user-defined tree radius.
5. Download the selected assemblies.
6. Retry failed downloads for a limited number of attempts.
7. Remove representatives whose assemblies remain unavailable.
8. Build the sketch database and retain metadata linking every representative to its SNP cluster.

This approach reduces redundancy in compact clusters while preserving multiple genomic backgrounds in large or diverse clusters.

## Installation

### Option 1: Conda or Mamba

The latest stable release can be installed from Bioconda:

```bash
conda create -n mashpit -c conda-forge -c bioconda mashpit
conda activate mashpit
```

### Option 2: pip

Mashpit requires the NCBI Datasets command-line tool for downloading assemblies.

For Linux:

```bash
curl -o datasets \
  'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
chmod +x datasets
export PATH="$PATH:/path/to/datasets"
```

Install the stable release:

```bash
python -m pip install mashpit
```

Install the development version from source:

```bash
git clone https://github.com/tongzhouxu/mashpit.git
cd mashpit
python -m pip install .
```

Verify the installation:

```bash
mashpit --help
```

## Database contents

A Mashpit database is stored in a directory containing the sketch database, representative metadata, and build information. Depending on the Mashpit version, the directory includes files such as:

```text
DATABASE_NAME/
├── DATABASE_NAME.db
├── DATABASE_NAME.sig
├── representative metadata
└── build logs
```

Keep the entire database directory together when moving it to another computer.

## Building a database

Mashpit supports two database types:

- **Taxon database:** Build from NCBI Pathogen Detection SNP clusters for a selected organism.
- **Accession database:** Build from a user-provided list of NCBI BioSample accessions.

### Taxon database

```bash
mashpit build taxon salmonella --species Salmonella
```

A specific NCBI Pathogen Detection release can be selected with `--pd_version`:

```bash
mashpit build taxon listeria_innocua \
  --species Listeria_innocua \
  --pd_version PDG000000091.9
```

Supported taxon directory names are listed in the [NCBI Pathogen Detection results directory](https://ftp.ncbi.nlm.nih.gov/pathogen/Results/).

### Accession database

Prepare a text file containing one BioSample accession per line:

```text
SAMN00000001
SAMN00000002
SAMN00000003
```

Then run:

```bash
mashpit build accession custom_database --list biosamples.txt
```

### Tree-radius option

Mashpit 1.0 can optionally use a tree-radius threshold during representative selection. A smaller radius generally selects more representatives and gives denser coverage of within-cluster diversity. A larger radius selects fewer representatives and produces a smaller database.

Use the radius option only after evaluating an appropriate value for the organism and intended surveillance resolution.

### Assembly download handling

Before representative selection and download, Mashpit removes rows without valid assembly accessions. Selected assemblies that fail to download are retried automatically. After the configured retry limit, unavailable assemblies are removed from the representative set and recorded in the build log.

This prevents missing assemblies from producing invalid database entries.

## Querying a database

Run a query using an assembled genome in FASTA format:

```bash
mashpit query sample.fasta path/to/database
```

Common query options include:

```text
--number       maximum number of database matches returned
--threshold    minimum similarity used when constructing the local result tree
--annotation   metadata field used to annotate tree tips
```

Example:

```bash
mashpit query sample.fasta path/to/salmonella \
  --number 200 \
  --threshold 0.85 \
  --annotation isolation_source
```

## Query output

A query produces:

1. **Result table**

   A CSV file containing the closest representatives and associated metadata. Results include fields such as BioSample accession, assembly accession, SNP-cluster accession, similarity score, and NCBI links.

2. **Local tree**

   Newick and image files showing the query genome together with its closest database representatives.

3. **Log file**

   A timestamped record of the query parameters and processing steps.

Mashpit results are intended for rapid screening and prioritization. They do not replace validated SNP-pipeline, cgMLST, or whole-genome phylogenetic analyses when formal outbreak confirmation is required.

## Graphical interface

Start the local Mashpit interface with:

```bash
mashpit gui
```

This launches a Streamlit app and opens it in your browser, by default at:

```text
http://localhost:8501
```

Use `--port` to run on a different port:

```bash
mashpit gui --port 8888
```

The interface lets users select a local database, upload a query assembly, run the search, and inspect the result table and tree - all without the genome or results leaving the local machine.

A pre-built Mashpit database is required.

## Intended use

Mashpit is intended for fast, lightweight genomic screening, including:

- preliminary placement of newly sequenced isolates;
- identification of the most relevant NCBI Pathogen Detection SNP clusters;
- local analysis of sensitive genomic data;
- decentralized surveillance in laboratories with limited computing infrastructure;
- prioritization of isolates for higher-resolution analysis.

Mashpit should not be used as the sole evidence for declaring or excluding an outbreak relationship.

## Citation

Xu T. et al. (2024). Mashpit: sketching out genomic epidemiology.  
*Journal of Open Source Software*, 9(104), 7306.  
https://doi.org/10.21105/joss.07306

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development and contribution guidelines.

## License

Mashpit is distributed under the GNU General Public License v1.0. See [LICENSE](LICENSE).
