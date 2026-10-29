# Mashpit 1.0 release validation

Checked on 2026-10-29. These are software compatibility and smoke checks,
not a biological accuracy benchmark or validation of outbreak assignment.

## Automated checks

The core offline suite passes 30 tests, including compatible merges, annotation
conflicts, refusal of unsupported tables, ranking, membership summaries, and
optional-tree failures. The release workflow installs the built wheel, checks
dependencies, and verifies that its version matches the release tag before
publishing to PyPI. Source distributions include the offline test fixtures.

## IsoScout compatibility

The consumer code was cloned directly from GitHub at IsoScout commit
`8e1250d8bd6fdb6e0b625a3a06a69f1efaf7de72`. This check used its candidate
parser and REP-table reader, not an older local checkout. It does not validate
the complete Docker, raw-read, or SNP-comparison workflow, and does not update
IsoScout's pinned Mashpit runtime.

The local Cronobacter `databases-v2` package has release `PDG000000043.415`,
665 representatives, k=31, 1,000 hashes, and tree radius 15. Its SQLite,
signature, and both membership-file checksums matched `database.json`.

| Check | Observed result |
| --- | --- |
| NCBI assembly `GCA_002107775.1`, a database representative | Top cluster `PDS000017232.2`, similarity 1.0; IsoScout parsed the expected candidate |
| Full membership context | Top cluster contained 4 clinical isolates, although only 1 representative was returned for it |
| Default tree versus `--no-tree` | Both CSVs were byte-identical; tree states were `generated` and `disabled` |
| Unrelated Salmonella assembly `ABOVSX010000001.1` / `SAMN38247565` | Best similarity 0.008; IsoScout retained ambiguity and set `below_threshold=true` |
| Merge two copies of the Cronobacter package | 665 output representatives; 665 duplicate signatures removed; full cluster membership retained |
| Add a synthetic REP label to one merge input | Label survived; current IsoScout read it for `PDS000017232.2`, but not `PDS000017232.3` |
| Input package manifests | Archived under `merge_sources/`; no stale manifest presented as the output's manifest |

The Cronobacter query FASTA SHA-256 was
`de8de44aa42179d0db16e9a8184b8e733323c3a6ab5c36e63c01bb0aa3b9244d`.
It was downloaded with NCBI Datasets:

```bash
datasets download genome accession GCA_002107775.1 --include genome --filename cronobacter.zip
unzip cronobacter.zip -d cronobacter-input
mkdir query-check
cd query-check
mashpit query ../cronobacter-input/ncbi_dataset/data/GCA_002107775.1/GCA_002107775.1_ASM210777v1_genomic.fna \
  /path/to/databases/cronobacter --number 200 --threshold 0.85 --tie-tolerance-hashes 2
```

The REP label in this check was deliberately synthetic and is not a CDC
designation. Self-retrieval and an unrelated-organism check establish basic
software behavior, not sensitivity, specificity, nearest-genome recovery, or
performance on held-out isolates. Raspberry Pi hardware validation remains
pending.
