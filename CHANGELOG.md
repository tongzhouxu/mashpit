# Changelog

## 1.0.0 — 2026-10-29

- Select adaptive representatives from NCBI Pathogen Detection SNP trees with optional radius control, bounded assembly retries, and representative reselection.
- Build custom databases from local FASTA files, annotate metadata with change history, and shard or reshard signatures for parallel local queries.
- Merge compatible custom or taxon snapshots atomically, with strict sketch/settings checks, metadata conflict reports, provenance, and rebuilt cluster summaries.
- Preserve IsoScout `REP_CLUSTER` annotations during taxon merges, reject conflicting annotations and unsupported SQLite tables, and archive input package manifests without reusing stale checksums.
- Retain exact-release metadata and complete cluster membership, including isolates without assemblies. Read bundled versioned TSVs from existing IsoScout `databases-v2` packages.
- Add standalone offline HTML reports and enriched Streamlit cluster reports with full membership statistics, dates, serotypes, and downloadable evidence.
- Document Mashpit's screening role alongside IsoScout's downstream SNP workflow, Python 3.8–3.11 support, and Git/venv installation on ARM64. Raspberry Pi hardware performance remains unvalidated.
- Include source-distribution test fixtures and a real-data IsoScout compatibility record; verify wheel installation and tag/version consistency before publishing.

- Constrain PyArrow to versions below 26 for compatibility with Mashpit's NumPy 1.24 dependency, fixing Streamlit report failures on fresh Python 3.11 installations.
- Add `mashpit query --no-tree` and the equivalent Streamlit setting. Search scores, ranks, ties, candidate counts, and cluster summaries are unchanged.
- Write required search CSVs before optional tree retrieval, pairwise distances, construction, and rendering. Optional exceptions now warn and preserve a successful query exit; required search/output failures remain fatal.
- Add `<query>_tree_status.json` with `generated`, `disabled`, `skipped_insufficient_hits`, `construction_failed`, and `rendering_failed` states, reasons, exception types, and valid artifact filenames. Retain Newick after rendering failure and remove partial/stale images; publish tree artifacts atomically.
- Replace TreeViz rendering with an iterative Matplotlib implementation. TreeViz 0.2.0's recursive `deepcopy` in `_set_uniq_innode_name` reproduced a `RecursionError` on a synthetic 201-tip ladder at Python's default recursion limit. BioPython traversal and Newick serialization also recurse; Mashpit now uses explicit stacks for these operations. The Newick parser remains BioPython's iterative parser. No global recursion-limit change is required.
- Keep Streamlit tables available after preview failure and honor nonzero CLI exits even when an earlier CSV exists.
- Extend the compact offline suite with synthetic stage-failure tests, tree-enabled/disabled search parity, and balanced 201-tip / ladder 1,201-tip annotation and rendering cases. No external data or services are used.

Large distance matrices or figures can still exceed available resources. Python exceptions in optional stages are reported in tree status, but forced process termination cannot be recovered. If the filesystem cannot write the status itself, a warning is logged; search results remain valid and tree availability must be treated as unknown.
