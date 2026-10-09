# Changelog

## Unreleased

- Add `mashpit query --no-tree` and the equivalent Streamlit setting. Search scores, ranks, ties, candidate counts, and cluster summaries are unchanged.
- Write required search CSVs before optional tree retrieval, pairwise distances, construction, and rendering. Optional exceptions now warn and preserve a successful query exit; required search/output failures remain fatal.
- Add `<query>_tree_status.json` with `generated`, `disabled`, `skipped_insufficient_hits`, `construction_failed`, and `rendering_failed` states, reasons, exception types, and valid artifact filenames. Retain Newick after rendering failure and remove partial/stale images; publish tree artifacts atomically.
- Replace TreeViz rendering with an iterative Matplotlib implementation. TreeViz 0.2.0's recursive `deepcopy` in `_set_uniq_innode_name` reproduced a `RecursionError` on a synthetic 201-tip ladder at Python's default recursion limit. BioPython traversal and Newick serialization also recurse; Mashpit now uses explicit stacks for these operations. The Newick parser remains BioPython's iterative parser. No global recursion-limit change is required.
- Keep Streamlit tables available after preview failure and honor nonzero CLI exits even when an earlier CSV exists.
- Extend the compact offline suite with synthetic stage-failure tests, tree-enabled/disabled search parity, and balanced 201-tip / ladder 1,201-tip annotation and rendering cases. No external data or services are used.

Large distance matrices or figures can still exceed available resources. Python exceptions in optional stages are reported in tree status, but forced process termination cannot be recovered. If the filesystem cannot write the status itself, a warning is logged; search results remain valid and tree availability must be treated as unknown.
