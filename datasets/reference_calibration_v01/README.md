# Published cross-channel reference matches

Unmodified source table from White, Lu, Ombelets and Cai,
[Untangling overlapping barcodes](https://pmc.ncbi.nlm.nih.gov/articles/PMC12262624/),
[author repository](https://github.com/CaiGroup/UntanglingBarcodes), path
`figure_generating_notebooks/Fig2F-H/alignment/cross_channel_matched_reference_dots_pos0.csv`.
Git blob 024fe5233ece13c91db1116c9aa43bbbef791cf0 fixes the exact source bytes.
Checksums are in source_checksums.json. The acquisition script retrieves that
immutable blob and source license, not an evolving branch version.

The included SOURCE_LICENSE.md retains the authors' academic/non-commercial
reuse conditions. This dataset is not relicensed under the repository's code
license. No source analysis code has been copied into our estimator.

There are 680 selected directed matches from one position. Duplicate reference
header names must be parsed explicitly. Source raw coordinates, author-aligned
coordinates and channels are present; individual identity truth, unmatched beads,
probe efficiency and an independent specimen cohort are absent. Our benchmark
reconstructs 131 source-identity graph components, refits only on training features
and keeps all repeats of a feature in one split. It does not use author-aligned
coordinates for held-out fitting or evaluation. See CROSS_CHANNEL_PROTOCOL.md
and SUBMISSION_RESEARCH_RESULTS.md in the repository root.
