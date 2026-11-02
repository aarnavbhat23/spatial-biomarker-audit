# Published pixel-to-syndrome comparator

The Julia Project/Manifest pin the actual SeqFISH_ADCG and
SeqFISHSyndromeDecoding implementations by Jonathan A. White and collaborators.
Their original academic/noncommercial license texts are preserved in licenses/;
they are not relicensed under this repository's license. No upstream code was
modified. GLPK, not a paid Gurobi job, was used for optimization.

`inputs.json.gz` contains the fixed controlled means and the prepared public
MERFISH development crop for SYNDROME_PIXEL_PROTOCOL.md. The measured source is
Tsui et al., DOI 10.5281/zenodo.15678603, field 000, z=14; source/license details
are in datasets/spatial_v01. The synthetic subset is not independent biological
truth. `parity_check.csv` is an exact GF(2) nullspace of the source codeword span.
The source codeword order is preserved in codebook.csv.

Reproduce after installing Julia 1.11.7:

```sh
julia --project=comparators/white_pixel -e 'using Pkg; Pkg.instantiate()'
python -c 'import gzip; from pathlib import Path; p=Path("comparators/white_pixel"); (p/"inputs.json").write_bytes(gzip.decompress((p/"inputs.json.gz").read_bytes()))'
julia --project=comparators/white_pixel run_syndrome_pixel.jl comparators/white_pixel results/syndrome_pixel.json
python verify_syndrome_pixel.py
```

The workspace-local macOS ARM runtime archive used for the recorded run was
`julia-1.11.7-macaarch64.tar.gz`, SHA256
74df9d4755a7740d141b04524a631e2485da9d65065d934e024232f7ba0790b6,
verified against the official Julia version metadata. Runtime binaries and the
dependency depot are not committed. Other platforms should use the official
runtime appropriate to their architecture and the same package manifest.

The fixed measured crop hit its original fitting cap in all bits. That failure
is retained; it is not an accuracy benchmark. The separately frozen
SYNDROME_CONVERGENCE_PROTOCOL.md follows up on this computational limitation.
That extended run did not complete a measured decode; its failed outcome and
sanitized log are retained separately. It does not supply a second point table.
To reproduce its attempted configuration, decompress convergence_inputs.json.gz
instead of inputs.json.gz and keep its output separate. The runner now saves
a .checkpoint.json file after each fitted bit.
