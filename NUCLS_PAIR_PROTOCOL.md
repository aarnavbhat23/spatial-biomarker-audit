# Coverage-driven two-expert extension, before pairwise outcomes

The three-senior-expert screen produced zero eligible fields: SP.1 contains
DidNotAnnotateFOV for all 1,569 anchors. Its result is unevaluable, not a successful
zero-disagreement finding. Coverage inspection (without spatial outcomes) shows
SP.2/SP.3 share 622 detected labeled anchors, with seven fields meeting the
unchanged minimum of 20. The other two senior-expert pairs have no common cells.

Freeze a separate external reference-choice extension using SP.2/SP.3 only.
Keep the source, biological mapping, shared bounding-box geometry, kNN6, field
eligibility, exact patient weighting, all five patients and pairwise compatibility
checks from NUCLS_REFERENCE_PROTOCOL.md. Select no fields or patients by sign,
disagreement or effect size. The observed expert coverage entirely determines
this choice. Preserve the failed three-expert artifact unchanged; save a separate
nucls_pair_references.json with this protocol hash and stage clearly identified.

This extension is not confirmation of a successful three-expert protocol and
not validation of original protein-marker claims. Observed-reference stability
is only stability over these two experts on their commonly labeled population.
It need not hold for an absent expert, undetected cells or biological truth.
A failed compatibility triangle demonstrates that both experts cannot be within
the specified error allowance of the same truth; it identifies neither as wrong.
