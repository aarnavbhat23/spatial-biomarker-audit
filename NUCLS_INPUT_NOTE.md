# Non-annotation sentinel clarification

The first NuCLS run stopped before computing spatial outcomes because the expert
columns also contain `DidNotAnnotateFOV`. This is a source non-annotation sentinel,
not a biological class. Treat it as missing, consistent with the frozen exclusion
of missing expert annotations. Keep all biological mappings, field eligibility,
graphs, expert comparisons and patient aggregation unchanged. This explicit input
clarification precedes the first evaluable result; do not describe it as an
unchanged preregistration. No class frequencies or outcomes informed this choice.
