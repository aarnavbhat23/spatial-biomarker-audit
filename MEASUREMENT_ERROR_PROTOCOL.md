# Marker-call robustness screen

Exploratory development sensitivity analysis on CRC02/CRC11/CRC19 only. Source
calls are the reference for injected errors, not independently verified truth.
Keep acquisition policy, batch 16, 1024-pixel hypothetical fields, joint alpha
.05 and widths .005 frozen. No new confirmation patients are used.

Derive a label-blind graph-degree enlargement of the existing confidence intervals
for externally supplied bounds on the number of wrong marker calls. Separately,
derive sharp full-census target bounds when at most b PanCK calls are wrong and
CD3e calls are exact. That special case tests whether full acquisition could
resolve the biological target under the stated error budget.

For reference data, examine b=ceil(n*f) at f=.0001/.001/.005/.01 for each of the
three radius claims. Then inject PanCK-only errors at fractions .001 and .005,
either uniformly at random or concentrated to maximize the absolute radius-80
target shift. Apply the same flipped cell set to all three claims. Five acquisition
repeats per patient, fraction and corruption type. Compare coverage of recorded
corrupted-call targets with coverage of reference targets, and with geometry-only
error-enlarged intervals. Record width, bias and source-call coverage separately.

The adversarial one-marker flip set is exact for the radius-80 objective when
CD3e is fixed; its effect at other radii is evaluated without re-selection.
Do not use unseen marker labels to set online interval margins. Geometry-only
degree bounds are known before acquisition; label-weighted sharp bounds are a
full-census diagnostic. The assumed error budgets are sensitivity parameters,
not estimates of actual assay error. Seed 20261114 is an arbitrary identifier.
