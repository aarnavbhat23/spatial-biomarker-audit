# Acquisition on a reusable camera raster

Exploratory sensitivity test, written after the cell-reuse confirmation. The
confirmed comparison treated each variable-size rectangle as an acquisition.
A real instrument may instead acquire fixed-size fields and reuse complete
previous fields. Compare logical query cores under that alternative physical model.

A logical rectangle is expanded to all aligned camera fields it intersects. It
reveals every cached target-frame cell in those fields. Each physical field is
charged once. Sampling and inference use the original logical ownership and the
same residual estimator; the extra visible cells are reused after each query.

Use hypothetical field widths 512, 1024 and 2048 native pixels, logical cores 128
and 2048, width .005, five repeats on development patients CRC02, CRC11 and CRC19,
seed 20261104+1000*patient_number+repeat. These are sensitivity values, not verified
instrument specifications. Count physical fields, logical query rounds, area and
precision. A cost model is (field_area+field_overhead)*fields + round_overhead*rounds.
Report break-even round overhead rather than assuming it is zero in deployment.
No untouched-patient confirmation is claimed for this new acquisition model.

The current manufacturer's public specifications provide whole-slide scan speed
and resolution but do not establish arbitrary ROI control, field size, caching,
or adaptive round latency. Source: https://www.rarecyte.com/orion-instruments.html
A prospective instrument benchmark is still required before claiming time savings.

## Stronger baseline extension

After the first raster screen, add frame-pair sampling before selecting a winner.
Partition target edges by the unordered pair of physical fields containing their
endpoints. Draw an unresolved group by residual PPS; acquire its one or two fields,
reuse every visible cell and every logically known edge. This avoids acquiring a
full rectangular halo when only two fields are necessary. Keep the same patient,
field-size, precision and seed grid. Retain the original raster results and append
this comparison. The method remains a classical sampling baseline; fewer fields
alone does not account for adaptive query-round latency.
