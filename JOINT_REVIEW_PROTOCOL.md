# Online acquisition/review prototype protocol

Development-only screen, locally written before outcomes. One uncertain marker
A, exact B and a fixed graph. Full-panel B observation is a common prerequisite;
its cost is excluded from incremental comparisons and must be charged in any
end-to-end claim. Reference A calls supply injected-error truth only. No lab,
independent adjudication, vendor prediction or external cohort is evaluated.

For each development patient CRC02/CRC11/CRC19, select the first eight occupied
1024-pixel fields in lexicographic geometry order. Include all directed source
graph pairs whose A endpoint lies in those fields, retaining B halo endpoints.
Use the existing three radii and a joint deterministic target width .005.
Weights count B-positive neighbors, known through the prerequisite B census.
Inject ceil(.005*n_A) A errors, either random (seed 20261115+patient number) or
the largest absolute radius-80 shift. Use one table per patient/corruption.

Compare one-step minimax width reduction per cost, full acquisition followed
by influence-based review, and full acquisition followed by random review.
Each new field costs 1+.1 feedback units; independent review costs .05 or .5.
These are arbitrary field-equivalent scenarios, not real instrument timings.
Review is an ideal exact oracle. Record intervals at every decision, attainment,
fields, reviews, hypothetical cost and software runtime. All six tables and
36 policy/cost scenarios must remain visible. Run one scenario through localhost
HTTP to verify the incremental observation interface. Bulk comparisons use the
same observation-service interface without network serialization overhead.

Continue the proposed interleaving policy only if it beats the best fixed-stage
control by at least 10% in at least four of six patient/corruption tables for
one review-cost scenario, with all intervals covering the reference targets.
Otherwise retain the valid interval/oracle and stop this policy. A successful
prototype still needs uncertain B, review error, real costs and fresh data.

Implementation correction: the initial run charged field feedback in reported
costs but omitted it from the policy ranking denominator. Correct the ranking
to divide field gain by 1.1, retain the initial results as an audit file and
rerun the unchanged protocol. This is a formula correction, with no retuning
of costs, precision, panels or success criteria.
