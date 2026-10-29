# Conserved molecular-position screen

Freeze on 2026-10-28 before outcome evaluation. This explores a change of
observable, not an improvement to the previous cell-type contrast. Nature
motivates conservation: reassigning detected molecules between segmentations
does not change their original coordinates or total measured mass. Detection
loss, wrong molecular identities and biased localization are separate errors.

Target: cross-channel kernel proximity S=mean K(x_A-x_B), with
K(d)=exp(-||d||^2/(2*h^2)). No cell assignments enter. This is an established
kernel statistic, not a new biological interaction or causal estimator. Report
raw proximity, not an enrichment p-value or between-cell interaction claim.

Each molecule can be displaced by at most rho. Triangle inequality gives pair
distance intervals [max(0,d-2rho),d+2rho]. Average kernel extrema for safe bounds.
This conventional pointwise enclosure is the strong control. Also report the
global Lipschitz bound S +/- 2rho*exp(-1/2)/h, clipped to [0,1]. Do not count
improvement over this weak uniform bound as algorithmic novelty.

Exact cases: three A and three B molecules, each allowed its recorded location
or one unit left/right; seeds 31,47,83; coordinates uniform in [0,20]^2 rounded
to integers; h=4,8. Enumerate all 3^6 position assignments. Compare the analytic
outer interval with exact extrema. Relabeling cell assignments must leave the
score unchanged exactly.

Scaling cases: 256 molecules per channel; seeds 101,202,303; domain [0,100]^2;
h=10,20; rho=1,3. Use two configurations for each seed: shared Gaussian cluster
centered at (50,50), SD10, and separated clusters at (30,50),(70,50), SD10.
Clip to the domain. No selection based on outcomes. Compare the two feasible
score intervals at identical h,rho. This is a controlled sensitivity screen,
not a valid null-hypothesis test or biological validation.

Continue the observable direction if exact containment holds and shared versus
separated clusters have disjoint score intervals in >=9/12 paired scaling cases.
Check one missing molecule per channel as a total-variation sensitivity and
show an explicit failure when the displacement budget is violated. Positive
results justify only a raw-coordinate public-data feasibility study. They do
not justify clinical claims, calibrated coverage, cell-type interactions or
originality; SSAM, spage2vec and spatial kernel methods are prior art.
