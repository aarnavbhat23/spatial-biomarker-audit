# Stronger greedy control, after the primary screen

The frozen exhaustive screen has now run. Primary optimal costs for E=0/E=1
are 2/8, 8/20 and 6/12; original greedy costs are 4/16, 8/24 and 7/16.
Some greedy plans retain routine reviews even after buying exact adjudications.
An advantage over that implementation alone is insufficient evidence of a
research contribution.

Add a separately labeled post-screen control: traverse the chosen greedy
actions in decreasing cost, breaking ties by increasing action index, and
remove an action whenever the remaining plan still separates every opposite
decision pair. One pass suffices: removing actions cannot make a previously
essential action redundant. Also compare to full exact measurement. Preserve
the original greedy outputs and exact optimum. This uses the complete finite
world set and says nothing about computational scaling. Freeze this rule
before evaluating its costs; do not tune deletion orders after seeing results.
