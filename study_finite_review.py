"""Exact tiny planning, comparing both-marker and fallible-review assumptions."""
import json,time
from pathlib import Path
from finite_review_oracle import FiniteReviewOracle


def main():
    report={'scope':'Two-cell exact planning with assumed error budgets; hypothetical costs.','rows':[]}
    for cost in (.05,.5):
        for label,r,cap,unit in [('perfect',0,1,'events'),('one_fallible',1,1,'events'),('repeat_fallible',1,3,'events'),('systematic_cell_errors',1,3,'cells'),('systematic_with_trusted_reference',1,3,'cells')]:
            oracle=FiniteReviewOracle(2,[[(0,1)],[(0,1),(1,0)]],[(0,),(1,)],review_budgets=(r,r),review_cap=cap,review_cost=cost,review_error_unit=unit,anchor_cost=1. if label=='systematic_with_trusted_reference' else None)
            start=time.perf_counter();greedy=oracle.policy_value(oracle.initial,'greedy');census=oracle.policy_value(oracle.initial,'census_first')
            optimal=oracle.optimal_value(oracle.initial)
            def finite(x):return None if x==float('inf') else x
            row={'review_model':label,'review_cost':cost,'review_error_unit':unit,'review_error_budget_each_marker':r,'review_cap_each_cell_marker':cap,'trusted_reference_cost':oracle.anchor_cost,
                 'greedy_minimax_cost':finite(greedy),'census_first_minimax_cost':finite(census),'optimal_minimax_cost':finite(optimal),
                 'precision_attainable':optimal<float('inf'),'optimal_root_action':oracle.decisions.get(oracle.initial),
                 'solved_belief_states':oracle.optimal_value.cache_info().currsize,'seconds':time.perf_counter()-start}
            report['rows'].append(row);print(row,flush=True)
            Path('results/finite_review_study.json').write_text(json.dumps(report,indent=2)+'\n')
            for method in (oracle.feasible,oracle.bounds,oracle.children,oracle.policy_value,oracle.optimal_value):method.cache_clear()


if __name__=='__main__':main()
