"""A classical two-population lower bound with realizable marker assignments."""
import json
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parent


def main():
    report={}
    for patient in ('CRC02','CRC11','CRC19'):
        d=dict(np.load(ROOT/'data'/f'{patient}-prepared.npz'))
        xy=d['coords'];edges=d['edges'];rectangles=d['rectangles'];population=len(rectangles)
        candidates=[]
        for tile,(x0,y0,x1,y1) in enumerate(rectangles):
            # Cells 80 pixels inside the unexpanded core cannot be observed
            # through any other coarse region's 80-pixel measurement halo.
            deep=(xy[:,0]>=x0+160)&(xy[:,0]<x1-160)&(xy[:,1]>=y0+160)&(xy[:,1]<y1-160)
            visible_elsewhere=False
            for other,(a,b,c,e) in enumerate(rectangles):
                if other==tile:continue
                if np.any(deep&(xy[:,0]>=a)&(xy[:,0]<c)&(xy[:,1]>=b)&(xy[:,1]<e)):
                    visible_elsewhere=True
            if visible_elsewhere:raise AssertionError('Hotspot leaks into another measurement region')
            positive=deep[edges[:,0]]&deep[edges[:,1]]
            candidates.append((int(positive.sum()),tile,int(deep.sum())))
        count,tile,cells=max(candidates)
        delta=count/len(edges)
        report[patient]={'tiles':population,'hotspot_tile':tile,'changed_cells':cells,
                         'realizable_positive_directed_pairs':count,'readout_separation':delta,
                         'cases':[{'sampled_tiles':size,'indistinguishable_probability':1-size/population,
                                   'minimax_mse_lower_bound':delta**2*(1-size/population)/4,
                                   'maximum_worst_population_coverage_if_always_shorter_than_separation':1-(1-size/population)/2}
                                  for size in (4,8) if size<=population]}
    (ROOT/'results/hotspot_information_bound.json').write_text(json.dumps({
      'scope':'Classical two-point bound under uniform coarse-region sampling; synthetic marker assignments on real fixed geometry.',
      'construction':'All markers zero versus both markers positive only in one physically unobserved deep core. Not a biological coexpression claim.',
      'patients':report},indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
