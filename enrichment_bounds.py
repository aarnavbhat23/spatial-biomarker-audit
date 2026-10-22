"""Sparse joint label-error bounds for a degree-adjusted edge covariance."""
import numpy as np
from general_marker_bounds import fast_bounds


def _top(values,budget):
    values=np.asarray(values,float)
    return float(np.sort(np.maximum(values,0))[-budget:].sum()) if budget else 0.


def _linear(weights,raw,budget):
    gain=weights*(1-2*raw)
    return [-_top(-gain,budget),_top(gain,budget)]


def _best(weights,raw,budget,maximize):
    gain=weights*(1-2*raw)*(1 if maximize else -1)
    ids=np.argsort(-gain,kind='stable')[:budget];ids=ids[gain[ids]>0]
    answer=raw.copy();answer[ids]^=1
    return answer


def contrast_bounds(edges,labels,budget_a,budget_b):
    """Conditional on a supplied simple graph and arbitrary bounded binary flips."""
    labels=np.asarray(labels);edges=np.asarray(edges);n=len(labels)
    if labels.shape!=(n,2) or not np.isin(labels,[0,1]).all():raise ValueError('Need binary two-marker labels')
    if not np.issubdtype(edges.dtype,np.integer):raise ValueError('Need integer vertex indices')
    labels=labels.astype(int);edges=edges.astype(int)
    if edges.ndim!=2 or edges.shape[1]!=2 or len(edges)==0:raise ValueError('Need nonempty directed graph')
    if edges.min()<0 or edges.max()>=n or len(np.unique(edges,axis=0))!=len(edges):
        raise ValueError('Need valid simple directed graph')
    if any(int(b)!=b or b<0 for b in (budget_a,budget_b)):raise ValueError('Need nonnegative integer budgets')
    ba,bb=min(int(budget_a),n),min(int(budget_b),n)
    u,v=edges.T;m=len(edges);a,b=labels.T.copy()
    d=np.bincount(u,minlength=n).astype(float);e=np.bincount(v,minlength=n).astype(float)
    def forward(x):return np.bincount(u,weights=x[v],minlength=n)/m-d*(e@x)/m**2
    def transpose(x):return np.bincount(v,weights=x[u],minlength=n)/m-e*(d@x)/m**2
    value=float(a@forward(b));ga,gb=forward(b),transpose(a)
    linear_a,linear_b=_linear(ga,a,ba),_linear(gb,b,bb)
    rank_cross=_top(d,ba)*_top(e,bb)/m**2
    sparse_cross=min(ba*bb,_top(d,ba),_top(e,bb))/m
    on_edge=d[u]*e[v]/m**2
    row_abs=d/m+np.bincount(u,weights=np.abs(1/m-on_edge)-on_edge,minlength=n)
    col_abs=e/m+np.bincount(v,weights=np.abs(1/m-on_edge)-on_edge,minlength=n)
    # Pad computed norms against roundoff in sparse/rank-one cancellation.
    pad=1e-12
    spectral_cross=float(np.sqrt(ba*bb*(max(0.,row_abs.max())+pad)*(max(0.,col_abs.max())+pad)))
    cross=min(sparse_cross+rank_cross,spectral_cross)
    direct=[max(-.25,value+linear_a[0]+linear_b[0]-cross-pad),
            min(.25,value+linear_a[1]+linear_b[1]+cross+pad)]
    contacts=fast_bounds(edges,labels,ba,bb)['outer_bounds']
    pa,pb=float(d@a/m),float(e@b/m)
    delta_pa,delta_pb=_linear(d/m,a,ba),_linear(e/m,b,bb)
    pa_interval=[max(0.,pa+delta_pa[0]),min(1.,pa+delta_pa[1])]
    pb_interval=[max(0.,pb+delta_pb[0]),min(1.,pb+delta_pb[1])]
    separate=[max(-.25,contacts[0]-pa_interval[1]*pb_interval[1]-pad),
              min(.25,contacts[1]-pa_interval[0]*pb_interval[0]+pad)]
    outer=[max(direct[0],separate[0]),min(direct[1],separate[1])]
    witnesses=[]
    for maximize in (False,True):
        solutions=[]
        for reverse in (False,True):
            x,y=a.copy(),b.copy()
            for _ in range(20):
                old_x,old_y=x.copy(),y.copy()
                if reverse:
                    y=_best(transpose(x),b,bb,maximize);x=_best(forward(y),a,ba,maximize)
                else:
                    x=_best(forward(y),a,ba,maximize);y=_best(transpose(x),b,bb,maximize)
                if np.array_equal(x,old_x) and np.array_equal(y,old_y):break
            solutions.append({'value':float(x@forward(y)),'a_flips':np.flatnonzero(x!=a).tolist(),
                              'b_flips':np.flatnonzero(y!=b).tolist()})
        witnesses.append((max if maximize else min)(solutions,key=lambda x:x['value']))
    attained=[x['value'] for x in witnesses]
    if not outer[0]-1e-10<=attained[0]<=attained[1]<=outer[1]+1e-10:
        raise AssertionError('Feasible witness outside outer bound')
    return {'source_contrast':value,'source_degree_abundances':[pa,pb],
            'linear_deltas':[linear_a,linear_b],'interaction_cap':cross,
            'interaction_caps':{'sparse_plus_rank_one':sparse_cross+rank_cross,'norm':spectral_cross},
            'absolute_norms':[float(row_abs.max()),float(col_abs.max())],
            'direct_bounds':direct,'separated_bounds':separate,'outer_bounds':outer,
            'attained_bounds':attained,'witnesses':witnesses}


def coupled_contrast_bounds(edges,labels,budget_a,budget_b,base=None):
    """Development tightening: combine row-specific interaction caps before sorting."""
    if base is None:base=contrast_bounds(edges,labels,budget_a,budget_b)
    labels=np.asarray(labels,int);edges=np.asarray(edges,int);n=len(labels);m=len(edges);u,v=edges.T
    a,b=labels.T;ba,bb=min(int(budget_a),n),min(int(budget_b),n)
    d=np.bincount(u,minlength=n).astype(float);e=np.bincount(v,minlength=n).astype(float)
    ga=np.bincount(u,weights=b[v],minlength=n)/m-d*(e@b)/m**2
    gb=np.bincount(v,weights=a[u],minlength=n)/m-e*(d@a)/m**2
    caps=[]
    for sigma in (-1,1):
        sparse_eligible=(a[u]==b[v]) if sigma==1 else (a[u]!=b[v])
        row_degree=np.bincount(u[sparse_eligible],minlength=n)
        col_degree=np.bincount(v[sparse_eligible],minlength=n)
        row_rank=np.array([_top(e[b==(1-bit if sigma==1 else bit)],bb) for bit in (0,1)])[a]
        col_rank=np.array([_top(d[a==(1-bit if sigma==1 else bit)],ba) for bit in (0,1)])[b]
        row_cap=np.minimum(bb,row_degree)/m+d*row_rank/m**2
        col_cap=np.minimum(ba,col_degree)/m+e*col_rank/m**2
        gain_a=sigma*(1-2*a)*ga;gain_b=sigma*(1-2*b)*gb
        caps.append(min(_top(gain_a+row_cap,ba)+_top(gain_b,bb),
                        _top(gain_a,ba)+_top(gain_b+col_cap,bb)))
    theta=base['source_contrast'];outer=[max(base['outer_bounds'][0],theta-caps[0]-1e-12),
                                      min(base['outer_bounds'][1],theta+caps[1]+1e-12)]
    return {'outer_bounds':outer,'signed_gain_caps':caps,'stage':'Post-screen development tightening'}
