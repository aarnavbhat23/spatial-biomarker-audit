"""Exact finite-population audit with errors in both markers and repeated review."""
from functools import lru_cache
import itertools
import numpy as np


class FiniteReviewOracle:
    def __init__(self,n,claims,fields,budgets=(1,1),review_budgets=(1,1),review_cap=3,
                 field_cost=1.1,review_cost=.1,width=0.,review_error_unit='events',anchor_cost=None):
        self.n=n;self.claims=[np.asarray(e,int) for e in claims];self.fields=tuple(tuple(f) for f in fields)
        self.budgets=np.array(budgets);self.review_budgets=np.array(review_budgets);self.review_cap=review_cap
        self.anchor_cost=anchor_cost;self.review_error_unit=review_error_unit
        if review_error_unit not in ('events','cells'):raise ValueError('Unknown review error unit')
        self.field_cost=field_cost;self.review_cost=review_cost;self.width=width
        self.populations=np.array(list(itertools.product((0,1),repeat=2*n)),dtype=np.int8).reshape(-1,n,2)
        self.targets=np.stack([(self.populations[:,e[:,0],0]*self.populations[:,e[:,1],1]).mean(axis=1) for e in self.claims],axis=1)
        self.initial=((-1,)*(2*n),(0,)*(2*n),(0,)*(2*n))
        if anchor_cost is not None:self.initial=self.initial+((-1,)*(2*n),)
        self.decisions={}

    @lru_cache(None)
    def feasible(self,state):
        raw,counts,ones=(np.array(x).reshape(self.n,2) for x in state[:3])
        acquisition=((self.populations!=raw)&(raw>=0)).sum(axis=1)
        mistakes=np.where(self.populations==1,counts-ones,ones)
        review=(mistakes if self.review_error_unit=='events' else mistakes>0).sum(axis=1)
        valid=np.all(acquisition<=self.budgets,axis=1)&np.all(review<=self.review_budgets,axis=1)
        if len(state)>3:
            anchors=np.array(state[3]).reshape(self.n,2)
            valid&=np.all((anchors<0)|(self.populations==anchors),axis=(1,2))
        return tuple(np.flatnonzero(valid).tolist())

    @lru_cache(None)
    def bounds(self,state):
        ids=self.feasible(state)
        if not ids:raise ValueError('No population satisfies the stated error budgets')
        values=self.targets[list(ids)]
        return tuple(zip(values.min(axis=0).tolist(),values.max(axis=0).tolist()))

    def done(self,state):return max(hi-lo for lo,hi in self.bounds(state))<=self.width+1e-12

    def actions(self,state):
        raw,counts,_=state[:3];actions=[]
        for index,cells in enumerate(self.fields):
            if any(raw[2*c]<0 for c in cells):actions.append(('field',index))
        for i in range(2*self.n):
            anchored=len(state)>3 and state[3][i]>=0
            if raw[i]>=0 and not anchored and counts[i]<self.review_cap:
                if self.review_error_unit=='events' or counts[i]==0:actions.append(('review',i))
            if self.anchor_cost is not None and raw[i]>=0 and not anchored:actions.append(('anchor',i))
        return actions

    @lru_cache(None)
    def children(self,state,action):
        raw,counts,ones=state[:3];kind,index=action;result=[]
        outcomes=itertools.product((0,1),repeat=2*len(self.fields[index])) if kind=='field' else ((0,),(1,))
        for outcome in outcomes:
            a=list(raw);b=list(counts);c=list(ones)
            if kind=='field':
                for position,cell in enumerate(self.fields[index]):a[2*cell:2*cell+2]=outcome[2*position:2*position+2]
            elif kind=='review':b[index]+=1;c[index]+=outcome[0]
            updated=(tuple(a),tuple(b),tuple(c))
            if len(state)>3:
                anchor=list(state[3])
                if kind=='anchor':anchor[index]=outcome[0]
                updated=updated+(tuple(anchor),)
            if self.feasible(updated):result.append((tuple(outcome),updated))
        return tuple(result)

    def cost(self,action):
        if action[0]=='anchor':return self.anchor_cost
        return self.field_cost if action[0]=='field' else self.review_cost

    def greedy(self,state,census_first=False):
        width=max(hi-lo for lo,hi in self.bounds(state));choices=[]
        actions=self.actions(state)
        if census_first and any(a[0]=='field' for a in actions):actions=[a for a in actions if a[0]=='field']
        for action in actions:
            worst=max(max(hi-lo for lo,hi in self.bounds(child)) for _,child in self.children(state,action))
            gain=(width-worst)/self.cost(action)
            choices.append((gain,action[0]=='review',-state[1][action[1]] if action[0]=='review' else 0,-action[1],action))
        return max(choices)[-1] if choices else None

    @lru_cache(None)
    def policy_value(self,state,policy):
        if self.done(state):return 0.
        action=self.greedy(state,policy=='census_first')
        if action is None:return float('inf')
        return self.cost(action)+max(self.policy_value(child,policy) for _,child in self.children(state,action))

    @lru_cache(None)
    def optimal_value(self,state):
        if self.done(state):return 0.
        actions=self.actions(state)
        if not actions:return float('inf')
        best=self.policy_value(state,'greedy');chosen=self.greedy(state)
        for action in sorted(actions,key=lambda a:(a!=chosen,self.cost(a),a)):
            worst=0.
            children=sorted(self.children(state,action),key=lambda c:max(hi-lo for lo,hi in self.bounds(c[1])),reverse=True)
            for _,child in children:
                worst=max(worst,self.optimal_value(child))
                if self.cost(action)+worst>=best-1e-12:break
            value=self.cost(action)+worst
            if value<best-1e-12:best=value;chosen=action
        self.decisions[state]=chosen
        return best
