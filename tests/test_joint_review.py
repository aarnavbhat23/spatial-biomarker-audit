import itertools,unittest
import numpy as np
from joint_review import intervals,worst_added_ambiguity


class JointReviewTests(unittest.TestCase):
    def test_partial_intervals_are_sharp_over_all_consistent_populations(self):
        w=np.array([[1,3,2,5],[2,1,4,2]],float);totals=np.array([15,12]);recorded=np.array([0,1,1,0],bool)
        for visible in itertools.product((False,True),repeat=4):
            seen=np.array(visible)
            for known in itertools.product((False,True),repeat=4):
                reviewed=np.array(known)
                if np.any(reviewed&~seen):continue
                verified=recorded.copy()
                if reviewed.any():verified[np.flatnonzero(reviewed)[0]]^=True
                for budget in (1,2):
                    feasible=[]
                    for truth in itertools.product((False,True),repeat=4):
                        truth=np.array(truth)
                        if np.any(truth[reviewed]!=verified[reviewed]):continue
                        # Unseen recorded calls are unobserved, and can equal the true call.
                        if np.sum((truth!=recorded)&seen)>budget:continue
                        feasible.append(w@truth/totals)
                    got=intervals(w,totals,recorded,seen,reviewed,verified,budget)
                    expected=np.stack([np.min(feasible,axis=0),np.max(feasible,axis=0)],axis=1)
                    np.testing.assert_allclose(got,expected,atol=1e-15)

    def test_field_lookahead_is_exact(self):
        positive=np.array([1,4,7]);negative=np.array([2,6]);fresh=np.array([3,5,8,9])
        for b in range(5):
            values=[]
            for labels in itertools.product((False,True),repeat=len(fresh)):
                labels=np.array(labels)
                p=np.r_[positive,fresh[labels]];n=np.r_[negative,fresh[~labels]]
                values.append(sum(sorted(p,reverse=True)[:b])+sum(sorted(n,reverse=True)[:b]))
            self.assertAlmostEqual(worst_added_ambiguity(positive,negative,fresh,b),max(values))


if __name__=='__main__':unittest.main()
