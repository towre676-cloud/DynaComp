import copy
import itertools
import random
import unittest
from fractions import Fraction as Q
from dynacomp import LinearModel,uniform_enrichment,replay_uniform,sharp_box_bound
from dynacomp.uniform import SparseLinearModel
from dynacomp.rational import matrix,mv,addv,identity
from dynacomp.bounds import box_horizon_bound


def trace(model,packet,x,inputs):
    x=tuple(Q(v) for v in x)
    e,g,h,d=map(matrix,(packet[k] for k in ("E","G","H","D")))
    z=mv(e,x)
    errors=[]
    for t in range(len(inputs)+1):
        errors.append(tuple(a-b for a,b in zip(mv(model.C,x),mv(d,z))))
        if t<len(inputs):
            x=addv(mv(model.A,x),mv(model.B,inputs[t]))
            z=addv(mv(g,z),mv(h,inputs[t]))
    return errors


class UniformTests(unittest.TestCase):
    def test_sharp_bounds_equal_exhaustive_corners_random_models(self):
        rng=random.Random(9405)
        for _ in range(16):
            model=LinearModel.from_rows([[Q(rng.randrange(-3,4),4) for j in range(2)] for i in range(2)],[[1],[2]],[[1,1],[1,-1]])
            v,e,g,h=matrix([[1],[0]]),matrix([[1,0]]),matrix([[model.A[0][0]]]),matrix([[1]])
            center=[Q(1,3),-1];radius=[Q(1,2),Q(2,3)]
            inputs=[[Q(1,4)],[Q(1,2)]];means=[[1],[-2]]
            cert=sharp_box_bound(model,v,e,g,h,center,radius,inputs,means)
            packet={"E":e,"G":g,"H":h,"D":((Q(1),),(Q(1),))}
            observed=[[Q(0)]*2 for _ in range(3)]
            for signs in itertools.product((-1,1),repeat=4):
                x=[c+s*r for c,r,s in zip(center,radius,signs[:2])]
                word=[[means[t][0]+signs[2+t]*inputs[t][0]] for t in range(2)]
                for t,error in enumerate(trace(model,packet,x,word)):
                    observed[t]=[max(a,abs(b)) for a,b in zip(observed[t],error)]
            self.assertEqual(cert["output_bounds"],[[str(v) for v in row] for row in observed])

    def test_general_direction_initial_cancellation_is_retained(self):
        model=LinearModel.from_rows([["1/2",0],[0,"1/2"]],[[1],[1]],[[1,1]])
        v,e,g,h=matrix([["1/2"],["1/2"]]),matrix([[1,1]]),matrix([["1/2"]]),matrix([[2]])
        sharp=sharp_box_bound(model,v,e,g,h,[0,0],[1,1],[[1]]*3)
        legacy=box_horizon_bound(model,v,e,g,h,[0,0],[1,1],[[1]]*3)
        self.assertEqual(sharp["bound"],["0"]*4)
        self.assertGreater(Q(legacy["bound"][0]),0)

    def test_uniform_enrichment_and_attaining_witness(self):
        model=LinearModel.from_rows([["1/2",1,0],[0,"1/3",1],[0,0,"1/4"]],[[0],[0],[1]],[[1,1,0]])
        for strategy in ("prefix","worst_kernel"):
            packet=uniform_enrichment(model,[1,0,-1],[1,2,1],[[1]]*4,"1/100",input_centers=[[2]]*4,strategy=strategy)
            self.assertEqual(packet["status"],"tolerance_met")
            self.assertTrue(replay_uniform(packet))
            self.assertEqual(packet["rank"],3)

    def test_prefix_closes_finite_horizon_without_all_future_closure(self):
        model=LinearModel.from_rows([[0,1,0],[0,0,1],[0,0,0]],[[0],[0],[1]],[[1,0,0]])
        packet=uniform_enrichment(model,[0]*3,[1]*3,[[1]],0,strategy="prefix")
        self.assertEqual(packet["rank"],2)
        self.assertTrue(replay_uniform(packet))
        self.assertEqual(packet["certificate"]["bound"],["0","0"])

    def test_budget_exhaustion_has_a_real_counterexample(self):
        model=LinearModel.from_rows([[1,2],[0,1]],[[0],[1]],[[1,0]])
        packet=uniform_enrichment(model,[0,0],[1,1],[[1]]*2,0,max_rank=1)
        self.assertEqual(packet["status"],"rank_budget_exhausted")
        self.assertTrue(replay_uniform(packet))
        witness=packet["certificate"]["attaining_witness"]
        errors=trace(model,packet,witness["x0"],witness["inputs"])
        self.assertEqual(abs(errors[witness["time"]][witness["output"]]),Q(packet["certificate"]["worst_bound"]))

    def test_signed_centers_negative_attaining_witness(self):
        model=LinearModel.from_rows([[1,1],[0,1]],[[0],[1]],[[1,0]])
        packet=uniform_enrichment(model,[0,-3],[0,1],[[0],[1]],0,input_centers=[[-2],[-2]],max_rank=1)
        self.assertLess(Q(packet["certificate"]["attaining_witness"]["signed_error"]),0)
        self.assertTrue(replay_uniform(packet))

    def test_tampered_packets_and_false_status_rejected(self):
        model=LinearModel.from_rows([[1,1],[0,1]],[[0],[1]],[[1,0]])
        original=uniform_enrichment(model,[0,0],[1,1],[[1]]*2,0,max_rank=1)
        mutations=[lambda p:p.__setitem__('status','tolerance_met'),lambda p:p['G'][0].__setitem__(0,'2'),lambda p:p['certificate']['bound'].__setitem__(1,'0'),lambda p:p['certificate']['attaining_witness']['x0'].__setitem__(1,'0'),lambda p:p['initial_radius'].__setitem__(1,'2'),lambda p:p.__setitem__('rank',2)]
        for mutate in mutations:
            p=copy.deepcopy(original);mutate(p)
            with self.assertRaises(ValueError):replay_uniform(p)

    def test_invalid_domains(self):
        model=LinearModel.from_rows([[1]],[[1]],[[1]])
        for kwargs in ({'tolerance':-1},{'max_rank':0},{'input_centers':[[0],[0]]},{'strategy':'bad'}):
            args={'tolerance':0};args.update(kwargs)
            with self.assertRaises(ValueError):uniform_enrichment(model,[0],[1],[[1]],**args)
        with self.assertRaises(ValueError):uniform_enrichment(model,[0],[-1],[[1]],0)

    def test_sparse_and_dense_agree_random_completion(self):
        rng=random.Random(510)
        for _ in range(12):
            model=LinearModel.from_rows([[Q(rng.randrange(-2,3),5) for j in range(3)] for i in range(3)],[[1],[0],[1]],[[1,0,0]])
            entries=[[i,j,str(v)] for i,row in enumerate(model.A) for j,v in enumerate(row) if v]
            sparse=SparseLinearModel.from_spec({'n':3,'entries':entries,'B':model.B,'C':model.C})
            p=uniform_enrichment(model,[0]*3,[1]*3,[[1]]*3,0)
            s=uniform_enrichment(sparse,[0]*3,[1]*3,[[1]]*3,0)
            self.assertEqual(p,s)
            self.assertTrue(replay_uniform(p))

    def test_zero_horizon_target_span_is_exact(self):
        model=LinearModel.from_rows([[1,1],[0,1]],[[0],[1]],[[1,1]])
        p=uniform_enrichment(model,[0,0],[1,1],[],0)
        self.assertEqual(p["rank"],1)
        self.assertTrue(replay_uniform(p))

    def test_full_output_ranges_equal_brute_force_corners(self):
        from dynacomp.uniform import output_box_range
        model=LinearModel.from_rows([[1,-1],[0,1]],[[1],[2]],[[1,1],[1,-1]])
        result=output_box_range(model,[2,-1],[1,2],[[1],[2]],[[3],[-2]])
        observed=[[] for _ in range(3)]
        for signs in itertools.product((-1,1),repeat=4):
            x=(Q(2+signs[0]),Q(-1+2*signs[1]))
            word=[[3+signs[2]],[-2+2*signs[3]]]
            for t in range(3):
                observed[t].append(mv(model.C,x))
                if t<2:x=addv(mv(model.A,x),mv(model.B,word[t]))
        for t in range(3):
            for k in range(2):
                self.assertEqual(Q(result['lower'][t][k]),min(v[k] for v in observed[t]))
                self.assertEqual(Q(result['upper'][t][k]),max(v[k] for v in observed[t]))
