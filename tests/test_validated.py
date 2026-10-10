import copy
from fractions import Fraction
from itertools import product
import random
import unittest
from dynacomp.rational import matrix,mv,addv,q
from dynacomp.uniform import SparseLinearModel,_project,sharp_box_bound
from dynacomp.validated import validated_box,replay_validated,BoxKernelCache,euclidean_contraction


class ValidatedTests(unittest.TestCase):
    def model(self,a,b,c):
        return SparseLinearModel.from_spec({'n':len(a),'entries':[[i,j,str(v)] for i,row in enumerate(a) for j,v in enumerate(row) if v], 'B':b,'C':c})

    def test_random_brackets_enclose_exact_sharp_bounds(self):
        rng=random.Random(812)
        for _ in range(25):
            n=rng.randrange(2,5);a=[[Fraction(rng.randrange(-4,5),7) for j in range(n)] for i in range(n)]
            b=[[Fraction(rng.randrange(-3,4),5)] for i in range(n)];c=[[1]+[0]*(n-1)]
            m=self.model(a,b,c);v,e,g,h,d=_project(m,c)
            center=[Fraction(rng.randrange(-3,4),5) for i in range(n)];radius=[Fraction(rng.randrange(0,4),5) for i in range(n)]
            means=[[Fraction(rng.randrange(-3,4),5)] for j in range(4)];radii=[[Fraction(rng.randrange(0,4),5)] for j in range(4)]
            exact=sharp_box_bound(m,v,e,g,h,center,radius,radii,means)
            p=validated_box(m,e,g,h,d,center,radius,radii,0,input_centers=means)
            cert=p['certificate'];truth=q(exact['worst_bound'])
            self.assertLessEqual(q(cert['worst_attained_lower']),truth);self.assertGreaterEqual(q(cert['worst_upper']),truth)
            for t,row in enumerate(exact['output_bounds']):
                self.assertLessEqual(q(cert['output_attained_lower'][t][0]),q(row[0]));self.assertGreaterEqual(q(cert['output_upper'][t][0]),q(row[0]))
            w=cert['witness'];x=tuple(q(v) for v in w['x0']);z=mv(e,x)
            for u in w['inputs'][:w['time']]:x=addv(m.operator.step(x),mv(m.B,u));z=addv(mv(g,z),mv(h,u))
            actual=abs(mv(m.C,x)[w['output']]-mv(d,z)[w['output']])
            self.assertGreaterEqual(actual,q(w['absolute_error_at_least']));self.assertTrue(replay_validated(p))

    def test_arbitrary_encoder_no_inverse_needed_vertex_crosscheck(self):
        m=self.model([[Fraction(1,3),1],[0,Fraction(1,2)]],[[1],[0]],[[1,0],[0,1]])
        e,g,h,d=map(matrix,([[2,1]],[[Fraction(1,4)]],[[Fraction(1,3)]],[[1],[2]]))
        p=validated_box(m,e,g,h,d,[1,-1],[1,2],[[1],[2]],0,input_centers=[[1],[-1]])
        worst=q(0)
        for signs in product((-1,1),repeat=4):
            x=(q(1)+signs[0],q(-1)+2*signs[1]);z=mv(e,x)
            for t in range(3):
                worst=max(worst,max(abs(a-b) for a,b in zip(mv(m.C,x),mv(d,z))))
                if t<2:
                    u=(q([1,-1][t])+[1,2][t]*signs[t+2],);x=addv(m.operator.step(x),mv(m.B,u));z=addv(mv(g,z),mv(h,u))
        self.assertLessEqual(q(p['certificate']['worst_attained_lower']),worst);self.assertGreaterEqual(q(p['certificate']['worst_upper']),worst)

    def test_precision_refinement_and_reused_source(self):
        m=self.model([[Fraction(1,3)]],[[0]],[[1]])
        cache=BoxKernelCache(m,10)
        args=(m,[[1]],[[Fraction(1,3)]],[[0]],[[1]],[0],[1],[[0]]*10,'1/10000000000')
        coarse=validated_box(*args,precisions=(8,),cache=cache,slack=0)
        self.assertEqual(coarse['status'],'undecided')
        refined=validated_box(*args,precisions=(8,64),cache=cache,slack=0)
        self.assertEqual(refined['status'],'tolerance_met');self.assertEqual(cache.source_passes,2)
        validated_box(*args,precisions=(64,),cache=cache);self.assertEqual(cache.source_passes,2)

    def test_zero_horizon_and_dyadic_exactness(self):
        m=self.model([[1]],[[0]],[[1]])
        p=validated_box(m,[[2]],[[0]],[[0]],[[1]],[0],[1],[],0)
        self.assertEqual(q(p['certificate']['worst_upper']),1);self.assertEqual(q(p['certificate']['worst_attained_lower']),1)

    def test_tamper_status_witness_and_stability(self):
        m=self.model([[Fraction(1,2)]],[[0]],[[1]])
        p=validated_box(m,[[1]],[[Fraction(1,2)]],[[0]],[[1]],[0],[1],[[0]]*3,0,require_contraction=True)
        self.assertEqual(p['status'],'tolerance_met');self.assertTrue(p['reduced_contraction_certified'])
        for key in ('status','rank','reduced_contraction_certified'):
            bad=copy.deepcopy(p);bad[key]=None
            with self.assertRaises(ValueError):replay_validated(bad)
        bad=copy.deepcopy(p);bad['certificate']['witness']['x0']=['999']
        with self.assertRaises(ValueError):replay_validated(bad)
        self.assertFalse(euclidean_contraction([[1,1],[0,1]]))
        unstable=validated_box(m,[[0]],[[2]],[[0]],[[0]],[0],[0],[[0]],0,require_contraction=True)
        self.assertEqual(unstable['status'],'stability_not_certified');self.assertTrue(replay_validated(unstable))

    def test_invalid_dimensions_precision_and_cache(self):
        m=self.model([[Fraction(1,2)]],[[0]],[[1]])
        args=(m,[[1]],[[Fraction(1,2)]],[[0]],[[1]],[0],[1],[[0]],0)
        for options in ({'precisions':(True,)},{'precisions':()},{'cache':BoxKernelCache(m,2)},{'require_contraction':1}):
            with self.assertRaises(ValueError):validated_box(*args,**options)


try:
    import numpy as np
    import scipy
    NUMERICAL=True
except ImportError:NUMERICAL=False


@unittest.skipUnless(NUMERICAL,'Optional numerical discovery')
class StableProposalTests(unittest.TestCase):
    def test_metric_context_reuses_nested_basis_and_rejects_wrong_source(self):
        from dynacomp.stable_proposals import StableProposalContext,stable_observable_proposal
        m=SparseLinearModel.from_spec({'n':3,'entries':[[0,0,'1/2'],[0,1,3],[1,1,'1/3'],[1,2,2],[2,2,'1/4']],'B':[[0],[0],[1]],'C':[[1,0,0]]})
        context=StableProposalContext(m,6)
        first=context.propose(1);context.propose(3);again=context.propose(1)
        self.assertEqual(first['E'],again['E']);self.assertEqual(first['G'],again['G']);self.assertEqual(context.metric_preparations,1)
        with self.assertRaises(ValueError):stable_observable_proposal(m,7,2,context=context)

    def test_nonnormal_projection_stable_and_accuracy_certified(self):
        from dynacomp.stable_proposals import stable_observable_proposal
        from dynacomp.stability import schur_stable
        m=SparseLinearModel.from_spec({'n':2,'entries':[[0,0,'1/2'],[0,1,8],[1,1,'1/2']],'B':[[0],[1]],'C':[[1,0]]})
        result=stable_observable_proposal(m,8,2)
        self.assertTrue(result['reduced_contraction_certified']);self.assertTrue(schur_stable(result['G'])['stable'])
        p=validated_box(m,*(result[k] for k in ('E','G','H','D')),[0,0],[1,1],[[1]]*8,'1/1000000',require_contraction=True)
        self.assertEqual(p['status'],'tolerance_met');self.assertTrue(replay_validated(p))

    def test_unstable_source_does_not_receive_a_metric_proposal(self):
        from dynacomp.stable_proposals import stable_observable_proposal
        m=SparseLinearModel.from_spec({'n':1,'entries':[[0,0,2]],'B':[[1]],'C':[[1]]})
        with self.assertRaises((ValueError,np.linalg.LinAlgError)):stable_observable_proposal(m,3,1)

    def test_prepared_reuse_isolated_packet_and_shared_field(self):
        from dynacomp.validated_numeric import PreparedValidated
        m=SparseLinearModel.from_spec({'n':1,'entries':[[0,0,'1/2']],'B':[[1]],'C':[[1]]})
        p=validated_box(m,[[1]],[[Fraction(1,2)]],[[1]],[[1]],[0],[1],[[1]]*4,0,require_contraction=True)
        reduced=PreparedValidated(p);full=PreparedValidated(p,'full');p['G']=[['999']]
        initial=np.array([.5]);inputs=np.ones((4,1,3));expected=full.simulate(initial,inputs)
        np.testing.assert_allclose(reduced.simulate(initial,inputs),expected,atol=1e-14)
        np.testing.assert_allclose(reduced.simulate(np.full((1,3),.5),inputs),expected,atol=1e-14)
        self.assertEqual(reduced.cache_stats,{'packet_replays':1,'batches':2,'initial_fields_encoded':4})
        with self.assertRaises(ValueError):reduced.simulate([2],inputs)
        with self.assertRaises(ValueError):reduced.simulate([.5],inputs[:3])

    def test_prepared_rejects_unaccepted_plan_and_rounded_outside_endpoint(self):
        from dynacomp.validated_numeric import PreparedValidated
        m=SparseLinearModel.from_spec({'n':1,'entries':[[0,0,'1/2']],'B':[[0]],'C':[[1]]})
        p=validated_box(m,[[0]],[[0]],[[0]],[[0]],[0],[1],[],0)
        with self.assertRaises(ValueError):PreparedValidated(p)
        p=validated_box(m,[[1]],[[Fraction(1,2)]],[[0]],[[1]],[0],['1/10'],[],0)
        runner=PreparedValidated(p)
        with self.assertRaises(ValueError):runner.simulate([.1],np.empty((0,1,1)))
