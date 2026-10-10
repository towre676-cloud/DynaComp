import copy
import itertools
import random
import unittest
from fractions import Fraction as Q
from dynacomp.robust import switching_enrichment,replay_switching,_models,_candidate,_sharp_for_word
from dynacomp.uniform import _domain
from dynacomp.rational import matrix,mv,addv


def specs(a,b,ba=((1,),(2,)),bb=((0,),(1,))):
    return {name:{'n':len(values),'entries':[[i,j,str(v)] for i,row in enumerate(values) for j,v in enumerate(row) if v],'B':B} for name,values,B in [('a',a,ba),('b',b,bb)]}


def errors(models,packet,initial,word,inputs):
    e,d=matrix(packet['E']),matrix(packet['D']);x=tuple(Q(v) for v in initial);z=mv(e,x);out=[]
    for t in range(len(word)+1):
        out.append(tuple(a-b for a,b in zip(mv(next(iter(models.values())).C,x),mv(d,z))))
        if t<len(word):
            s,u=word[t],inputs[t];x=addv(models[s].operator.step(x),mv(models[s].B,u));z=addv(mv(matrix(packet['G'][s]),z),mv(matrix(packet['H'][s]),u))
    return out


class RobustTests(unittest.TestCase):
    def test_count_lattice_matches_all_words_and_box_corners(self):
        rng=random.Random(606)
        for _ in range(12):
            f=[[Q(rng.randrange(-2,3),4) for j in range(2)] for i in range(2)]
            a=[[Q(i==j)+f[i][j] for j in range(2)] for i in range(2)]
            b=[[Q(i==j,2)-f[i][j] for j in range(2)] for i in range(2)]
            modes=specs(a,b);C=[[1,1]];center=[1,-1];radius=[1,2];means=[[1],[-1],[2]];radii=[[1],[2],[1]]
            packet=switching_enrichment(modes,C,center,radius,radii,0,input_centers=means,max_rank=1)
            models=_models(modes,C);observed=[Q(0)]*4
            for word in itertools.product(modes,repeat=3):
                for signs in itertools.product((-1,1),repeat=5):
                    x=[c+s*r for c,r,s in zip(center,radius,signs)]
                    inputs=[[means[j][0]+signs[2+j]*radii[j][0]] for j in range(3)]
                    for t,row in enumerate(errors(models,packet,x,word,inputs)):observed[t]=max(observed[t],abs(row[0]))
            self.assertEqual(packet['certificate']['bound'],[str(x) for x in observed])
            self.assertIn(packet['certificate']['engine'],('commuting_count_lattice','exact_common_closure'))
            self.assertTrue(replay_switching(packet))

    def test_mode_order_matters_for_forcing_even_when_operators_commute(self):
        modes=specs([[1,1],[0,1]],[[1,2],[0,1]],((1,),(0,)),((0,),(1,)))
        models=_models(modes,[[1,0]]);candidate=_candidate(models,[[1,0]])
        domain=_domain(next(iter(models.values())),[0,0],[0,0],[[0],[0]],[[1],[1]])
        ab=_sharp_for_word(models,candidate,domain,('a','b'),0)[0]
        ba=_sharp_for_word(models,candidate,domain,('b','a'),0)[0]
        self.assertEqual((ab,ba),(Q(0),Q(1)))
        packet=switching_enrichment(modes,[[1,0]],[0,0],[0,0],[[0],[0]],0,input_centers=[[1],[1]],max_rank=1)
        self.assertEqual(packet['certificate']['worst_bound'],'2')
        self.assertTrue(replay_switching(packet))

    def test_noncommuting_enumeration_matches_words(self):
        modes=specs([[1,1],[0,1]],[[1,0],[1,1]])
        packet=switching_enrichment(modes,[[1,0]],[1,-1],[1,2],[[1]]*3,0,max_rank=1)
        self.assertEqual(packet['certificate']['engine'],'endpoint_word_enumeration')
        self.assertEqual(packet['certificate']['symbolic_states'],15)
        models=_models(modes,[[1,0]]);candidate=_candidate(models,packet['E']);domain=_domain(next(iter(models.values())),[1,-1],[1,2],[[1]]*3,None)
        for t,bound in enumerate(packet['certificate']['bound']):
            self.assertEqual(Q(bound),max(_sharp_for_word(models,candidate,domain,w,0)[0] for w in itertools.product(modes,repeat=t)))
        self.assertTrue(replay_switching(packet))

    def test_interior_convex_switches_are_covered(self):
        modes=specs([[1,1],[0,1]],[[1,2],[0,1]])
        p=switching_enrichment(modes,[[1,0]],[0,0],[1,1],[[1]]*4,0,max_rank=1)
        models=_models(modes,[[1,0]]);e,d=matrix(p['E']),matrix(p['D'])
        x=(Q(1),Q(-1));z=mv(e,x)
        for t,lam in enumerate((Q(1,3),Q(2,5),Q(4,7),Q(1,2))):
            u=(Q((-1)**t),)
            full=[addv(m.operator.step(x),mv(m.B,u)) for m in models.values()]
            small=[addv(mv(matrix(p['G'][s]),z),mv(matrix(p['H'][s]),u)) for s in models]
            x=tuple((1-lam)*a+lam*b for a,b in zip(*full));z=tuple((1-lam)*a+lam*b for a,b in zip(*small))
            self.assertLessEqual(abs(x[0]-mv(d,z)[0]),Q(p['certificate']['bound'][t+1]))

    def test_exact_common_closure_bypasses_long_word_budget(self):
        modes=specs([[1,1],[0,1]],[[1,0],[1,1]])
        packet=switching_enrichment(modes,[[1,0]],[0,0],[1,1],[[1]]*80,0,max_rank=2,max_symbolic_states=1,seed_rows=[[0,1]])
        self.assertEqual(packet['certificate']['engine'],'exact_common_closure')
        self.assertEqual(packet['certificate']['endpoint_words_at_horizon'],str(2**80))
        self.assertTrue(replay_switching(packet))

    def test_resource_limits_never_accept_partial_bound(self):
        for modes in (specs([[1,1],[0,1]],[[1,2],[0,1]]),specs([[1,1],[0,1]],[[1,0],[1,1]])):
            with self.assertRaises(OverflowError):switching_enrichment(modes,[[1,0]],[0,0],[1,1],[[1]]*10,0,max_rank=1,max_symbolic_states=10)

    def test_tampered_witness_bound_mode_and_projection_rejected(self):
        modes=specs([[1,1],[0,1]],[[1,2],[0,1]])
        original=switching_enrichment(modes,[[1,0]],[0,0],[1,1],[[1]]*3,0,max_rank=1)
        mutations=[lambda p:p['certificate'].__setitem__('worst_bound','0'),lambda p:p['certificate']['attaining_witness']['modes'].__setitem__(0,'forged'),lambda p:p['G']['a'][0].__setitem__(0,'99'),lambda p:p.__setitem__('status','tolerance_met'),lambda p:p['initial_radius'].__setitem__(1,'2')]
        for mutate in mutations:
            p=copy.deepcopy(original);mutate(p)
            with self.assertRaises(ValueError):replay_switching(p)

    def test_multioutput_bounds_and_zero_horizon(self):
        modes=specs([[1,1],[0,1]],[[1,2],[0,1]])
        p=switching_enrichment(modes,[[1,0],[0,1]],[0,0],[1,1],[],0,max_rank=2)
        self.assertEqual(p['certificate']['output_bounds'],[['0','0']]);self.assertTrue(replay_switching(p))

    def test_reduced_noncommutativity_forces_enumeration(self):
        # Diagonal full operators commute, but two-dimensional projection can
        # destroy commutation. The algorithm must check both levels.
        modes=specs([[1,0,0],[0,2,0],[0,0,3]],[[4,0,0],[0,1,0],[0,0,2]],((1,),(0,),(1,)),((0,),(1,),(0,)))
        p=switching_enrichment(modes,[[1,1,0],[0,1,1]],[0]*3,[1]*3,[[1]]*2,0,max_rank=2)
        self.assertEqual(p['certificate']['engine'],'endpoint_word_enumeration');self.assertTrue(replay_switching(p))
