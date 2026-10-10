import copy
import random
import unittest
from fractions import Fraction as Q
from dynacomp import compile_sparse_family,replay_sparse_family
from dynacomp.switched import compile_switched
from dynacomp.rational import matrix,row_basis,mul,mv,addv


class SparseFamilyTests(unittest.TestCase):
    def test_matches_dense_switched_random_families(self):
        rng=random.Random(450)
        for _ in range(12):
            matrices={s:[[Q(rng.randrange(-1,2),3) for j in range(3)] for i in range(3)] for s in ('a','b')}
            inputs={s:[[1],[0],[1]] for s in matrices}
            specs={s:{'n':3,'entries':[[i,j,str(v)] for i,row in enumerate(a) for j,v in enumerate(row) if v],'B':inputs[s]} for s,a in matrices.items()}
            sparse=compile_sparse_family(specs,[[1,0,0]],3)
            dense=compile_switched(matrices,[[1,0,0]],inputs)
            self.assertEqual(row_basis(matrix(sparse['E'])),row_basis(matrix(dense['E'])))
            self.assertTrue(replay_sparse_family(sparse))

    def test_arbitrary_interpolation_and_switching_inputs(self):
        modes={s:{'n':3,'entries':entries,'B':[[1],[0],[1]]} for s,entries in {'a':[[0,0,1],[1,1,1],[2,2,1]],'b':[[0,0,1],[0,1,1],[1,1,1],[2,2,2]]}.items()}
        p=compile_sparse_family(modes,[[1,0,0]],3)
        e,d=matrix(p['E']),matrix(p['D'])
        from dynacomp.uniform import SparseLinearModel
        models={s:SparseLinearModel.from_spec({**spec,'C':[[1,0,0]]}) for s,spec in modes.items()}
        x=(Q(1),Q(2),Q(3));z=mv(e,x)
        for lam in [Q(1,3),Q(2,5),Q(-2),Q(1)]:
            u=(Q(2),)
            a,b=(addv(models[s].operator.step(x),mv(models[s].B,u)) for s in ('a','b'))
            x=tuple((1-lam)*v+lam*w for v,w in zip(a,b))
            a,b=(addv(mv(matrix(p['G'][s]),z),mv(matrix(p['H'][s]),u)) for s in ('a','b'))
            z=tuple((1-lam)*v+lam*w for v,w in zip(a,b))
            self.assertEqual(mv(d,z),(x[0],))

    def test_budget_and_tampering(self):
        modes={'a':{'n':2,'entries':[[0,1,1],[1,0,1]],'B':[[1],[0]]}}
        with self.assertRaises(OverflowError):compile_sparse_family(modes,[[1,0]],1)
        p=compile_sparse_family(modes,[[1,0]],2)
        for key in ('G','H'):
            bad=copy.deepcopy(p);bad[key]['a'][0][0]='99'
            with self.assertRaises(ValueError):replay_sparse_family(bad)

    def test_mismatched_input_width_rejected(self):
        modes={'a':{'n':1,'entries':[[0,0,1]],'B':[[1]]},'b':{'n':1,'entries':[[0,0,1]],'B':[[1,1]]}}
        with self.assertRaises(ValueError):compile_sparse_family(modes,[[1]],1)
