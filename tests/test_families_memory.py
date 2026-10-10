import copy
import random
import unittest
from fractions import Fraction as Q
from dynacomp import LinearModel, compile_linear, compile_family, replay_family, closure_conditions
from dynacomp.families import specialize
from dynacomp.memory import compile_output_memory, replay_memory, simulate_memory
from dynacomp.rational import matrix, mv, addv


class FamilyTests(unittest.TestCase):
    def models(self):
        return [LinearModel.from_rows([["1/2",0,0],[0,"1/3",0],[0,0,2]],[[1],[0],[0]],[[1,0,0]]),
                LinearModel.from_rows([[0,1,0],[0,0,0],[0,0,0]],[[0],[1],[0]],[[1,0,0]])]

    def test_specialization_accident_vs_family(self):
        models=self.models()
        self.assertEqual(len(compile_linear(models[0]).E),1)
        p=compile_family(models)
        self.assertEqual(len(p['E']),2)
        self.assertTrue(replay_family(p))
        conditions=closure_conditions([m.A for m in models],[[1,0,0]])
        self.assertEqual(conditions['equations'],[['0','1']])

    def test_arbitrary_switching_and_coefficients(self):
        models=self.models()
        p=compile_family(models)
        x=(Q(2),Q(3),Q(4)); z=mv(matrix(p['E']),x)
        for parameters,u in (([1,0],[1]),([2,3],[-1]),([-2,5],[2]),([0,1],[3])):
            def combine(key):
                items=[matrix(v) for v in p[key]]
                return matrix([[sum(Q(k)*item[i][j] for k,item in zip(parameters,items))
                                for j in range(len(items[0][0]))] for i in range(len(items[0]))])
            x=addv(mv(combine('A'),x),mv(combine('B'),u))
            z=addv(mv(combine('G'),z),mv(combine('H'),u))
            self.assertEqual(mv(matrix(p['C']),x),mv(matrix(p['D']),z))

    def test_nonclosing_specialized_encoder_detected(self):
        models=self.models()
        a=specialize([m.A for m in models],[2])
        m=LinearModel.from_rows(a,[[0],[0],[0]],[[1,0,0]])
        self.assertEqual(len(compile_linear(m).E),2)

    def test_mutated_family_rejected(self):
        p=compile_family(self.models())
        p['G'][1][0][0]='3'
        with self.assertRaises(ValueError):
            replay_family(p)


class MemoryTests(unittest.TestCase):
    def test_integrator_input_memory(self):
        m=LinearModel.from_rows([[1,1],[0,1]],[[0],[1]],[[1,0]])
        p=compile_output_memory(m)
        self.assertEqual(p['order'],2)
        self.assertEqual(p['output_coefficients'],['-1','2'])
        self.assertEqual(p['input_kernels'],[['1'],['0']])
        self.assertTrue(replay_memory(p))
        self.assertEqual(simulate_memory(p,[2,5],[[1]]*6),(Q(2),Q(5),Q(9),Q(14),Q(20),Q(27),Q(35)))

    def test_random_scalar_outputs_with_two_inputs(self):
        rng=random.Random(31)
        for n in range(1,5):
            for _ in range(6):
                m=LinearModel.from_rows([[Q(rng.randrange(-2,3),3) for _ in range(n)] for _ in range(n)],
                    [[rng.randrange(-2,3),rng.randrange(-2,3)] for _ in range(n)],[[1]+[0]*(n-1)])
                p=compile_output_memory(m)
                self.assertTrue(replay_memory(p))
                x=tuple(Q(i+1) for i in range(n))
                controls=[[rng.randrange(-2,3),rng.randrange(-2,3)] for _ in range(9)]
                ys=[mv(m.C,x)[0]]
                for u in controls:
                    x=addv(mv(m.A,x),mv(m.B,u)); ys.append(mv(m.C,x)[0])
                self.assertEqual(simulate_memory(p,ys[:p['order']],controls),tuple(ys))
                self.assertEqual(p['order'],len(compile_linear(m).E))

    def test_mutated_input_memory_rejected(self):
        m=LinearModel.from_rows([[1,1],[0,1]],[[0],[1]],[[1,0]])
        p=compile_output_memory(m); p['input_kernels'][0][0]='0'
        with self.assertRaises(ValueError): replay_memory(p)

    def test_non_scalar_target_rejected(self):
        m=LinearModel.from_rows([[1,0],[0,1]],[[0],[0]],[[1,0],[0,1]])
        with self.assertRaises(ValueError): compile_output_memory(m)
