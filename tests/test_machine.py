import copy
import random
import unittest
from fractions import Fraction as Q
from dynacomp.machine import (certify_execution, replay_execution, compose_execution,
                             replay_machine_composition, Binary64Runtime)
from dynacomp.binary64 import bits, exact, INF
from dynacomp.rational import matrix, mv, addv
from dynacomp.uniform import SparseLinearModel
from dynacomp.validated import validated_box


class MachineTests(unittest.TestCase):
    def test_random_traces_enclosed_in_uniform_rounding_budget(self):
        rng=random.Random(805)
        for _ in range(30):
            n,r,m,p=3,2,2,2
            e,g,h,d=[matrix([[Q(rng.randrange(-5,6),7) for j in range(b)] for i in range(a)])
                     for a,b in [(r,n),(r,r),(r,m),(p,r)]]
            packet=certify_execution(e,g,h,d,[0]*n,[1]*n,[[1]*m]*10,1)
            runtime=Binary64Runtime(packet)
            self.assertTrue(replay_execution(packet))
            for scenario in range(6):
                x=[rng.uniform(-1,1) for i in range(n)]
                inputs=[[rng.uniform(-1,1) for i in range(m)] for j in range(10)]
                words=runtime.simulate_bits(list(map(bits,x)),[list(map(bits,row)) for row in inputs])
                state=mv(e,tuple(Q.from_float(v) for v in x))
                for t,row in enumerate(words):
                    ideal=mv(d,state)
                    for actual,reference,bound in zip(row,ideal,packet['output_error_upper'][t]):
                        self.assertLessEqual(abs(exact(actual)-reference),Q(bound))
                    if t<10:state=addv(mv(g,state),mv(h,tuple(Q.from_float(v) for v in inputs[t])))

    def test_guarded_domain_horizon_and_packet_copy(self):
        packet=certify_execution([[1]],[[Q(1,3)]],[[1]],[[1]],[0],[1],[[1]]*3,'1/100000000')
        runtime=Binary64Runtime(packet);packet['loaded_words']['G'][0][0]=INF
        exposed=runtime.packet;exposed['initial_radius']=['1000'];exposed['loaded_words']['E'][0][0]=INF
        out=runtime.simulate([1.],[[0.]]*3)
        self.assertEqual(out[0],(1.,))
        with self.assertRaises(ValueError):runtime.simulate([1.01],[])
        with self.assertRaises(ValueError):runtime.simulate([1.],[[0.]]*4)
        with self.assertRaises(ValueError):runtime.simulate([1.],[[2.]])
        with self.assertRaises(ValueError):runtime.simulate_bits([INF],[])
        with self.assertRaises(TypeError):runtime.simulate([1],[])

    def test_rational_endpoint_guard_no_float_widening(self):
        packet=certify_execution([[1]],[[1]],[[0]],[[1]],['1/10'],[0],[],1)
        runtime=Binary64Runtime(packet)
        with self.assertRaises(ValueError):runtime.simulate([0.1],[])

    def test_subnormal_and_exact_zero(self):
        p=certify_execution([[1]],[[Q(1,2)]],[[0]],[[1]],[0],[Q(1,1 << 1074)],[[0]]*2,1)
        result=Binary64Runtime(p).simulate_bits([1],[[0],[0]])
        self.assertEqual(result,((1,),(0,),(0,)))
        p=certify_execution([[0]],[[0]],[[0]],[[0]],[0],[1],[[0]],0)
        self.assertEqual(p['worst_execution_upper'],'0')
        self.assertEqual(Binary64Runtime(p).simulate([1.],[[0.]]),((0.,),(0.,)))

    def test_overflow_fail_closed_and_nonaccepted_budget(self):
        with self.assertRaises(ValueError):certify_execution([[2]],[[1]],[[0]],[[1]],[0],[Q(1 << 1023)],[],1)
        with self.assertRaises(ValueError):certify_execution([[Q(1 << 1024)]],[[1]],[[0]],[[1]],[0],[0],[],1)
        p=certify_execution([[1]],[[Q(1,3)]],[[0]],[[1]],[0],[1],[[0]],0)
        self.assertEqual(p['status'],'not_certified')
        with self.assertRaises(ValueError):Binary64Runtime(p)
        self.assertTrue(Binary64Runtime(p,require_accepted=False).simulate([1.],[[0.]]))

    def test_composed_source_model_and_machine_error(self):
        model=SparseLinearModel.from_spec({'n':2,'entries':[[0,0,'1/3'],[1,1,'1/2']],
                                          'B':[['1'],['1']],'C':[['1','1/1000']]})
        source=validated_box(model,[[1,0]],[[Q(1,3)]],[[1]],[[1]],[0,0],[1,1],[[1]]*5,'1/100')
        packet=compose_execution(source,'1/100')
        self.assertTrue(replay_machine_composition(packet))
        runtime=Binary64Runtime(packet);x=(Q(1),Q(-1));inputs=[[1.],[-1.],[1.],[-1.],[0.]]
        actual=runtime.simulate([1.,-1.],inputs)
        for t,row in enumerate(actual):
            self.assertLessEqual(abs(Q.from_float(row[0])-mv(model.C,x)[0]),Q(packet['combined_output_upper'][t][0]))
            if t<5:x=addv(model.operator.step(x),mv(model.B,[Q.from_float(inputs[t][0])]))

    def test_all_receipt_parts_replay_and_tampering(self):
        p=certify_execution([[1]],[[Q(1,3)]],[[1]],[[1]],[0],[1],[[1]],1)
        for key in ('loaded_words','output_error_upper','state_error_upper','status','contract'):
            bad=copy.deepcopy(p);bad[key]=None
            with self.assertRaises(ValueError):replay_execution(bad)

    def test_spectral_source_composition_and_tamper(self):
        from dynacomp.spectral import compile_spectral,compose_spectral
        source=compile_spectral({'kind':'power_moments','start':2},3,1)
        e=[[1,0,0]];g=[['1']];h=[[0]];d=[[1]]
        spectral=compose_spectral(source,e,g,h,d,3,1)
        packet=compose_execution(spectral,1)
        self.assertTrue(replay_machine_composition(packet))
        packet['combined_upper']='0'
        with self.assertRaises(ValueError):replay_machine_composition(packet)


if __name__ == '__main__':unittest.main()
