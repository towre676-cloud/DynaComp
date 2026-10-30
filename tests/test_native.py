import copy
import ctypes
from contextlib import contextmanager
import sys
from fractions import Fraction as Q
import random
import os
import shutil
import unittest

from dynacomp.binary64 import bits,from_bits,exact,round_bits
from dynacomp.native import (certify_sequential,replay_sequential,compose_sequential,
                             replay_sequential_composition,PreparedSequential,sequential_budget,
                             _compiler, _compile_command)
from dynacomp.rational import q
from dynacomp.uniform import SparseLinearModel
from dynacomp.validated import validated_box

try:
    import numpy as np
    AVAILABLE = bool(os.environ.get('DYNACOMP_CC') or shutil.which('cc') or shutil.which('gcc') or shutil.which('clang') or shutil.which('x86_64-w64-mingw32-gcc'))
except ImportError:
    AVAILABLE = False


def reference_dot(coefficients,words):
    total=bits(0.0)
    for a,b in zip(coefficients,words):
        product=round_bits(exact(a)*exact(b))
        total=round_bits(exact(total)+exact(product))
    return total


@contextmanager
def environment_control():
    """Compile fenv constants on the target, without assuming Linux values."""
    import subprocess
    import tempfile
    from pathlib import Path
    with tempfile.TemporaryDirectory() as folder:
        source=Path(folder)/'control.c'
        target=Path(folder)/('control.dll' if sys.platform=='win32' else 'control.so')
        source.write_text("""
#include <fenv.h>
#if defined(__x86_64__) || defined(_M_X64)
#include <xmmintrin.h>
#else
static unsigned _mm_getcsr(void) { return 0; }
static void _mm_setcsr(unsigned v) { (void)v; }
#endif
#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif
EXPORT int get_round(void) { return fegetround(); }
EXPORT int set_round(int v) { return fesetround(v); }
EXPORT int set_downward(void) { return fesetround(FE_DOWNWARD); }
EXPORT unsigned get_control(void) { return _mm_getcsr(); }
EXPORT void set_control(unsigned v) { _mm_setcsr(v); }
""")
        subprocess.run(_compile_command(_compiler(),source,target),check=True,capture_output=True)
        control=ctypes.CDLL(str(target))
        control.get_control.restype=ctypes.c_uint
        control.set_control.argtypes=[ctypes.c_uint]
        control.set_round.argtypes=[ctypes.c_int]
        try:yield control
        finally:
            if sys.platform=='win32':
                from _ctypes import FreeLibrary
                FreeLibrary(control._handle)


@unittest.skipUnless(AVAILABLE,'NumPy and a C compiler required')
class NativeTests(unittest.TestCase):
    def test_random_traces_match_independent_integer_rounding(self):
        rng=random.Random(1011)
        for _ in range(16):
            e=[[Q(rng.randrange(-6,7),8) for j in range(3)] for i in range(2)]
            g=[[Q(rng.randrange(-2,3),8) for j in range(2)] for i in range(2)]
            h=[[Q(rng.randrange(-3,4),8)] for i in range(2)]
            d=[[Q(rng.randrange(-5,6),7) for j in range(2)]]
            packet=certify_sequential(e,g,h,d,[0]*3,[1]*3,[[1]]*5,'1/1000000')
            runner=PreparedSequential(packet)
            initial=np.asarray([[rng.uniform(-1,1) for j in range(4)] for i in range(3)])
            inputs=np.asarray([[[rng.uniform(-1,1) for j in range(4)]] for t in range(5)])
            result=runner.simulate(initial,inputs)
            words=packet['loaded_words']
            for j in range(4):
                x=[bits(float(v)) for v in initial[:,j]]
                z=[reference_dot(row,x) for row in words['E']]
                rational=[sum(a*exact(b) for a,b in zip(row,x)) for row in e]
                for t in range(6):
                    expected=reference_dot(words['D'][0],z)
                    self.assertEqual(bits(float(result[t,0,j])),expected)
                    ideal=sum(a*b for a,b in zip(d[0],rational))
                    self.assertLessEqual(abs(exact(expected)-ideal),q(packet['output_error_upper'][t][0]))
                    if t<5:
                        u=[bits(float(inputs[t,0,j]))]
                        z=[reference_dot(a+b,z+u) for a,b in zip(words['G'],words['H'])]
                        rational=[sum(a*b for a,b in zip(row,rational))+forcing[0]*exact(u[0]) for row,forcing in zip(g,h)]

    def test_sequential_cancellation_is_not_exact_dot(self):
        packet=certify_sequential([[1,1,1]],[[0]],[[0]],[[1]],[0]*3,[10**16]*3,[],'100')
        runner=PreparedSequential(packet)
        result=runner.simulate([1e16,1.0,-1e16],np.empty((0,1,1)))
        self.assertEqual(result[0,0,0],0.0)
        self.assertLessEqual(Q(1),q(packet['worst_execution_upper']))

    def test_subnormal_transport_and_copy_guards(self):
        tiny=Q(1,2**1074)
        packet=certify_sequential([[1]],[[1]],[[0]],[[1]],[tiny],[0],[[0]],'1/100')
        runner=PreparedSequential(packet)
        output=runner.simulate([float(tiny)],np.zeros((1,1,1)))
        self.assertEqual(bits(float(output[0,0,0])),1)
        self.assertEqual(bits(float(output[1,0,0])),1)
        packet['loaded_words']['E'][0][0]=0
        public=runner.packet;public['loaded_words']['E'][0][0]=0
        self.assertEqual(bits(float(runner.simulate([float(tiny)],np.zeros((1,1,1)))[0,0,0])),1)
        with self.assertRaises(ValueError):runner.simulate([0.0],np.zeros((1,1,1)))
        with self.assertRaises(ValueError):runner.simulate([float(tiny)],np.zeros((2,1,1)))
        with self.assertRaises(ValueError):runner.simulate([float('nan')],np.zeros((1,1,1)))

    def test_environment_refuses_changed_rounding(self):
        packet=certify_sequential([[1]],[[0]],[[0]],[[1]],[0],[1],[],'1/100')
        runner=PreparedSequential(packet)
        with environment_control() as control:
            old=control.get_round()
            try:
                self.assertEqual(control.set_downward(),0)
                with self.assertRaises(RuntimeError):runner.simulate([0.25],np.empty((0,1,1)))
            finally:control.set_round(old)
        self.assertEqual(runner.simulate([0.25],np.empty((0,1,1)))[0,0,0],0.25)

    def test_replay_and_composition_tampering(self):
        model=SparseLinearModel.from_spec({'n':1,'entries':[[0,0,'1/2']],'B':[[1]],'C':[[1]]})
        source=validated_box(model,[[1]],[[Q(1,2)]],[[1]],[[1]],[0],[1],[[1]]*3,'1/100')
        composed=compose_sequential(source,'1/100')
        self.assertTrue(replay_sequential_composition(composed))
        self.assertTrue(replay_sequential(composed['execution']))
        bad=copy.deepcopy(composed);bad['combined_upper']='0'
        with self.assertRaises(ValueError):replay_sequential_composition(bad)
        bad=copy.deepcopy(composed['execution']);bad['contract']='BLAS'
        with self.assertRaises(ValueError):replay_sequential(bad)

    def test_flush_and_denormals_modes_are_refused(self):
        import platform
        if platform.machine().lower() not in ('x86_64','amd64'):
            self.skipTest('x86 control-register test')
        with environment_control() as control:
            packet=certify_sequential([[1]],[[0]],[[0]],[[1]],[0],[1],[],'1/100')
            runner=PreparedSequential(packet);old=control.get_control()
            for mask in (0x40,0x8000,0x8040):
                try:
                    control.set_control(old|mask)
                    with self.assertRaises(RuntimeError):runner.simulate([0.25],np.empty((0,1,1)))
                finally:control.set_control(old)

    def test_halfway_products_additions_and_subnormal_generation(self):
        # Exact integer rounding is independent of the C compiler and host fenv.
        # Half-ulp sums, half-minsubnormal ties, and the normal boundary.
        tiny=Q(1,2**1074)
        cases=[([1,1],[1,Q(1,2**53)]),
               ([1,1],[1+Q(1,2**52),Q(1,2**53)]),
               ([Q(1,2)],[tiny]), ([Q(3,2)],[tiny]),
               ([Q(1,2)],[Q(1,2**1022)]),
               ([1+Q(1,2**52)],[1+Q(1,2**52)])]
        for coefficients, values in cases:
            packet=certify_sequential([coefficients],[[0]],[[0]],[[1]],
                                      values,[0]*len(values),[],'1/100')
            runner=PreparedSequential(packet)
            expected=reference_dot(packet['loaded_words']['E'][0],
                                   [round_bits(v) for v in values])
            result=runner.simulate([float(v) for v in values],np.empty((0,1,1)))
            self.assertEqual(bits(float(result[0,0,0])),expected)
            ideal=sum(a*b for a,b in zip(coefficients,values))
            self.assertLessEqual(abs(exact(expected)-ideal),q(packet['output_error_upper'][0][0]))
            self.assertTrue(np.array_equal(result,runner.simulate([float(v) for v in values],np.empty((0,1,1)))))

    def test_overflow_fails_closed(self):
        with self.assertRaises(ValueError):
            certify_sequential([[2]],[[0]],[[0]],[[1]],[0],[Q.from_float(float.fromhex('0x1.fffffffffffffp+1023'))],[],'1')
        self.assertEqual(sequential_budget(Q(0),3),0)
        with self.assertRaises(ValueError):sequential_budget(Q(1),2**53)


if __name__=='__main__':unittest.main()
