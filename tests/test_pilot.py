import copy
import gzip
import json
from pathlib import Path
import os
import shutil
import tempfile
import unittest
from dynacomp.pilot import workload,cost_comparison,render_report
try:
    import scipy
    import numpy as np
    import threadpoolctl
    AVAILABLE=bool(os.environ.get('DYNACOMP_CC') or shutil.which('cc') or shutil.which('gcc') or shutil.which('clang') or shutil.which('x86_64-w64-mingw32-gcc'))
except ImportError:AVAILABLE=False

SPEC={'schema':'dynacomp.workload/1','name':'<sample>','model':{'A':[['1/2',0],[0,'1/4']],'B':[[1],[1]],'C':[[1,'1/1000']]},
      'initial_center':[0,0],'initial_radius':[1,1],'input_radii':[[1]]*4,'tolerance':'1/100','min_rank':1,'max_rank':2}


class PilotTests(unittest.TestCase):
    def test_costs_include_setup_and_no_false_break_even(self):
        result=cost_comparison(12,2,3,1,5)
        self.assertEqual(result['break_even_batches'],5)
        self.assertEqual(result['total_speedup'],1)
        result=cost_comparison(12,2,1,3,5)
        self.assertIsNone(result['break_even_batches'])
        self.assertLess(result['total_speedup'],1)
        with self.assertRaises(ValueError):cost_comparison(1,1,1,1,0)

    def test_invalid_workloads(self):
        for change in ({'tolerance':0},{'min_rank':0},{'max_rank':3},{'initial_radius':[-1,1]}, {'input_radii':[[1,2]]}):
            bad=copy.deepcopy(SPEC);bad.update(change)
            with self.assertRaises(ValueError):workload(bad)
        bad=copy.deepcopy(SPEC);bad['model']['A'][0][0]=0.5
        with self.assertRaises(TypeError):workload(bad)

    @unittest.skipUnless(AVAILABLE,'Numerical extras and compiler required')
    def test_saved_plan_cli_workflow_and_contract_guards(self):
        from dynacomp.pilot import compile_workload,evaluate_plan,simulate_plan
        with tempfile.TemporaryDirectory() as folder:
            plan=compile_workload(SPEC,folder)
            self.assertTrue(all(result['accepted'] for result in plan['methods'].values()))
            result=evaluate_plan(folder,queries=8,trials=3,batches=10)
            self.assertEqual(result['states'],2)
            self.assertEqual(result['steps'],4)
            self.assertEqual(len(result['full_blas_raw_batch_seconds']),3)
            self.assertIn('&lt;sample&gt;',render_report(result))
            plan_file=Path(folder)/'plan.json'
            (Path(folder)/'plan.json.gz').write_bytes(gzip.compress(plan_file.read_bytes(),mtime=0))
            plan_file.unlink()
            supplied=simulate_plan(folder,'observable',np.zeros(2),np.zeros((4,1,8)))
            self.assertEqual(supplied.shape,(5,1,8))
            self.assertTrue(np.array_equal(supplied,np.zeros_like(supplied)))
            with self.assertRaises(ValueError):simulate_plan(folder,'observable',np.zeros(2),np.ones((4,1,8))*2)
            for method in result['methods'].values():
                self.assertEqual(method['rank'],1)
                self.assertLess(float(__import__('fractions').Fraction(method['combined_upper'])),0.01)
            artifact=Path(folder)/plan['methods']['observable']['packet']
            artifact.write_bytes(artifact.read_bytes()+b'changed')
            with self.assertRaises(ValueError):evaluate_plan(folder,8,3,10)

    @unittest.skipUnless(AVAILABLE,'Numerical extras and compiler required')
    def test_mat_sampling_and_refusal(self):
        from scipy.io import savemat
        from dynacomp.pilot import import_mat
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'input.mat'
            data={'A':np.asarray([[-1.0]]),'B':np.asarray([[1.0]]),'C':np.asarray([[2.0]])}
            savemat(path,data)
            spec=import_mat(path,'1/10',True)
            model,_,_,_,_=workload(spec)
            self.assertEqual(model.C,((1,),))
            self.assertAlmostEqual(float(model.operator.entries[0][2]),np.exp(-0.1))
            self.assertIn('sampling/physical',spec['provenance']['scope'])
            for extra in ({'E':np.asarray([[2.0]])},{'D':np.asarray([[1.0]])}):
                savemat(path,{**data,**extra})
                with self.assertRaises(ValueError):import_mat(path,'1/10')


if __name__=='__main__':unittest.main()
