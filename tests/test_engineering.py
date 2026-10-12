import copy
import tempfile
import unittest
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
from scipy.io import savemat
from scipy.linalg import expm
from dynacomp.pilot import import_mat,workload,compile_workload,simulate_plan,evaluate_plan
from dynacomp.engineering import replay_import,enclosure,output_contract
from dynacomp.rational import matrix


class EngineeringTests(unittest.TestCase):
    def fixture(self,folder,**changes):
        data={'A':np.array([[-1.,.25],[0.,-2.]]),'B':np.array([[1.],[2.]]),'C':np.array([[2.,-1.]]),
              'E':np.array([[2.,.5],[0.,3.]]),'D':np.array([[.75]])}
        data.update(changes);path=Path(folder)/'source.mat';savemat(path,data);return path,data

    def test_exact_descriptor_discrete_and_unit_preservation(self):
        with tempfile.TemporaryDirectory() as folder:
            path,data=self.fixture(folder)
            spec=import_mat(path);self.assertEqual(replay_import(spec)['status'],'exact')
            a=np.array(spec['model']['A'],dtype=object)
            a=np.array([[float(Q(v)) for v in row] for row in a])
            np.testing.assert_allclose(data['E']@a,data['A'],atol=1e-15)
            self.assertEqual(spec['model']['C'],[['2','-1']]);self.assertEqual(spec['feedthrough'],[['3/4']])
            normalized=import_mat(path,normalize_output=True)
            self.assertEqual(normalized['feedthrough'],[['1/4']])
            self.assertEqual(normalized['import_certificate']['output_divisor'],'3')

    def test_zoh_enclosure_vs_independent_high_precision(self):
        try: import mpmath as mp
        except ImportError: self.skipTest('mpmath numerical test extra required')
        with tempfile.TemporaryDirectory() as folder:
            path,data=self.fixture(folder)
            spec=import_mat(path,'1/7',require_bounded=True);cert=replay_import(spec)
            with mp.workdps(80):
                e=mp.matrix(data['E'].tolist());a=e**-1*mp.matrix(data['A'].tolist());b=e**-1*mp.matrix(data['B'].tolist())
                z=mp.matrix(3);z[:2,:2]=a/7;z[:2,2:3]=b/7;reference=mp.expm(z)
                for field,start,width in [('A',0,2),('B',2,1)]:
                    error=max(sum(abs(reference[i,start+j]-mp.mpf(Q(spec['model'][field][i][j]).numerator)/Q(spec['model'][field][i][j]).denominator) for j in range(width)) for i in range(2))
                    bound=Q(cert[field+'_error_inf']);self.assertLessEqual(error,mp.mpf(bound.numerator)/bound.denominator)
            # Adversarial stored coefficient discrepancies are bounded too.
            bad=copy.deepcopy(spec);bad['model']['A'][0][0]='10'
            with self.assertRaises(ValueError):workload(bad)
            bad['import_certificate']['transformation']=enclosure(matrix([[str(Q(v)) for v in row] for row in np.linalg.solve(data['E'],data['A'])]),matrix([[str(Q(v)) for v in row] for row in np.linalg.solve(data['E'],data['B'])]),Q(1,7),matrix(bad['model']['A']),matrix(bad['model']['B']))
            # The exact mass identities, not a numerical solve, govern replay.
            with self.assertRaises(ValueError):workload(bad)

    def test_refuse_malformed_singular_and_nonfinite(self):
        with tempfile.TemporaryDirectory() as folder:
            for change in ({'E':np.array([[1.,2.],[2.,4.]])},{'D':np.ones((2,1))},{'E':np.array([[np.nan,0],[0,1]])},{'D':np.array([[1j]])},{'B':np.ones((3,1))}):
                path,_=self.fixture(folder,**change)
                with self.assertRaises(ValueError):import_mat(path,'1/10')
            path,_=self.fixture(folder,A=np.eye(2)*1000)
            spec=import_mat(path,'1/10');self.assertEqual(replay_import(spec)['status'],'uncertified')
            with self.assertRaises(ValueError):import_mat(path,'1/10',require_bounded=True)

    def test_feedthrough_terminal_semantics_and_certificate_guards(self):
        with tempfile.TemporaryDirectory() as folder:
            path,_=self.fixture(folder,E=np.eye(2),A=np.diag([-1.,-2.]),C=np.array([[1.,0.]]),D=np.array([[2.]]),B=np.array([[1.],[0.]]))
            spec=import_mat(path,'1/10',require_bounded=True);spec['input_radii']=[[1]]*3
            plan=compile_workload(spec,Path(folder)/'plan');self.assertTrue(plan['methods']['observable']['accepted'])
            u=np.array([[[.25]],[[.5]],[[-.75]]]);y=simulate_plan(Path(folder)/'plan','observable',np.zeros(2),u)
            state=0.;reference=[]
            for t in range(4):
                reference.append(state+2*u[min(t,2),0,0])
                if t<3:state=np.exp(-.1)*state+(1-np.exp(-.1))*u[t,0,0]
            np.testing.assert_allclose(y[:,0,0],reference,atol=1e-14)
            result=evaluate_plan(Path(folder)/'plan',8,3,10)
            self.assertEqual(result['engineering_reference'],'continuous ZOH samples')
            self.assertTrue(result['continuous_sampling_bounded'])
            bad=copy.deepcopy(spec);bad['feedthrough'][0][0]='3'
            with self.assertRaises(ValueError):workload(bad)
            bad=copy.deepcopy(spec);bad['output_input_semantics']='left'
            with self.assertRaises(ValueError):workload(bad)

    def test_feedthrough_overflow_and_zero_horizon_refused(self):
        spec={'schema':'dynacomp.workload/1','model':{'A':[[0]],'B':[[0]],'C':[[1]]},
              'feedthrough':[['17976931348623157'+'0'*292]],'initial_center':[0],'initial_radius':[1],
              'input_radii':[[1]],'tolerance':1,'min_rank':1,'max_rank':1}
        model,domain,*_=workload(spec)
        packet={'combined_output_upper':[['0'],['0']]}
        with self.assertRaises(ValueError):output_contract(spec,model,domain,packet)
        spec['input_radii']=[]
        with self.assertRaises(ValueError):workload(spec)

    def test_large_dominant_mass_residual_replay(self):
        from dynacomp.engineering import descriptor_transform
        n=33
        e=matrix([[3 if i==j else -1 if abs(i-j)==1 else 0 for j in range(n)] for i in range(n)])
        a=matrix([[int(i==j) for j in range(n)] for i in range(n)]);b=matrix([[1] for _ in range(n)])
        sa,sb,cert=descriptor_transform(a,b,e)
        self.assertEqual(cert['status'],'bounded')
        self.assertEqual(cert['mass_dominance_margin'],'1')
        self.assertEqual(descriptor_transform(a,b,e,(sa,sb))[2],cert)
        bad=matrix([[v+int(i==0 and j==0) for j,v in enumerate(row)] for i,row in enumerate(sa)])
        self.assertGreater(Q(descriptor_transform(a,b,e,(bad,sb))[2]['A_error_inf']),Q(cert['A_error_inf']))
        with self.assertRaises(ValueError):descriptor_transform(a,b,matrix([[int(i==j)*0 for j in range(n)] for i in range(n)]))

    def test_native_feedthrough_matches_independently_rounded_operations(self):
        import ctypes
        from dynacomp.native import _library
        from dynacomp.binary64 import round_bits,exact,bits,from_bits
        lib=_library();rng=np.random.default_rng(17)
        f=np.ascontiguousarray(rng.normal(size=(3,4)));u=np.ascontiguousarray(rng.normal(size=(2,4,7)))
        base=np.ascontiguousarray(rng.normal(size=(3,3,7)));out=base.copy();pointer=ctypes.POINTER(ctypes.c_double)
        status=lib.feedthrough(4,3,2,7,*(v.ctypes.data_as(pointer) for v in (f,u,out)));self.assertEqual(status,0)
        for t in range(3):
            for i in range(3):
                for j in range(7):
                    value=Q(0)
                    for k in range(4):
                        product=exact(round_bits(Q(float(f[i,k]))*Q(float(u[min(t,1),k,j]))))
                        value=exact(round_bits(value+product))
                    expected=round_bits(Q(float(base[t,i,j]))+value)
                    self.assertEqual(bits(float(out[t,i,j])),expected)

    def test_horizon_propagation_and_uncertified_scope(self):
        with tempfile.TemporaryDirectory() as folder:
            path,_=self.fixture(folder)
            spec=import_mat(path,'1/10');model,domain,*_=workload(spec)
            bound=output_contract(spec,model,domain)
            self.assertEqual(bound['import_output_upper'][0],['0'])
            self.assertGreater(Q(bound['import_output_upper'][-1][0]),0)
            altered=copy.deepcopy(spec);altered['input_radii']=[[10]]*20
            model,domain,*_=workload(altered)
            self.assertGreater(Q(output_contract(altered,model,domain)['import_output_upper'][-1][0]),Q(bound['import_output_upper'][-1][0]))

if __name__=='__main__':unittest.main()
