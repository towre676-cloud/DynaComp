import importlib.util
import unittest
from dynacomp.thermal import plate_model
from dynacomp import uniform_enrichment


@unittest.skipUnless(importlib.util.find_spec('numpy') and importlib.util.find_spec('scipy'),'optional numerical dependencies')
class ThermalNumericTests(unittest.TestCase):
    def test_batch_difference_respects_bound_with_roundoff_margin(self):
        import numpy as np
        from fractions import Fraction
        from dynacomp.thermal_numeric import run_thermal_batch
        model=plate_model(6)
        packet=uniform_enrichment(model,[20]*36,[2]*36,[[1]]*8,'1/10',input_centers=[[1]]*8,strategy='prefix')
        rng=np.random.default_rng(510)
        initial=rng.uniform(18,22,(36,12));inputs=rng.uniform(0,2,(8,1,12))
        full=run_thermal_batch(packet,initial,inputs,'full')
        reduced=run_thermal_batch(packet,initial,inputs)
        bounds=np.asarray([float(Fraction(v)) for v in packet['certificate']['bound']])
        self.assertTrue((np.max(np.abs(full-reduced),axis=(1,2))<=bounds+1e-10).all())

    def test_out_of_domain_rejected(self):
        import numpy as np
        from dynacomp.thermal_numeric import run_thermal_batch
        packet=uniform_enrichment(plate_model(4),[20]*16,[2]*16,[[1]]*2,0,input_centers=[[1]]*2)
        for initial,inputs in ((np.full((16,1),23),np.ones((2,1,1))),(np.full((16,1),20),np.full((2,1,1),3)),(np.full((16,1),20),np.ones((3,1,1)))):
            with self.assertRaises(ValueError):run_thermal_batch(packet,initial,inputs)

    def test_binary64_rounded_endpoint_does_not_admit_outside_value(self):
        import numpy as np
        from dynacomp import LinearModel
        from dynacomp.thermal_numeric import run_thermal_batch
        packet=uniform_enrichment(LinearModel.from_rows([[1]],[[1]],[[1]]),[0],['1/10'],[],0)
        # binary64 0.1 is strictly larger than the rational endpoint 1/10.
        with self.assertRaises(ValueError):run_thermal_batch(packet,[[0.1]],np.empty((0,1,1)))
