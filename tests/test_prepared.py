import importlib.util
import unittest
from fractions import Fraction as Q
from dynacomp.thermal import plate_model
from dynacomp.robust import switching_enrichment


@unittest.skipUnless(importlib.util.find_spec('numpy') and importlib.util.find_spec('scipy'),'optional numerical dependencies')
class PreparedTests(unittest.TestCase):
    def packet(self):
        modes={s:plate_model(6,k).spec() for s,k in [('low','1/16'),('high','1/8')]}
        return switching_enrichment(modes,modes['low']['C'],[20]*36,[2]*36,[[1]]*6,'1/10',input_centers=[[1]]*6,max_rank=6)

    def test_reuse_shared_initial_and_interpolated_controls(self):
        import numpy as np
        from dynacomp.prepared import PreparedSwitching
        packet=self.packet();reduced=PreparedSwitching(packet);full=PreparedSwitching(packet,'full')
        rng=np.random.default_rng(610)
        initial=np.full(36,20.);inputs=rng.uniform(0,2,(6,1,12));mixing=rng.uniform(0,1,(6,12))
        bounds=np.asarray([float(Q(v)) for v in packet['certificate']['bound']])
        for _ in range(3):
            r=reduced.simulate(initial,inputs,mixing);f=full.simulate(initial,inputs,mixing)
            self.assertTrue((np.max(np.abs(f-r),axis=(1,2))<=bounds+1e-10).all())
            expanded=reduced.simulate(np.repeat(initial[:,None],12,axis=1),inputs,mixing)
            self.assertTrue(np.allclose(r,expanded,rtol=1e-12,atol=1e-12))
        self.assertEqual(reduced.cache_stats,{'packet_replays':1,'batches':6,'initial_fields_encoded':39})

    def test_source_mutation_does_not_change_prepared_plan(self):
        import numpy as np
        from dynacomp.prepared import PreparedSwitching
        packet=self.packet();runner=PreparedSwitching(packet)
        packet['initial_radius'][0]='0';packet['G']['low'][0][0]='99'
        runner.simulate(np.full(36,21),np.ones((6,1,1)),np.zeros(6))
        self.assertEqual(runner.cache_stats['packet_replays'],1)

    def test_invalid_domain_modes_dimensions_rejected(self):
        import numpy as np
        from dynacomp.prepared import PreparedSwitching
        runner=PreparedSwitching(self.packet())
        for initial,inputs,mixing in ((np.full(36,23),np.ones((6,1,1)),np.zeros(6)),(np.full(36,20),np.full((6,1,1),3),np.zeros(6)),(np.full(36,20),np.ones((6,1,1)),np.full(6,1.1)),(np.full(36,20),np.ones((5,1,1)),np.zeros(5))):
            with self.assertRaises(ValueError):runner.simulate(initial,inputs,mixing)

    def test_false_packet_rejected_before_preparation(self):
        from dynacomp.prepared import PreparedSwitching
        packet=self.packet();packet['certificate']['worst_bound']='99'
        with self.assertRaises(ValueError):PreparedSwitching(packet)

    def test_directed_binary64_membership_matches_rational_comparison(self):
        import numpy as np
        from dynacomp.prepared import PreparedSwitching
        for lo,hi in ((Q(-1,10),Q(1,10)),(Q(1,3),Q(2,3)),(Q(0),Q(0)),(Q(10)**400,Q(10)**401),(-Q(10)**401,-Q(10)**400)):
            lower,upper=PreparedSwitching._float_endpoints(lo,hi)
            candidates=[0.,float(np.finfo(float).max),-float(np.finfo(float).max)]
            if abs(lo)<Q(10)**300:
                for v in (float(lo),float(hi)):
                    candidates.extend([v,float(np.nextafter(v,-np.inf)),float(np.nextafter(v,np.inf))])
            for v in candidates:
                self.assertEqual(lower<=v<=upper,lo<=Q(v)<=hi)
