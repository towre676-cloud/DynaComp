import importlib.util
import sys
import unittest
from pathlib import Path
from dynacomp.network import compile_network,run_scenarios
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"applications"/"network_scenarios"))
from make_demo import model,scenarios


@unittest.skipUnless(importlib.util.find_spec("numpy") and importlib.util.find_spec("scipy"),"optional numerical dependencies")
class NumericNetworkTests(unittest.TestCase):
    def test_production_batch_backend_matches_standard_mixed_horizons(self):
        plan=compile_network(model(3))
        batch=scenarios(3,10)
        batch["scenarios"][1]["replenishment"]=[]
        batch["scenarios"][2]["replenishment"]=batch["scenarios"][2]["replenishment"][:3]
        standard=run_scenarios(plan,batch)
        for engine in ("full","reduced"):
            numeric=run_scenarios(plan,batch,engine,backend="numpy")
            self.assertEqual([s["name"] for s in numeric["scenarios"]],[s["name"] for s in standard["scenarios"]])
            for a,b in zip(standard["scenarios"],numeric["scenarios"]):
                for row,other in zip(a["stock_trace"]+a["dispatch_trace"],b["stock_trace"]+b["dispatch_trace"]):
                    for x,y in zip(row,other):
                        self.assertAlmostEqual(x,y,places=10)

    def test_sparse_no_compression_and_zero_transition(self):
        for source in (model(8,heterogeneous=True),{"schema":"dynacomp.network-model/1","n":2,"region_names":["one"],"node_regions":[0,0],"transitions":[],"inputs":[{"name":"supply","distribution":[[0,1]]}]}):
            plan=compile_network(source,max_groups=32)
            batch={"schema":"dynacomp.network-scenarios/1","scenarios":[{"name":"one","initial_stock":[1]*source["n"],"replenishment":[[1]*len(source["inputs"])]*2}]}
            full=run_scenarios(plan,batch,"full",backend="numpy")
            reduced=run_scenarios(plan,batch,backend="numpy")
            self.assertEqual(full["scenarios"],reduced["scenarios"])

    def test_numeric_dimensions_and_negative_values_rejected(self):
        import numpy as np
        from dynacomp.network import Network
        from dynacomp.network_numeric import NumpyEvolution
        source=model(2)
        runner=NumpyEvolution(Network.load(source),compile_network(source),"reduced")
        for initial,inputs in ((np.ones((7,1)),np.zeros((2,4,1))),(-np.ones((8,1)),np.zeros((2,4,1))),(np.ones((8,1)),np.zeros((2,3,1))),(np.full((8,1),np.nan),np.zeros((2,4,1)))):
            with self.assertRaises(ValueError):
                runner.simulate(initial,inputs)

    def test_exact_backend_requirement(self):
        with self.assertRaises(ValueError):
            run_scenarios(compile_network(model(2)),scenarios(2,1),exact=True,backend="numpy")
