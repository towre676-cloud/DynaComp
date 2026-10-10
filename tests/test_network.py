import copy
import json
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest
from dynacomp.network import compile_network, replay_network, run_scenarios, Network
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"applications"/"network_scenarios"))
from make_demo import model, scenarios


class NetworkTests(unittest.TestCase):
    def test_exact_full_reduced_arbitrary_inputs(self):
        # Compare every output, not just final stock, on 12 independent batches.
        for seed in range(12):
            source = model(3)
            plan = compile_network(source)
            self.assertEqual(plan["rank"], 4)
            batch = scenarios(3, 6, 3, seed)
            full = run_scenarios(plan, batch, "full", exact=True)
            reduced = run_scenarios(plan, batch, exact=True)
            self.assertEqual(full["scenarios"], reduced["scenarios"])
            for item in reduced["scenarios"]:
                self.assertEqual(item["mass_balance_residual"], "0")

    def test_dispatch_timing_and_input_distribution(self):
        source = {"schema":"dynacomp.network-model/1", "n":2, "region_names":["A"], "node_regions":[0,0], "transitions":[[0,0,"1/2"],[1,1,"1/2"]], "inputs":[{"name":"supply","distribution":[[0,"1/3"],[1,"2/3"]]}]}
        batch = {"schema":"dynacomp.network-scenarios/1", "scenarios":[{"name":"one", "initial_stock":[4,6], "replenishment":[[2]], "stock_floor":[8]}]}
        item = run_scenarios(compile_network(source), batch, exact=True)["scenarios"][0]
        self.assertEqual(item["stock_trace"], [["10"],["7"]])
        self.assertEqual(item["dispatch_trace"], [["5"]])
        self.assertEqual(item["cumulative_dispatch"], ["5"])
        self.assertEqual(item["stock_floor_violations"], [1])

    def test_zero_horizon(self):
        batch = scenarios(2, 0)
        for item in run_scenarios(compile_network(model(2)), batch, exact=True)["scenarios"]:
            self.assertEqual(item["cumulative_dispatch"], ["0"]*4)
            self.assertEqual(item["dispatch_trace"], [])
            self.assertEqual(item["mass_balance_residual"], "0")

    def test_local_defect_refines_instead_of_assuming_symmetry(self):
        source = model(4, defect=True)
        plan = compile_network(source, 16)
        self.assertGreater(plan["rank"], 4)
        self.assertTrue(plan["refinement_rounds"])
        batch = scenarios(4, 4)
        self.assertEqual(run_scenarios(plan,batch,exact=True)["scenarios"], run_scenarios(plan,batch,"full",True)["scenarios"])

    def test_no_compression_case_and_budget(self):
        source = model(4, heterogeneous=True)
        self.assertEqual(compile_network(source,16)["rank"],16)
        with self.assertRaises(OverflowError):
            compile_network(source,8)
        with self.assertRaises(ValueError):
            compile_network(model(),0)

    def test_tampered_plan_rejected(self):
        original = compile_network(model(2))
        changes = [("G",lambda p:p["G"][0].__setitem__(0,"0")),
                   ("H",lambda p:p["H"][0].__setitem__(0,"0")),
                   ("D",lambda p:p["D"][0].__setitem__(0,"0")),
                   ("source",lambda p:p["model"]["transitions"][0].__setitem__(2,"1/4")),
                   ("class",lambda p:p["classes"].__setitem__(0,1)),
                   ("rank",lambda p:p.__setitem__("rank",99999999))]
        for label, mutate in changes:
            with self.subTest(label=label):
                packet = copy.deepcopy(original)
                mutate(packet)
                with self.assertRaises(ValueError):
                    replay_network(packet)

    def test_negative_or_nonconservative_models_rejected(self):
        for value in ("-1/10","2"):
            source=model(2)
            source["transitions"][0][2]=value
            with self.assertRaises(ValueError):
                Network.load(source)
        source=model(2)
        source["inputs"][0]["distribution"][0][1]="1/7"
        with self.assertRaises(ValueError):
            Network.load(source)

    def test_scenario_validation(self):
        plan=compile_network(model(2))
        for key,value in (("initial_stock",[0]),("initial_stock",[-1]*8),("replenishment",[[0,0]]),("stock_floor",[-1]*4)):
            batch=scenarios(2,1)
            batch["scenarios"][0][key]=value
            with self.assertRaises(ValueError):
                run_scenarios(plan,batch)
        with self.assertRaises(ValueError):
            run_scenarios(plan,scenarios(2,1),"bogus")

    def test_binary64_agreement_and_conservation(self):
        plan=compile_network(model(5))
        batch=scenarios(5,20)
        full=run_scenarios(plan,batch,"full")
        reduced=run_scenarios(plan,batch)
        for a,b in zip(full["scenarios"],reduced["scenarios"]):
            for row,other in zip(a["stock_trace"]+a["dispatch_trace"],b["stock_trace"]+b["dispatch_trace"]):
                for x,y in zip(row,other):
                    self.assertAlmostEqual(x,y,places=10)
            self.assertLess(abs(b["mass_balance_residual"]),1e-9)

    def test_binary64_nonfinite_rejected(self):
        plan=compile_network(model(2))
        batch=scenarios(2,1)
        batch["scenarios"][0]["initial_stock"]=["1"+"0"*400]*8
        with self.assertRaises(OverflowError):
            run_scenarios(plan,batch)
