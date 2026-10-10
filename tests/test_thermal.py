import unittest
from fractions import Fraction as Q
from dynacomp.thermal import plate_model,workbench
from dynacomp import replay_uniform,replay_sparse_family


class ThermalTests(unittest.TestCase):
    def test_workbench_uniform_population_and_switching_family(self):
        approximate,family,report=workbench(6,8,'1/10')
        self.assertTrue(replay_uniform(approximate))
        self.assertTrue(replay_sparse_family(family))
        self.assertEqual(report['approximate_rank'],5)
        self.assertEqual(report['robust_exact_family_rank'],6)
        self.assertEqual(report['robust_exact_switching_checks'],8)
        self.assertLessEqual(Q(report['uniform_worst_reduction_error']),Q(1,10))
        for item in report['scenarios']:
            for t,row in enumerate(item['full_corner_mean']):
                self.assertLessEqual(Q(report['sharp_fixed_model_temperature_ranges']['lower'][t][0]),Q(row[0]))
                self.assertGreaterEqual(Q(report['sharp_fixed_model_temperature_ranges']['upper'][t][0]),Q(row[0]))

    def test_insufficient_rank_is_reported_without_false_acceptance(self):
        approximate,_,report=workbench(4,3,'1/100',max_rank=1)
        self.assertEqual(report['approximate_status'],'rank_budget_exhausted')
        self.assertTrue(replay_uniform(approximate))

    def test_physical_stencil_conserves_except_ambient_loss(self):
        model=plate_model(6)
        totals=[Q(0)]*36
        for i,j,v in model.operator.entries:
            self.assertGreaterEqual(v,0)
            totals[j]+=v
        self.assertEqual(totals,[Q(99,100)]*36)
        self.assertEqual(sum(row[0] for row in model.B),1)
        self.assertEqual(sum(model.C[0]),1)

    def test_invalid_grid_or_unsafe_timestep_rejected(self):
        for kwargs in ({'grid':3},{'grid':0},{'conductance':'1/3'},{'loss':-1}):
            with self.assertRaises(ValueError):plate_model(**kwargs)
