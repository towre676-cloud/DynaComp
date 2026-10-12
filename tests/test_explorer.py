import unittest
from dynacomp.explorer import scenario

class ExplorerTests(unittest.TestCase):
    def test_examples_preserve_outputs(self):
        for kind in ('warmth','supplies'):
            for output in ('whole','one'):
                result=scenario({'example':kind,'output':output,'steps':10,'strength':3})
                self.assertEqual(result['full'],result['reduced'])
                self.assertEqual(result['max_actual_difference'],'0')
                self.assertLess(result['retained_values'],result['original_values'])
    def test_approximation_meets_declared_bound(self):
        result=scenario({'example':'warmth','output':'one','tolerance':'5','steps':8})
        from fractions import Fraction
        self.assertLessEqual(Fraction(result['max_error_bound']),Fraction(5))
        self.assertLessEqual(Fraction(result['max_actual_difference']),Fraction(result['max_error_bound']))
    def test_custom(self):
        result=scenario({'example':'custom','model':{'A':[['1/2',0],[0,'1/4']],'B':[[1],[1]],'C':[[1,0]],'initial':[2,3]},'steps':4})
        self.assertEqual(result['retained_values'],1)
        self.assertEqual(result['full'],result['reduced'])
    def test_limits(self):
        with self.assertRaises(ValueError):scenario({'steps':100})
        with self.assertRaises(ValueError):scenario({'example':'unknown'})
        with self.assertRaises(ValueError):scenario({'tolerance':-1})
        with self.assertRaises(ValueError):scenario({'steps':'3/2'})
        with self.assertRaises(ValueError):scenario({'example':'custom','model':{'A':[[1]],'B':[[1]],'C':[[1]],'time':'continuous'}})
    def test_descriptor_and_feedthrough_preserve_output_coordinates(self):
        result=scenario({'example':'custom','steps':2,'strength':2,'model':{'A':[[1]],'E':[[2]],'B':[[1]],'C':[[3]],'D':[[4]],'initial':[6]}})
        self.assertEqual(result['full'],[['26'],['20'],['17']])
        self.assertEqual(result['full'],result['reduced'])
        self.assertEqual(result['descriptor']['status'],'exact')
    def test_workload_feedthrough_is_not_discarded(self):
        result=scenario({'example':'custom','steps':1,'strength':2,'model':{'model':{'A':[[1]],'B':[[1]],'C':[[1]]},'feedthrough':[[3]],'initial_center':[4]}})
        self.assertEqual(result['full'],[['10'],['12']])
