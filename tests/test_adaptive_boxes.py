import itertools
import copy
import random
import unittest
from fractions import Fraction as Q
from dynacomp import LinearModel, box_horizon_bound, adaptive_coordinates
from dynacomp.adaptive import coordinate_candidate
from dynacomp.adaptive import replay_adaptive
from dynacomp.adaptive import optimal_coordinates
from dynacomp.rational import matrix, mv, addv


class BoxTests(unittest.TestCase):
    def test_initial_box_and_all_input_corners(self):
        m = LinearModel.from_rows([["1/2",1,0],[0,"1/3",1],[0,0,"1/4"]], [[1],[2],[-1]], [[1,0,0],[0,1,0]])
        p = coordinate_candidate(m, [0])
        center, radius = [1,-2,0], [1,2,1]
        bound = box_horizon_bound(m,p['V'],p['E'],p['G'],p['H'],center,radius,[[2]]*3)
        for signs in itertools.product((-1,1), repeat=3):
            initial = tuple(Q(c+s*r) for c,s,r in zip(center,signs,radius))
            for words in itertools.product((-2,2),repeat=3):
                x, z = initial, mv(matrix(p['E']), initial)
                for t in range(4):
                    y, yr = mv(m.C,x), mv(matrix(p['D']),z)
                    for i in range(2):
                        self.assertLessEqual(abs(y[i]-yr[i]),Q(bound['output_bounds'][t][i]))
                    if t<3:
                        x=addv(mv(m.A,x),mv(m.B,[words[t]]))
                        z=addv(mv(matrix(p['G']),z),mv(matrix(p['H']),[words[t]]))

    def test_general_rational_lift_and_nonzero_center(self):
        from dynacomp.proposals import reduced_candidate
        m=LinearModel.from_rows([["1/2",2],[0,"1/3"]],[[1],[1]],[[1,1]])
        p=reduced_candidate(m,[[1],[2]])
        b=box_horizon_bound(m,p['V'],p['E'],p['G'],p['H'],[2,-1],[1,1],[[1]]*3)
        for initial in itertools.product((1,3),(-2,0)):
            x,z=initial,mv(matrix(p['E']),initial)
            for t in range(4):
                error=abs(mv(m.C,x)[0]-mv(matrix(p['D']),z)[0])
                self.assertLessEqual(error,Q(b['bound'][t]))
                if t<3:
                    x=addv(mv(m.A,x),mv(m.B,[1]))
                    z=addv(mv(matrix(p['G']),z),mv(matrix(p['H']),[1]))

    def test_invalid_boxes_and_encoder(self):
        m=LinearModel.from_rows([[1]],[[1]],[[1]])
        for args in (([[1]],[[2]],[0],[1],[[1]]),([[1]],[[1]],[0],[-1],[[1]]),([[1]],[[1]],[0],[1],[[-1]])):
            v,e,c,r,u=args
            with self.assertRaises(ValueError):
                box_horizon_bound(m,v,e,[[1]],[[1]],c,r,u)


class AdaptiveTests(unittest.TestCase):
    def model(self):
        return LinearModel.from_rows([["1/2",1,0],[0,"1/2",0],[0,0,7]],[[0],[1],[1]],[[1,0,0]])

    def test_missing_driver_added_irrelevant_unstable_state_removed(self):
        p=adaptive_coordinates(self.model(),[0,0,0],[1,1,1],[[1]]*4,0)
        self.assertEqual(p['status'],'tolerance_met')
        self.assertEqual(p['candidate']['selected'],[0,1])
        self.assertEqual(p['box_bound']['bound'],['0']*5)
        self.assertEqual(p['trace'][0]['chosen_coordinate'],1)

    def test_budget_exhaustion_is_not_impossibility(self):
        p=adaptive_coordinates(self.model(),[0,0,0],[1,1,1],[[1]]*4,0,max_rank=1)
        self.assertEqual(p['status'],'rank_budget_exhausted')
        self.assertEqual(len(p['trace']),1)

    def test_loose_tolerance_stops_at_lower_rank(self):
        p=adaptive_coordinates(self.model(),[0,0,0],[1,1,1],[[1]]*2,100)
        self.assertEqual(p['candidate']['selected'],[0])

    def test_full_rank_guarantees_zero_for_random_models(self):
        rng=random.Random(42)
        for _ in range(12):
            m=LinearModel.from_rows([[Q(rng.randrange(-2,3),3) for _ in range(3)] for _ in range(3)],
                                   [[1],[0],[1]],[[1,0,0]])
            p=adaptive_coordinates(m,[0]*3,[1]*3,[[1]]*3,0)
            self.assertEqual(p['status'],'tolerance_met')

    def test_coordinate_cost_changes_greedy_choice(self):
        m=LinearModel.from_rows([[0,2,1],[0,0,0],[0,0,0]],[[0],[0],[0]],[[1,0,0]])
        p=adaptive_coordinates(m,[0]*3,[1]*3,[[0]],0,coordinate_costs=[1,100,1])
        self.assertEqual(p['trace'][0]['chosen_coordinate'],2)
        self.assertEqual(p['cost'],'102')

    def test_invalid_budget_and_costs(self):
        with self.assertRaises(ValueError):
            adaptive_coordinates(self.model(),[0]*3,[1]*3,[[1]],0,max_rank=0)
        with self.assertRaises(ValueError):
            adaptive_coordinates(self.model(),[0]*3,[1]*3,[[1]],0,coordinate_costs=[1,0,1])

    def test_adaptive_packet_replay_and_mutation_rejection(self):
        packet=adaptive_coordinates(self.model(),[0]*3,[1]*3,[[1]]*3,0)
        self.assertTrue(replay_adaptive(packet))
        for target in ('bound','candidate','status','domain'):
            p=copy.deepcopy(packet)
            if target=='bound': p['box_bound']['bound'][1]='99'
            if target=='candidate': p['candidate']['G'][0][0]='99'
            if target=='status': p['status']='rank_budget_exhausted'
            if target=='domain': p['domain']['initial_radius'][0]='-1'
            with self.assertRaises(ValueError): replay_adaptive(p)

    def test_least_cost_search_avoids_unnecessary_greedy_coordinate(self):
        # Reduction of bound per cost favors coordinate 2 first, but adding
        # coordinate 1 alone already meets tolerance at a smaller total cost.
        m=LinearModel.from_rows([[0,3,1],[0,0,0],[0,0,0]],[[0],[0],[0]],[[1,0,0]])
        greedy=adaptive_coordinates(m,[0]*3,[1]*3,[[0]],1,coordinate_costs=[1,5,1])
        optimal=optimal_coordinates(m,[0]*3,[1]*3,[[0]],1,coordinate_costs=[1,5,1])
        self.assertEqual(greedy['candidate']['selected'],[0,1,2])
        self.assertEqual(optimal['candidate']['selected'],[0,1])
        self.assertEqual(optimal['cost'],'6')
        self.assertEqual(optimal['status'],'optimal_for_bound_criterion')

    def test_search_budget_and_no_accepted_subset(self):
        p=optimal_coordinates(self.model(),[0]*3,[1]*3,[[1]]*3,0,max_evaluations=1)
        self.assertEqual(p['status'],'evaluation_budget_exhausted')
        p=optimal_coordinates(self.model(),[0]*3,[1]*3,[[1]]*3,0,max_rank=1)
        self.assertEqual(p['status'],'no_bound_accepted_subset')
