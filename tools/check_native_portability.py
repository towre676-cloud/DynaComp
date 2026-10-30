"""Reproduce the supplied pilot, checking native outputs against rational source.

Run with an output folder, or --compare Linux.json Windows.json. Timings are
retained by the pilot CLI but never compared between machines.
"""
import json
from fractions import Fraction as Q
from pathlib import Path
import sys

import numpy as np
from dynacomp.native import PreparedSequential
from dynacomp.pilot import compile_workload, read_json, _load_packet


def check(folder):
    spec=read_json('applications/customer_pilot/two_state.json')
    plan_path=Path(folder)/'plan.json'
    plan=read_json(plan_path) if plan_path.exists() else compile_workload(spec,folder)
    assert plan['workload']==spec
    # Fixed stored dyadics, covering endpoints and interior scenarios.
    initial=np.array([[-1,1,0.25,-0.5],[1,-1,-0.75,0.5]])
    inputs=np.array([[[-1,1,0.5,-0.25]],[[1,-1,-0.5,0.25]],
                     [[-1,1,0.25,-0.5]],[[1,-1,-0.25,0.5]]])
    a,b,c=[[[Q(v) for v in row] for row in spec['model'][key]] for key in ('A','B','C')]
    ideal=[]
    x=[[Q(float(v)) for v in row] for row in initial]
    for t in range(5):
        ideal.extend(sum(c[0][i]*x[i][j] for i in range(2)) for j in range(4))
        if t<4:
            x=[[sum(a[i][k]*x[k][j] for k in range(2))+b[i][0]*Q(float(inputs[t,0,j]))
                for j in range(4)] for i in range(2)]
    result={'schema':'dynacomp.portability-check/1','ideal':[str(v) for v in ideal], 'methods':{}}
    for name,method in plan['methods'].items():
        if not method['accepted']:
            raise AssertionError(f'Supplied pilot unexpectedly rejected {name}')
        packet=_load_packet(folder,method)
        runner=PreparedSequential(packet)
        output=runner.simulate(initial,inputs)
        assert np.array_equal(output,runner.simulate(initial,inputs))
        values=[Q(float(v)) for v in output.ravel()]
        bounds=[Q(row[0]) for row in packet['combined_output_upper'] for _ in range(4)]
        assert all(abs(v-ref)<=bound for v,ref,bound in zip(values,ideal,bounds))
        assert max(bounds)<=Q(spec['tolerance'])
        result['methods'][name]={'values':[str(v) for v in values],
                                 'bounds':[str(v) for v in bounds],
                                 'build':runner.build_info}
    Path(folder,'portability.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


def compare(left,right):
    left,right=read_json(left),read_json(right)
    assert left['ideal']==right['ideal']
    assert left['methods'].keys()==right['methods'].keys()
    for name,l in left['methods'].items():
        r=right['methods'][name]
        assert len(l['values'])==len(r['values'])==len(left['ideal'])
        for lv,rv,lb,rb in zip(l['values'],r['values'],l['bounds'],r['bounds']):
            assert abs(Q(lv)-Q(rv))<=Q(lb)+Q(rb)
    print('Linux/Windows pilot outputs agree within the sum of their certified bounds; timings not compared.')


if __name__=='__main__':
    if sys.argv[1]=='--compare':compare(*sys.argv[2:])
    else:check(sys.argv[1])
