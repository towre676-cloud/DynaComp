"""Production thermal batch executor versus full SciPy sparse evolution.

Kernel benchmark on a synthetic stencil, not physical installation data.
Exact rational projection error and observed binary64 difference are distinct.
"""
import argparse
import json
import platform
import sys
import time
from pathlib import Path
from statistics import median
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from dynacomp.thermal import plate_model
from dynacomp.uniform import uniform_enrichment,replay_uniform


def main():
    import numpy as np
    import scipy
    from dynacomp.thermal_numeric import ThermalEvolution
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grid',type=int,default=12)
    parser.add_argument('--steps',type=int,default=20)
    parser.add_argument('--queries',type=int,default=256)
    parser.add_argument('--trials',type=int,default=5)
    parser.add_argument('--output',type=Path,default=ROOT/'receipts'/'thermal_batch_benchmark.json')
    args=parser.parse_args()
    if min(args.steps,args.queries,args.trials)<1:parser.error('Positive budgets required')
    model=plate_model(args.grid);n=model.operator.n
    started=time.perf_counter()
    packet=uniform_enrichment(model,[20]*n,[2]*n,[[1]]*args.steps,'1/10',input_centers=[[1]]*args.steps,strategy='prefix')
    replay_uniform(packet)
    construction=time.perf_counter()-started
    setup={};executors={}
    for engine in ('full','reduced'):
        started=time.perf_counter();executors[engine]=ThermalEvolution(model,packet,engine);setup[engine]=time.perf_counter()-started
    rng=np.random.default_rng(50510)
    initial=rng.uniform(18,22,(n,args.queries))
    inputs=rng.uniform(0,2,(args.steps,1,args.queries))
    for runner in executors.values():runner.simulate(initial,inputs)
    timings={engine:[] for engine in executors};outputs={}
    for trial in range(args.trials):
        order=list(executors) if trial%2==0 else list(executors)[::-1]
        for engine in order:
            started=time.perf_counter();outputs[engine]=executors[engine].simulate(initial,inputs);timings[engine].append(time.perf_counter()-started)
    medians={engine:median(values) for engine,values in timings.items()}
    differences=np.max(np.abs(outputs['full']-outputs['reduced']),axis=(1,2))
    bounds=np.asarray([float(__import__('fractions').Fraction(v)) for v in packet['certificate']['bound']])
    assert np.all(differences<=bounds+1e-9)
    f,r=medians['full'],medians['reduced']
    result={'schema':'dynacomp.thermal-benchmark/1','grid':args.grid,'full_states':n,'compiled_states':packet['rank'],'steps':args.steps,'queries':args.queries,
            'construction_and_replay_seconds':construction,'numerical_setup_seconds':setup,'median_batch_seconds':medians,'raw_timings_seconds':timings,
            'execution_speedup':f/r,'first_batch_speedup_including_construction':(f+setup['full'])/(r+setup['reduced']+construction),
            'break_even_batches_estimate':max(0,construction+setup['reduced']-setup['full'])/(f-r) if f>r else None,
            'sharp_rational_uniform_error_bound':packet['certificate']['worst_bound'],
            'max_observed_binary64_output_difference':float(max(differences)),'difference_by_time':differences.tolist(),
            'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),
            'scope':'fixed conductance synthetic stencil; rational bound excludes rounding; numerical QR basis is an approximation to the certified subspace',
            'excluded_costs':'model generation, JSON I/O, scenario parsing, installation; validation of each batch domain excluded from timed prepared executor'}
    args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
