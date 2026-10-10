"""Reusable production switching batches, including exact membership checks.

Short approximate plan and long common-exact plan use different horizons.
Generation/JSON/installation excluded; compile/replay/preparation reported.
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
from dynacomp.robust import switching_enrichment
from dynacomp.sparse_family import compile_sparse_family


def measure(steps,queries,trials,exact):
    import numpy as np
    from dynacomp.prepared import PreparedSwitching
    modes={name:plate_model(12,k).spec() for name,k in [('low','1/16'),('high','1/8')]}
    started=time.perf_counter();seeds=None
    if exact:seeds=compile_sparse_family(modes,modes['low']['C'],144)['E']
    packet=switching_enrichment(modes,modes['low']['C'],[20]*144,[2]*144,[[1]]*steps,0 if exact else '1/10',input_centers=[[1]]*steps,max_rank=17 if exact else 12,seed_rows=seeds)
    construction=time.perf_counter()-started
    runners={};setup={}
    for engine in ('full','reduced'):
        started=time.perf_counter();runners[engine]=PreparedSwitching(packet,engine);setup[engine]=time.perf_counter()-started
    rng=np.random.default_rng(60610+steps)
    initial=np.full(144,20.);inputs=rng.uniform(0,2,(steps,1,queries));mixing=rng.uniform(0,1,(steps,queries))
    for runner in runners.values():runner.simulate(initial,inputs,mixing)
    timings={engine:[] for engine in runners};outputs={}
    for trial in range(trials):
        for engine in (('full','reduced') if trial%2==0 else ('reduced','full')):
            started=time.perf_counter();outputs[engine]=runners[engine].simulate(initial,inputs,mixing);timings[engine].append(time.perf_counter()-started)
    medians={engine:median(values) for engine,values in timings.items()}
    differences=np.max(np.abs(outputs['full']-outputs['reduced']),axis=(1,2))
    bounds=np.asarray([float(__import__('fractions').Fraction(v)) for v in packet['certificate']['bound']])
    assert (differences<=bounds+1e-9).all()
    f,r=medians['full'],medians['reduced']
    return {'case':('common exact long horizon' if steps>20 else 'common exact short horizon') if exact else 'sharp robust approximate short horizon','states':144,'rank':packet['rank'],'steps':steps,'queries':queries,
            'construction_seconds':construction,'setup_seconds':setup,'median_batch_seconds':medians,'raw_timings_seconds':timings,
            'execution_speedup_with_domain_checks':f/r,'first_batch_speedup_including_construction':(f+setup['full'])/(r+setup['reduced']+construction),
            'break_even_batches_estimate':max(0,construction+setup['reduced']-setup['full'])/(f-r) if f>r else None,
            'sharp_rational_worst_error':packet['certificate']['worst_bound'],'max_observed_binary64_difference':float(max(differences)),
            'bound_engine':packet['certificate']['engine'],'symbolic_states':packet['certificate']['symbolic_states'],
            'cache_stats':{engine:runner.cache_stats for engine,runner in runners.items()}}


def main():
    import numpy as np
    import scipy
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--queries',type=int,default=256);parser.add_argument('--trials',type=int,default=5)
    parser.add_argument('--long-steps',type=int,default=400)
    parser.add_argument('--output',type=Path,default=ROOT/'receipts'/'switching_batch_benchmark.json')
    args=parser.parse_args()
    if min(args.queries,args.trials,args.long_steps)<1:parser.error('Positive budgets required')
    cases=[]
    for steps,exact in ((20,False),(20,True),(args.long_steps,True)):
        cases.append(measure(steps,args.queries,args.trials,exact));print(json.dumps(cases[-1]),flush=True)
    result={'schema':'dynacomp.switching-benchmark/1','cases':cases,'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),
            'measurement':'five warmed alternating trials; production prepared executor, exact rational endpoint domain checks, shared initial field encoded once per batch, every-time outputs',
            'scope':'synthetic stencil and time-varying per-scenario conductance; binary64 QR/rounding excluded from rational certificates; JSON/installation/data generation excluded',
            'baseline':'full sparse source/domain validation and execution; no reduction-certificate replay charged to full baseline'}
    args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
