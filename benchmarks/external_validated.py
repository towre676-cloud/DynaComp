"""Recheck pinned external models with fast brackets and stable discovery.

Uses the already-hashed SLICOT-derived source in the v0.7.1 packets. Repeat
the original acquisition benchmark first if starting from only raw MAT files.
"""
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import platform
from statistics import median
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'benchmarks'))
from dynacomp.uniform import SparseLinearModel
from dynacomp.validated import validated_box,BoxKernelCache
from dynacomp.stable_proposals import stable_observable_proposal,StableProposalContext
from dynacomp.validated_numeric import PreparedValidated
from external_slicot import proposal,balanced


def bind(x):return [[str(Fraction(float(v))) for v in row] for row in x]


def run_case(previous,out):
    import numpy as np
    name=previous['name'];n=previous['states'];steps=previous['steps'];queries=previous['queries']
    with gzip.open(ROOT/'receipts'/'external_slicot'/(name+'_rank8.json.gz'),'rt') as f:old=json.load(f)
    source_digest=hashlib.sha256(json.dumps(old['model'],sort_keys=True,separators=(',',':')).encode()).hexdigest()
    manifest=json.loads((ROOT/'benchmarks'/'slicot_manifest.json').read_text())
    if source_digest!=manifest['sampled_model_sha256'][name]:raise ValueError('Retained sampled source differs from benchmark pin')
    source_started=time.perf_counter();model=SparseLinearModel.from_spec(old['model']);source_import=time.perf_counter()-source_started
    a=np.zeros((n,n))
    for i,j,v in model.operator.entries:a[i,j]=float(v)
    b=np.asarray(model.B,dtype=float);c=np.asarray(model.C,dtype=float)
    cache=BoxKernelCache(model,steps);packets={};costs={};rank_history=[]
    def certify(candidate,stable):
        return validated_box(model,*(candidate[k] for k in ('E','G','H','D')),old['initial_center'],old['initial_radius'],old['input_radii'],old['tolerance'],input_centers=old['input_centers'],cache=cache,require_contraction=stable)
    started=time.perf_counter();basis=proposal(a,b,c,steps,8)
    seeds=np.vstack([c,np.rint(basis.T[1:]*(2**20))/(2**20)])
    q=np.linalg.qr(seeds.T,mode='reduced')[0]
    candidate={k:bind(v) for k,v in zip(('E','G','H','D'),(q.T,q.T@a@q,q.T@b,c@q))}
    packets['orthogonal_8']=certify(candidate,False);costs['orthogonal_8']=time.perf_counter()-started
    started=time.perf_counter()
    proposal_context=StableProposalContext(model,steps)
    for rank in range(4,min(16,n)+1):
        candidate=stable_observable_proposal(model,steps,rank,context=proposal_context);packet=certify(candidate,True)
        rank_history.append({'rank':packet['rank'],'status':packet['status'],'lower':packet['certificate']['worst_attained_lower'],'upper':packet['certificate']['worst_upper'],'contraction':packet['reduced_contraction_certified']})
        if packet['status']=='tolerance_met':break
    packets['stable_search']=packet;costs['stable_search']=time.perf_counter()-started
    started=time.perf_counter();v,w,_=balanced(a,b,c,12)
    candidate={k:bind(vv) for k,vv in zip(('E','G','H','D'),(w,w@a@v,w@b,c@v))}
    packets['balanced_12']=certify(candidate,False);costs['balanced_12']=time.perf_counter()-started
    started=time.perf_counter();balanced_history=[]
    max_v,max_w,_=balanced(a,b,c,min(16,n))
    for rank in range(4,min(16,n)+1):
        v,w=max_v[:,:rank],max_w[:rank]
        candidate={k:bind(vv) for k,vv in zip(('E','G','H','D'),(w,w@a@v,w@b,c@v))}
        packet=certify(candidate,False)
        balanced_history.append({'rank':rank,'status':packet['status'],'lower':packet['certificate']['worst_attained_lower'],'upper':packet['certificate']['worst_upper']})
        if packet['status']=='tolerance_met':break
    packets['balanced_search']=packet;costs['balanced_search']=time.perf_counter()-started
    runners={};setup={}
    started=time.perf_counter();runners['full']=PreparedValidated(packets['orthogonal_8'],'full');setup['full']=time.perf_counter()-started
    for key,p in packets.items():
        started=time.perf_counter();runners[key]=PreparedValidated(p,allow_unaccepted=True);setup[key]=time.perf_counter()-started
        with gzip.GzipFile(filename=str(out/(name+'_'+key+'.json.gz')),mode='wb',mtime=0) as f:f.write(json.dumps(p,separators=(',',':')).encode())
    rng=np.random.default_rng(1010+n);initial=rng.uniform(-1,1,(n,queries));inputs=rng.uniform(-1,1,(steps,b.shape[1],queries))
    initial[:,0]=0;inputs[:,:,0]=1;initial[:,1]=1;inputs[:,:,1]=0
    for index,p in enumerate(packets.values(),2):
        w=p['certificate']['witness'];initial[:,index]=[float(Fraction(v)) for v in w['x0']]
        inputs[:,:,index]=[[float(Fraction(v)) for v in row] for row in w['inputs']]
    outputs={key:runner.simulate(initial,inputs) for key,runner in runners.items()};timings={key:[] for key in runners};keys=list(runners)
    for trial in range(5):
        for key in (keys if trial%2==0 else list(reversed(keys))):
            started=time.perf_counter();outputs[key]=runners[key].simulate(initial,inputs);timings[key].append(time.perf_counter()-started)
    full=median(timings['full']);results={}
    for index,(key,p) in enumerate(packets.items(),2):
        med=median(timings[key]);difference=abs(outputs[key]-outputs['full']);lo=float(Fraction(p['certificate']['worst_attained_lower']));hi=float(Fraction(p['certificate']['worst_upper']))
        allowance=1e-8*(1+float(abs(outputs['full']).max()))
        assert float(difference.max())<=hi+allowance
        assert float(difference[:,:,index].max())+allowance>=lo
        startup=costs[key]+setup[key]-setup['full']
        results[key]={'rank':p['rank'],'status':p['status'],'reduced_contraction_certified':p['reduced_contraction_certified'],
                      'rigorous_lower':lo,'rigorous_upper':hi,'bracket_width':hi-lo,'bits':p['certificate']['bits'],
                      'discovery_and_certificate_seconds':costs[key],'prepared_setup_with_independent_replay_seconds':setup[key],
                      'median_batch_seconds':med,'raw_batch_seconds':timings[key],'warm_speedup':full/med,
                      'break_even_batches':max(0,startup)/(full-med) if full>med else None,
                      'first_batch_speedup_including_setup':(full+setup['full'])/(med+setup[key]+costs[key]),
                      'max_observed_error':float(difference.max()),'own_witness_error':float(difference[:,:,index].max()),
                      'persistent_input_error':float(difference[:,:,0].max()),'coherent_initial_error':float(difference[:,:,1].max()),
                      'cache_stats':runners[key].cache_stats,'packet':str((out/(name+'_'+key+'.json.gz')).relative_to(ROOT))}
    old_cost=previous['methods']['dynacomp_8']['construction_seconds']+previous['methods']['dynacomp_8']['independent_replay_seconds']
    new_cost=costs['orthogonal_8']+setup['orthogonal_8']
    return {'name':name,'states':n,'steps':steps,'queries':queries,'source_mat_sha256':previous['mat_sha256'],
            'sampled_model_sha256':source_digest,
            'source_import_seconds_common':source_import,'source_passes_by_precision':cache.source_passes,
            'stable_metric_preparations':proposal_context.metric_preparations,'balanced_search_gramian_factorizations':1,
            'full_setup_seconds':setup['full'],'full_median_batch_seconds':full,'full_raw_batch_seconds':timings['full'],
            'old_exact_build_and_replay_seconds':old_cost,'new_orthogonal_discovery_certification_and_setup_seconds':new_cost,
            'startup_cost_ratio_old_over_new':old_cost/new_cost,'stable_search_history':rank_history,'balanced_search_history':balanced_history,'methods':results}


def main():
    import numpy as np
    import scipy
    from threadpoolctl import threadpool_limits
    previous=json.loads((ROOT/'receipts'/'external_slicot'/'summary.json').read_text())
    out=ROOT/'receipts'/'external_validated';out.mkdir(parents=True,exist_ok=True)
    report={'schema':'dynacomp.external-validated/1','upstream_commit':previous['upstream_commit'],'cases':[],
            'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),
            'scope':'same pinned sampled source and [-1,1] coordinate boxes as v0.7.1; finite twenty-step guarantee; binary64 evolution rounding excluded',
            'measurement':'five warmed alternating trials, one BLAS thread, cached matrices, production prepared runners with domain checks; reduced setup includes independent replay; source quantized pullbacks shared in ordered candidate search; first orthogonal candidate pays first source pass',
            'comparison':'previous exact rational versus outward dyadic acceptance of separately bound execution matrices, not bit-identical certificates; no broad superiority over balanced truncation claimed'}
    with threadpool_limits(1):
        for case in previous['cases']:
            report['cases'].append(run_case(case,out));print(json.dumps(report['cases'][-1]),flush=True)
            (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
