"""Pinned external models: exact box certification, full and balanced baselines.

Run with OPENBLAS_NUM_THREADS=1. Requires numerical extras and threadpoolctl.
Raw MAT files remain upstream; acquisition verifies every retained SHA-256.
"""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import platform
from statistics import median
import sys
import time
import urllib.request

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
PIN='98c1747030c0b1bf48d9353b2dcc9c4583867ca1'
CASES={'heat-cont':(.1,4),'pde':(.0005,4),'build':(.01,4)}


def dense(x):
    from scipy.sparse import issparse
    return x.toarray() if issparse(x) else x


def proposal(a,b,c,steps,rank):
    """Floating proposal only; exact acceptance below determines its contract."""
    import numpy as np
    pulls=[c]
    for _ in range(steps):pulls.append(pulls[-1]@a)
    q=np.linalg.qr(c.T,mode='reduced')[0]
    while q.shape[1]<rank:
        g=q.T@a@q;h=q.T@b;d=c@q;dg=d;errors=[];forcing=[]
        for ca in pulls:
            forcing.append(ca@b-dg@h)
            residual=ca-dg@q.T
            errors.append(residual);dg=dg@g
        scores=[np.sum(abs(e),axis=1)+sum((np.sum(abs(f),axis=1) for f in forcing[:t]),start=np.zeros(len(c))) for t,e in enumerate(errors)]
        t,k=max(((t,k) for t in range(len(scores)) for k in range(len(c))),key=lambda pair:scores[pair[0]][pair[1]])
        v=errors[t][k].copy();v-=q@(q.T@v);v-=q@(q.T@v)
        if np.linalg.norm(v)<1e-13:
            found=False
            for ca in pulls:
                for row in ca:
                    v=row-q@(q.T@row);v-=q@(q.T@v)
                    if np.linalg.norm(v)>1e-13:found=True;break
                if found:break
            if not found:break
        q=np.column_stack([q,v/np.linalg.norm(v)])
    return q


def balanced(a,b,c,rank):
    """Standard square-root balanced truncation of the sampled stable system."""
    import numpy as np
    from scipy.linalg import solve_discrete_lyapunov
    p=solve_discrete_lyapunov(a,b@b.T);o=solve_discrete_lyapunov(a.T,c.T@c)
    def root(x):
        values,vectors=np.linalg.eigh((x+x.T)/2)
        return vectors*np.sqrt(np.maximum(values,0))
    rp,ro=root(p),root(o)
    u,s,vh=np.linalg.svd(ro.T@rp)
    if s[rank-1]<=0:raise ValueError('Insufficient numerical Hankel rank')
    v=rp@vh[:rank].T/np.sqrt(s[:rank])
    w=(u[:,:rank].T@ro.T)/np.sqrt(s[:rank,None])
    return v,w,{'hankel_singular_values':s[:rank+1].tolist(),'left_inverse_residual':float(np.max(abs(w@v-np.eye(rank))))}


def diagnostic(a,b,c,v,w,steps):
    import numpy as np
    g=w@a@v;h=w@b;d=c@v;ca=c;dg=d;state=[];forcing=[]
    for _ in range(steps+1):
        state.append(ca-dg@w);forcing.append(ca@b-dg@h)
        ca=ca@a;dg=dg@g
    bounds=[float(max(np.sum(abs(state[t]),axis=1)+sum((np.sum(abs(f),axis=1) for f in forcing[:t]),start=np.zeros(len(c))))) for t in range(steps+1)]
    return max(bounds)


def timed(a,b,c,plans,steps,queries,trials,seed,witnesses):
    import numpy as np
    rng=np.random.default_rng(seed);n=len(a);m=b.shape[1]
    initial=rng.uniform(-1,1,(n,queries));inputs=rng.uniform(-1,1,(steps,m,queries))
    # A persistent input, a coherent initial field, then exact attaining corners.
    initial[:,0]=0;inputs[:,:,0]=1
    initial[:,1]=1;inputs[:,:,1]=0
    for k,w in enumerate(witnesses,2):
        initial[:,k]=[float(Fraction(v)) for v in w['x0']]
        inputs[:,:,k]=[[float(Fraction(v)) for v in row] for row in w['inputs']]
    def simulate(plan):
        v,w=plan
        if v is None:g,h,d,x=a,b,c,initial.copy()
        else:g,h,d,x=w@a@v,w@b,c@v,w@initial
        outputs=np.empty((steps+1,len(c),queries));outputs[0]=d@x
        # Per-batch domain checks included equally for every engine.
        if not np.isfinite(initial).all() or not np.isfinite(inputs).all() or (abs(initial)>1).any() or (abs(inputs)>1).any():raise ValueError('Outside declared binary64 domain')
        for t,u in enumerate(inputs):x=g@x+h@u;outputs[t+1]=d@x
        if not np.isfinite(outputs).all():raise OverflowError('Nonfinite output')
        return outputs
    all_plans={'full_dense':(None,None),**plans}
    outputs={key:simulate(plan) for key,plan in all_plans.items()}
    times={key:[] for key in all_plans}
    keys=list(all_plans)
    for trial in range(trials):
        for key in (keys if trial%2==0 else list(reversed(keys))):
            started=time.perf_counter();outputs[key]=simulate(all_plans[key]);times[key].append(time.perf_counter()-started)
    full=median(times['full_dense']);result={}
    for key,raw in times.items():
        med=median(raw);difference=abs(outputs[key]-outputs['full_dense'])
        result[key]={'median_batch_seconds':med,'raw_batch_seconds':raw,'execution_speedup':full/med,
                     'max_observed_error':float(difference.max()),'max_full_output':float(abs(outputs['full_dense']).max()),
                     'persistent_input_error':float(difference[:,:,0].max()),'coherent_initial_error':float(difference[:,:,1].max()),
                     'witness_errors':[float(difference[:,:,k].max()) for k in range(2,2+len(witnesses))]}
    return result


def run_case(name,dt,data_dir,out,steps,queries,trials):
    import numpy as np
    from scipy.io import loadmat
    from scipy.linalg import expm
    from dynacomp.uniform import SparseLinearModel,uniform_enrichment,replay_uniform
    from dynacomp.sparse import SparseOperator
    path=data_dir/(name+'.mat');digest=hashlib.sha256(path.read_bytes()).hexdigest()
    data=loadmat(path);a,b,c=map(dense,(data['A'],data['B'],data['C']))
    n,m=len(a),b.shape[1];mass=dense(data.get('E',np.eye(n)))
    if not np.array_equal(mass,np.eye(n)):raise ValueError('This continuous benchmark path requires identity E')
    if 'D' in data and np.any(dense(data['D'])):raise ValueError('Nonzero feedthrough unsupported in these cases')
    started=time.perf_counter();aug=np.zeros((n+m,n+m));aug[:n,:n]=a*dt;aug[:n,n:]=b*dt
    sampled=expm(aug);a,b=sampled[:n,:n],sampled[:n,n:];scale=float(abs(c).sum(axis=1).max());c=c/scale
    sampling=time.perf_counter()-started
    # No rounding of source coefficients after binary64 sampling. Every imported
    # finite binary64 coefficient is interpreted as its exact dyadic rational.
    entries=[(int(i),int(j),Fraction(float(a[i,j]))) for i,j in zip(*np.nonzero(a))]
    model=SparseLinearModel(SparseOperator.from_entries(n,entries),tuple(tuple(Fraction(float(x)) for x in row) for row in b),tuple(tuple(Fraction(float(x)) for x in row) for row in c))
    plans={};certified={};witnesses=[]
    for rank in (4,8):
        started=time.perf_counter();basis=proposal(a,b,c,steps,rank)
        seed=[tuple(Fraction(int(np.rint(value*(2**20))),2**20) for value in column) for column in basis.T[1:]]
        packet=uniform_enrichment(model,[0]*n,[1]*n,[[1]*m]*steps,'1/100',max_rank=rank,seed_rows=seed)
        construction=time.perf_counter()-started
        started=time.perf_counter();assert replay_uniform(packet);replay=time.perf_counter()-started
        artifact=out/(name+'_rank'+str(rank)+'.json.gz')
        with gzip.GzipFile(filename=str(artifact),mode='wb',mtime=0) as f:f.write(json.dumps(packet,separators=(',',':')).encode())
        e=np.asarray([[float(Fraction(x)) for x in row] for row in packet['E']]);q=np.linalg.qr(e.T,mode='reduced')[0]
        key='dynacomp_'+str(rank);plans[key]=(q,q.T);witnesses.append(packet['certificate']['attaining_witness'])
        certified[key]={'rank':packet['rank'],'status':packet['status'],'sharp_rational_worst_error_approx':float(Fraction(packet['certificate']['worst_bound'])),
                        'construction_seconds':construction,'independent_replay_seconds':replay,'packet':str(artifact.relative_to(ROOT)),
                        'basis_quantization':'dyadic denominator 2^20, shared target retained exactly; exact acceptance after quantization',
                        'binary64_box_diagnostic':diagnostic(a,b,c,q,q.T,steps)}
        print(name,key,json.dumps(certified[key]),flush=True)
        started=time.perf_counter();v,w,info=balanced(a,b,c,rank);bt_setup=time.perf_counter()-started
        bt_key='balanced_'+str(rank);plans[bt_key]=(v,w)
        certified[bt_key]={'rank':rank,'construction_seconds':bt_setup,'binary64_box_diagnostic':diagnostic(a,b,c,v,w,steps),
                           'diagnostic_is_not_a_rational_certificate':True,**info}
    # Higher-rank competitor prevents comparing only with an undersized ROM.
    started=time.perf_counter();v,w,info=balanced(a,b,c,12);bt_setup=time.perf_counter()-started
    plans['balanced_12']=(v,w)
    certified['balanced_12']={'rank':12,'construction_seconds':bt_setup,'binary64_box_diagnostic':diagnostic(a,b,c,v,w,steps),
                              'diagnostic_is_not_a_rational_certificate':True,**info}
    timings=timed(a,b,c,plans,steps,queries,trials,1010+n,witnesses)
    full=timings['full_dense']['median_batch_seconds']
    for key,entry in certified.items():
        entry.update(timings[key]);reduced=entry['median_batch_seconds'];setup=entry['construction_seconds']+entry.get('independent_replay_seconds',0)
        entry['first_batch_speedup_including_build_and_replay']=full/(setup+reduced)
        entry['break_even_batches']=setup/(full-reduced) if full>reduced else None
        if key.startswith('dynacomp'):
            assert entry['max_observed_error']<=entry['sharp_rational_worst_error_approx']+1e-8*(1+entry['max_full_output'])
            own_witness=0 if key.endswith('_4') else 1
            assert abs(entry['witness_errors'][own_witness]-entry['sharp_rational_worst_error_approx'])<1e-8*(1+entry['max_full_output'])
    return {'name':name,'mat_sha256':digest,'states':n,'inputs':m,'outputs':len(c),'sampling':'binary64 zero-order-hold augmented expm; continuous integration/rounding excluded from rational guarantee',
            'sample_interval':dt,'physical_duration':dt*steps,'steps':steps,'queries':queries,'sampling_seconds':sampling,
            'output_scale_original_units':scale,'domain':'every initial component and each input independently in [-1,1]; normalized outputs; tolerance 0.01',
            'full_dense':timings['full_dense'],'methods':certified}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir',type=Path,default=ROOT.parent/'external_slicot')
    parser.add_argument('--output-dir',type=Path,default=ROOT/'receipts'/'external_slicot')
    parser.add_argument('--steps',type=int,default=20);parser.add_argument('--queries',type=int,default=256);parser.add_argument('--trials',type=int,default=5)
    parser.add_argument('--case',choices=list(CASES),action='append')
    parser.add_argument('--download',action='store_true')
    args=parser.parse_args()
    if args.steps<1 or args.queries<4 or args.trials<1:parser.error('Positive horizon/trials and at least four scenarios required')
    # Trusted pinned benchmark source only: rational receipts can exceed Python's
    # decimal formatting cap. This CLI setting does not alter library defaults.
    if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
    manifest=json.loads((ROOT/'benchmarks'/'slicot_manifest.json').read_text())
    args.data_dir.mkdir(parents=True,exist_ok=True);args.output_dir.mkdir(parents=True,exist_ok=True)
    selected=args.case or list(CASES)
    for name in selected:
        path=args.data_dir/(name+'.mat')
        if args.download and not path.exists():
            url=f'https://raw.githubusercontent.com/SLICOT/Benchmark-ModelReduction/{PIN}/{name}.mat'
            path.write_bytes(urllib.request.urlopen(url,timeout=60).read())
        if hashlib.sha256(path.read_bytes()).hexdigest()!=manifest['sha256'][name+'.mat']:raise ValueError('Source digest mismatch: '+name)
    import numpy as np
    import scipy
    from threadpoolctl import threadpool_limits,threadpool_info
    with threadpool_limits(limits=1):
        results=[]
        for name in selected:
            results.append(run_case(name,CASES[name][0],args.data_dir,args.output_dir,args.steps,args.queries,args.trials))
            report={'schema':'dynacomp.external-slicot/1','upstream_commit':PIN,'cases':results,'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
                    'platform':platform.platform(),'threadpools':threadpool_info(),'measurement':'single BLAS thread, warmed alternating trials, identical full output histories and checked boxes; includes per-batch matrix setup; build/replay and common sampling reported separately',
                    'scope':'external published engineering models with controlled generated scenarios, not measured field data; no claim of beating balanced truncation'}
            (args.output_dir/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
