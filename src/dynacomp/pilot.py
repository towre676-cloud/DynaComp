"""A saved-plan workflow and matched-accuracy external workload evaluation."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import html
import json
from pathlib import Path
import platform
from statistics import median
import time

from .rational import q
from .uniform import SparseLinearModel, _domain
from .validated import validated_box, BoxKernelCache
from .stable_proposals import StableProposalContext
from .native import compose_sequential, PreparedSequential, _library


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def read_json(path):
    path = Path(path)
    if not path.exists() and path.suffix != '.gz':
        path = path.with_suffix(path.suffix+'.gz')
    if path.suffix == '.gz':
        with gzip.open(path,'rt') as stream:
            return json.load(stream)
    return json.loads(path.read_text())


def _bind(array):
    return [[str(Fraction(float(v))) for v in row] for row in array]


def workload(spec):
    """Canonical rational source and explicit finite uncertainty contract."""
    from .sparse import SparseOperator
    from .rational import matrix
    if spec.get('schema') != 'dynacomp.workload/1':
        raise ValueError('Expected dynacomp.workload/1')
    from .engineering import replay_import
    replay_import(spec)
    source = spec['model']
    if 'D' in source or 'E' in source:
        raise ValueError('Use import-mat for descriptor E and top-level feedthrough for engineering D')
    if 'A' in source:
        a = matrix(source['A']); n = len(a)
        if len(a[0]) != n:
            raise ValueError('Square A required')
        source = {'n':n,'entries':[[i,j,str(v)] for i,row in enumerate(a) for j,v in enumerate(row) if v],
                  'B':source['B'],'C':source['C']}
    model = SparseLinearModel.from_spec(source)
    domain = _domain(model,spec['initial_center'],spec['initial_radius'],spec['input_radii'],spec.get('input_centers'))
    tol = q(spec['tolerance'])
    if tol <= 0:
        raise ValueError('Positive error tolerance required')
    lower,upper = spec.get('min_rank',1),spec.get('max_rank',min(16,model.operator.n))
    if type(lower) is not int or type(upper) is not int or not 1 <= lower <= upper <= model.operator.n:
        raise ValueError('Rank range must lie in [1,n]')
    if len(domain[2]) > 4096:
        raise ValueError('Horizon must not exceed 4096')
    from .engineering import output_contract
    output_contract(spec,model,domain)
    return model,domain,tol,lower,upper


class BalancedContext:
    """One shared square-root balanced truncation factorization per search.

    This is our transparent SciPy implementation, not the SLICOT binary.
    Ordinary zero-state Hankel criteria propose; our supplied initial/input
    box and finite-horizon machine bound accept or reject both methods.
    """
    def __init__(self,model):
        import numpy as np
        from scipy.linalg import solve_discrete_lyapunov
        n = model.operator.n
        a = np.zeros((n,n))
        for i,j,v in model.operator.entries:
            a[i,j] = float(v)
        b,c = np.asarray(model.B,dtype=float),np.asarray(model.C,dtype=float)
        if any(not np.isfinite(x).all() for x in (a,b,c)):
            raise ValueError('Source conversion is nonfinite')
        if max(abs(np.linalg.eigvals(a))) >= 1:
            raise ValueError('Balanced proposal requires numerically Schur-stable source')
        def root(x):
            values,vectors = np.linalg.eigh((x+x.T)/2)
            scale = max(1,float(abs(values).max()))
            if float(values.min()) < -1e-8*scale:
                raise ValueError('Numerical Gramian is materially indefinite')
            return vectors*np.sqrt(np.maximum(values,0))
        rp = root(solve_discrete_lyapunov(a,b@b.T))
        ro = root(solve_discrete_lyapunov(a.T,c.T@c))
        u,self.s,vh = np.linalg.svd(ro.T@rp)
        self.a,self.b,self.c,self.rp,self.ro,self.u,self.vh = a,b,c,rp,ro,u,vh

    def propose(self,rank):
        import numpy as np
        if self.s[rank-1] <= 0:
            raise ValueError('Insufficient numerical Hankel rank')
        v = self.rp@self.vh[:rank].T/np.sqrt(self.s[:rank])
        w = (self.u[:,:rank].T@self.ro.T)/np.sqrt(self.s[:rank,None])
        return {k:_bind(x) for k,x in zip(('E','G','H','D'),(w,w@self.a@v,w@self.b,self.c@v))}


def compile_workload(spec,output):
    """Search both methods under one combined model+machine tolerance."""
    import numpy as np
    from threadpoolctl import threadpool_limits
    from scipy.linalg import LinAlgError
    output = Path(output); output.mkdir(parents=True,exist_ok=True)
    started = time.perf_counter()
    model,domain,tol,lower,upper = workload(spec)
    source_seconds = time.perf_counter()-started
    center,radius,radii,means = domain
    cache = BoxKernelCache(model,len(radii))
    report = {'schema':'dynacomp.pilot-plan/1','name':str(spec.get('name','custom workload')),
              'workload_sha256':digest(spec),'source_sha256':digest(model.spec()),'workload':spec,
              'source_import_seconds':source_seconds,'methods':{},
              'scope':'finite sampled rational source, independent initial/input boxes; sequential native execution includes rounding',
              'selection':'first accepted rank in supplied ascending range; not a global minimum-rank proof',
              'thread_limit':1}
    with threadpool_limits(1):
        for name,factory in (('observable',lambda:StableProposalContext(model,len(radii))),
                             ('balanced',lambda:BalancedContext(model))):
            started = time.perf_counter(); history = []; packet = None
            try:
                context = factory()
                for rank in range(lower,upper+1):
                    if name == 'observable' and rank < context.target_rank:
                        continue
                    candidate = context.propose(rank)
                    finite = validated_box(model,*(candidate[k] for k in ('E','G','H','D')),
                                           center,radius,radii,tol,input_centers=means,
                                           cache=cache,require_contraction=name=='observable')
                    # The cheaper model bracket rejects obvious failures before
                    # replay and native rounding composition are needed.
                    composed = compose_sequential(finite,tol) if finite['status']=='tolerance_met' else None
                    history.append({'rank':finite['rank'],'model_status':finite['status'],
                                    'model_lower':finite['certificate']['worst_attained_lower'],
                                    'model_upper':finite['certificate']['worst_upper'],
                                    'combined_upper':composed['combined_upper'] if composed else None,
                                    'status':composed['status'] if composed else finite['status']})
                    if composed:
                        from .engineering import output_contract
                        engineering = output_contract(spec,model,domain,composed)
                        engineering_upper = max(q(a)+q(b) for row,extra in zip(composed['combined_output_upper'],engineering['additional_output_upper']) for a,b in zip(row,extra))
                        engineering_status = 'tolerance_met' if engineering_upper <= tol and composed['status']=='tolerance_met' else 'not_certified'
                        # Do not alter the replayable native packet; the plan binds
                        # the separate engineering contract to this exact workload.
                        history[-1]['engineering_upper'] = str(engineering_upper)
                        history[-1]['status'] = engineering_status
                    if composed and composed['status']=='tolerance_met' and engineering_status=='tolerance_met':
                        packet = composed
                        break
                failure = None
            except (ValueError,OverflowError,LinAlgError) as error:
                failure = str(error)
            elapsed = time.perf_counter()-started
            result = {'accepted':packet is not None,'search_seconds':elapsed,
                      'history':history,'failure':failure}
            if packet:
                path = output/(name+'.json.gz')
                with gzip.GzipFile(filename=str(path),mode='wb',mtime=0) as stream:
                    stream.write(json.dumps(packet,separators=(',',':')).encode())
                result.update({'rank':len(packet['execution']['G']), 'combined_upper':packet['combined_upper'],
                               'model_upper':packet['source']['certificate']['worst_upper'],
                               'execution_upper':packet['execution']['worst_execution_upper'],
                               'engineering':engineering,'engineering_upper':str(engineering_upper),
                               'packet':path.name,'packet_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
            report['methods'][name] = result
    report['shared_source_passes'] = cache.source_passes
    report['search_cost_policy'] = 'observable search pays first shared source pass; balanced reuses it; per-method standalone setup also reported at evaluation'
    (output/'plan.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def _load_packet(folder,result):
    folder = Path(folder).resolve(); path = (folder/result['packet']).resolve()
    if path.parent != folder:
        raise ValueError('Plan packet path must remain inside plan directory')
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != result['packet_sha256']:
        raise ValueError('Plan artifact checksum changed')
    with gzip.open(path,'rt') as stream:
        packet = json.load(stream)
    if (result['rank'] != len(packet['execution']['G']) or result['combined_upper'] != packet['combined_upper'] or
        result['model_upper'] != packet['source']['certificate']['worst_upper'] or
        result['execution_upper'] != packet['execution']['worst_execution_upper']):
        raise ValueError('Plan summary differs from retained certificate')
    return packet


def cost_comparison(setup,full_setup,full_batch,reduced_batch,batches):
    if any(v < 0 for v in (setup,full_setup,full_batch,reduced_batch)) or type(batches) is not int or batches < 1:
        raise ValueError('Nonnegative costs and positive batch count required')
    saving = full_batch-reduced_batch
    break_even = max(0,setup-full_setup)/saving if saving > 0 else None
    full_total = full_setup+batches*full_batch
    reduced_total = setup+batches*reduced_batch
    return {'break_even_batches':break_even,'requested_batches':batches,
            'full_total_seconds':full_total,'reduced_total_seconds':reduced_total,
            'total_speedup':full_total/reduced_total if reduced_total else None,
            'warm_speedup':full_batch/reduced_batch if reduced_batch else None}


def evaluate_plan(folder,queries=256,trials=7,batches=1000,seed=20261010):
    """Time saved plans, full dense BLAS, and equal-protocol reduced methods."""
    import numpy as np
    import scipy
    from threadpoolctl import threadpool_limits
    from .validated_numeric import PreparedValidated
    folder = Path(folder)
    plan = read_json(folder/'plan.json')
    if plan.get('schema') != 'dynacomp.pilot-plan/1' or digest(plan['workload']) != plan['workload_sha256']:
        raise ValueError('Invalid saved workload plan')
    if type(queries) is not int or queries < 4 or type(trials) is not int or trials < 3 or type(batches) is not int or batches < 1:
        raise ValueError('At least four scenarios, three trials and one batch required')
    spec = plan['workload']; model,domain,tol,_,_ = workload(spec)
    if digest(model.spec()) != plan['source_sha256']:
        raise ValueError('Saved source identity differs')
    center,radius,radii,means = domain
    accepted = {name:result for name,result in plan['methods'].items() if result['accepted']}
    if not accepted:
        raise ValueError('No reduction met the target; inspect plan histories')
    packets = {name:_load_packet(folder,result) for name,result in accepted.items()}
    for packet in packets.values():
        source = packet['source']
        if (digest(source['model']) != plan['source_sha256'] or
            tuple(map(q,source['initial_center'])) != center or tuple(map(q,source['initial_radius'])) != radius or
            tuple(tuple(map(q,row)) for row in source['input_radii']) != radii or
            tuple(tuple(map(q,row)) for row in source['input_centers']) != means or q(packet['tolerance']) != tol):
            raise ValueError('Saved packet does not belong to workload contract')
    rng = np.random.default_rng(seed)
    # Random points are constructed then clipped to inward binary64 endpoints.
    from .prepared import PreparedSwitching
    def endpoints(cs,rs):
        pairs = [PreparedSwitching._float_endpoints(c-r,c+r) for c,r in zip(cs,rs)]
        return np.asarray([lo for lo,hi in pairs]),np.asarray([hi for lo,hi in pairs])
    lo,hi = endpoints(center,radius)
    initial = lo[:,None]+(hi-lo)[:,None]*rng.random((model.operator.n,queries))
    initial = np.clip(initial,lo[:,None],hi[:,None])
    inputs = np.empty((len(radii),len(model.B[0]),queries))
    for t,(cs,rs) in enumerate(zip(means,radii)):
        lower,upper = endpoints(cs,rs)
        inputs[t] = np.clip(lower[:,None]+(upper-lower)[:,None]*rng.random((len(lower),queries)),lower[:,None],upper[:,None])
    initial[:,0] = np.clip(np.asarray(list(map(float,center))),lo,hi)
    initial[:,1] = hi
    for t,(cs,rs) in enumerate(zip(means,radii)):
        lower,upper = endpoints(cs,rs)
        inputs[t,:,0] = upper
        inputs[t,:,1] = np.clip(np.asarray(list(map(float,cs))),lower,upper)
    for index,packet in enumerate(packets.values(),2):
        witness = packet['source']['certificate']['witness']
        initial[:,index] = np.clip([float(q(v)) for v in witness['x0']],lo,hi)
        for t,row in enumerate(witness['inputs']):
            lower,upper = endpoints(means[t],radii[t])
            inputs[t,:,index] = np.clip([float(q(v)) for v in row],lower,upper)
    setups = {}; runners = {}
    with threadpool_limits(1):
        started = time.perf_counter(); lib = _library(); compiler_wait = time.perf_counter()-started
        compiler_seconds = lib.build_seconds  # each standalone use pays compilation
        started = time.perf_counter()
        runners['full_blas'] = PreparedValidated(next(iter(packets.values()))['source'],'full')
        full_setup = time.perf_counter()-started
        for name,packet in packets.items():
            started = time.perf_counter(); runners[name] = PreparedSequential(packet)
            setups[name] = time.perf_counter()-started
        from .engineering import output_contract, EngineeringRunner
        contracts={name:output_contract(spec,model,domain,packet) for name,packet in packets.items()}
        for name,contract in contracts.items():
            started=time.perf_counter()
            if accepted[name].get('engineering',contract)!=contract: raise ValueError('Engineering contract does not replay')
            upper=max(q(a)+q(b) for row,extra in zip(packets[name]['combined_output_upper'],contract['additional_output_upper']) for a,b in zip(row,extra))
            if upper>tol or q(accepted[name].get('engineering_upper',upper))!=upper: raise ValueError('Engineering target does not replay')
            runners[name]=EngineeringRunner(runners[name],contract)
            setups[name]+=time.perf_counter()-started
        runners['full_blas']=EngineeringRunner(runners['full_blas'],next(iter(contracts.values())),False)
        outputs = {name:runner.simulate(initial,inputs) for name,runner in runners.items()}
        raw = {name:[] for name in runners}; keys = list(runners)
        for trial in range(trials):
            for name in (keys if trial%2==0 else list(reversed(keys))):
                started = time.perf_counter(); outputs[name] = runners[name].simulate(initial,inputs)
                raw[name].append(time.perf_counter()-started)

    full = median(raw['full_blas']); methods = {}
    for name,packet in packets.items():
        sample_error = float(abs(outputs[name]-outputs['full_blas']).max())
        # BLAS full output is a numerical diagnostic, NOT a proof oracle.
        # No tolerance pass/fail decision is made from sampled differences.
        reduced = median(raw[name]); combined = float(q(accepted[name].get('engineering_upper',packet['combined_upper'])))
        standalone_search = plan['methods'][name]['search_seconds']
        # Recompute a standalone acceptance with an empty source cache, which
        # charges each method for its own source pass, without redoing its
        # proposal search. Its incremental source cost is measured separately.
        source = packet['source']; started = time.perf_counter()
        from .validated import replay_validated
        replay_validated(source)
        standalone_replay = time.perf_counter()-started
        shared_adjustment = standalone_replay if name=='balanced' else 0.0
        setup = plan['source_import_seconds']+standalone_search+shared_adjustment+setups[name]+compiler_seconds
        methods[name] = {**plan['methods'][name],
                         'state_output_combined_upper':packet['combined_upper'],
                         'combined_upper':str(q(accepted[name].get('engineering_upper',packet['combined_upper']))),
                         'native_preparation_and_replay_seconds':setups[name],
                         'standalone_replay_adjustment_seconds':shared_adjustment,
                         'cost_policy':'conservative standalone estimate; balanced charged an extra empty-cache replay to cover the shared source pass',
                         'total_setup_seconds':setup,'median_batch_seconds':reduced,'raw_batch_seconds':raw[name],
                         'max_observed_difference_from_full_blas':sample_error,
                         'observed_difference_exceeds_bound':sample_error > combined,
                         'diagnostic_scope':'full BLAS reference has uncertified rounding; sampled difference is diagnostic only',
                         **cost_comparison(setup,plan['source_import_seconds']+full_setup,full,reduced,batches)}
    report = {'schema':'dynacomp.pilot-results/1','name':plan['name'],
              'source_sha256':plan['source_sha256'],'workload_sha256':plan['workload_sha256'],
              'states':model.operator.n,'steps':len(radii),'queries':queries,'trials':trials,'seed':seed,
              'tolerance':str(tol),'methods':methods,
              'rejected_methods':{name:result for name,result in plan['methods'].items() if not result['accepted']},
              'full_blas_setup_seconds':full_setup,
              'full_blas_median_batch_seconds':full,'full_blas_raw_batch_seconds':raw['full_blas'],
              'full_total_setup_seconds':plan['source_import_seconds']+full_setup,
              'native_compiler_seconds':compiler_seconds,'actual_compiler_wait_seconds':compiler_wait,
              'native_source_sha256':lib.source_sha256,'compiler_version':lib.compiler_version,
              'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
              'platform':platform.system()+' '+platform.machine(),'thread_limit':1,
              'measurement':'warmed alternating trials; one BLAS thread; domain checks and initial encoding included; complete output materialization',
              'cost_scope':'in-process sampled-source conversion, proposal search, independent replay, compiler and preparation; interpreter/dependency imports, file I/O, upstream MAT acquisition/sampling excluded',
              'guarantee':'accepted reduced native outputs including feedthrough versus the replayed engineering reference for every admitted stored input in boxes; finite horizon',
              'engineering_reference':next(iter(contracts.values()))['reference'],
              'continuous_sampling_bounded':next(iter(contracts.values()))['continuous_sampling_bounded'],
              'baseline':'SciPy square-root balanced truncation, accepted by same combined box/native error test; full dense NumPy BLAS reference',
              'projection':'requested total costs extrapolate measured batch medians; they are not a sustained run',
              'provenance':spec.get('provenance',{})}
    (folder/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    (folder/'report.html').write_text(render_report(report))
    return report


def render_report(report):
    esc = lambda value:html.escape(str(value))
    cards = ''
    for name,result in report['methods'].items():
        break_even = result['break_even_batches']
        cards += ('<article><h2>'+esc(name)+'</h2><p class="big">'+esc(result['rank'])+' states</p>'
                  '<p>Uniform combined error ≤ '+esc(f"{float(q(result['combined_upper'])):.6g}")+'</p>'
                  '<p>Warm speedup versus full BLAS: '+esc(f"{result['warm_speedup']:.2f}")+'×</p>'
                  '<p>Total setup: '+esc(f"{result['total_setup_seconds']:.3f}")+' s</p>'
                  '<p>Break-even: '+esc(f"{break_even:.1f} batches" if break_even is not None else 'no measured break-even')+'</p>'
                  '<p>Projected total speedup at '+esc(result['requested_batches'])+' batches: '+esc(f"{result['total_speedup']:.2f}")+'×</p>'
                  '<p class="projection" data-method="'+esc(name)+'"></p></article>')
    for name,result in report.get('rejected_methods',{}).items():
        reason = result['failure'] or 'No candidate passed within the supplied rank range.'
        cards += '<article><h2>'+esc(name)+'</h2><p class="big">Target unmet</p><p>'+esc(reason)+'</p><p>Candidate results are retained in plan.json.</p></article>'
    payload = json.dumps(report,separators=(',',':')).replace('<','\\u003c')
    interaction = ('<label>Explore total cost at <input id="batch-count" type="number" min="1" value="1000"> batches</label>'
                   '<script id="evidence" type="application/json">'+payload+'</script>'
                   '<script>const r=JSON.parse(document.getElementById("evidence").textContent);const input=document.getElementById("batch-count");'
                   'function update(){const n=Number(input.value);if(!Number.isInteger(n)||n<1)return;document.querySelectorAll(".projection").forEach(el=>{'
                   'const m=r.methods[el.dataset.method];const f=r.full_total_setup_seconds+n*r.full_blas_median_batch_seconds;'
                   'const c=m.total_setup_seconds+n*m.median_batch_seconds;el.textContent="Projected cost: full "+f.toFixed(2)+" s; reduced "+c.toFixed(2)+" s ("+(f/c).toFixed(2)+"×).";});}'
                   'input.addEventListener("input",update);update();</script>')
    data = esc(json.dumps(report,indent=2))
    return ('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width">'
            '<title>DynaComp workload results</title><style>body{font:18px system-ui;background:#101923;color:#e5eef6;max-width:1100px;margin:40px auto;padding:24px}'
            'h1{font-size:42px}main{display:flex;flex-wrap:wrap;gap:20px}article{background:#1a2b3c;border:1px solid #39546c;border-radius:16px;padding:24px;flex:1;min-width:260px}'
            '.big{font-size:36px;color:#6fe4b4}pre{font-size:13px;overflow:auto}small{color:#afc1d1}</style>'
            '<h1>'+esc(report['name'])+'</h1><p>'+esc(report['states'])+' source states · '+esc(report['steps'])+' steps · '
            +esc(report['queries'])+' scenarios per batch · target '+esc(report['tolerance'])+'</p><main>'+cards+'</main>'+interaction
            +'<p><strong>What is proved:</strong> '+esc(report['guarantee'])+'</p>'
            '<p><strong>What is measured:</strong> '+esc(report['measurement'])+'</p>'
            '<p><strong>Cost scope:</strong> '+esc(report['cost_scope'])+'</p>'
            '<p><small>Balanced truncation is a transparent SciPy implementation, not a timing of the SLICOT library. '
            'Projected totals extrapolate batch medians. Timings depend on the machine. Random traces are diagnostics; bounds cover the entire declared box. '
            'No physical measurement validation or universal performance claim.</small></p>'
            '<details><summary>Full measured evidence</summary><pre>'+data+'</pre></details></html>')


def import_mat(path,dt=None,normalize_output=False,require_bounded=False):
    """Import real nonsingular descriptors; preserve coordinates and ZOH outputs.

    Exact rational elimination; bounded sampling when enclosure limits permit.
    Outside those limits the stored-source scope is explicitly uncertified.
    """
    import numpy as np
    from scipy.io import loadmat
    from scipy.linalg import expm
    from scipy.sparse import issparse
    from .rational import matrix,encoded
    from .engineering import descriptor_transform,enclosure,replay_import
    path=Path(path); data=loadmat(path)
    def dense(value):
        value=value.toarray() if issparse(value) else np.asarray(value)
        if value.ndim!=2 or not np.isrealobj(value) or not np.isfinite(value).all():
            raise ValueError('Finite real two-dimensional MAT matrices required')
        return value
    if 'A' not in data or 'B' not in data or 'C' not in data:
        raise ValueError('Explicit A, B, C required to preserve output meaning')
    a,b,c=map(dense,(data[k] for k in ('A','B','C')));n=len(a)
    if not n or a.shape!=(n,n) or b.shape[0]!=n or not b.shape[1] or c.shape[1]!=n or not len(c):
        raise ValueError('MAT matrices have inconsistent dimensions')
    e=dense(data.get('E',np.eye(n)));f=dense(data.get('D',np.zeros((len(c),b.shape[1]))))
    if e.shape!=(n,n) or f.shape!=(len(c),b.shape[1]): raise ValueError('Descriptor or feedthrough dimensions differ')
    def exact(value):
        return matrix([[str(Fraction(v.item())) for v in row] for row in value])
    a,b,c,e,f=map(exact,(a,b,c,e,f))
    sa,sb,descriptor=descriptor_transform(a,b,e)
    divisor=max(sum(abs(v) for v in row) for row in c) if normalize_output else q(1)
    if divisor<=0: raise ValueError('Cannot normalize zero output')
    interval=None
    if dt is not None:
        interval=q(dt)
        if interval<=0 or not np.isfinite(float(interval)): raise ValueError('Positive finite sampling interval required')
        augmented=np.zeros((n+len(b[0]),n+len(b[0])))
        augmented[:n,:n]=np.asarray(sa,dtype=float)*float(interval)
        augmented[:n,n:]=np.asarray(sb,dtype=float)*float(interval)
        if not np.isfinite(augmented).all(): raise ValueError('Sampling input overflows')
        sampled=expm(augmented)
        if not np.isfinite(sampled).all(): raise ValueError('Sampling produced nonfinite coefficients')
        stored_a,stored_b=exact(sampled[:n,:n]),exact(sampled[:n,n:])
        transformation=enclosure(sa,sb,interval,stored_a,stored_b)
        if descriptor['status']!='exact':
            transformation={'status':'uncertified','reason':'continuous sampling after residual-bounded descriptor elimination is not yet enclosed'}
    else:
        stored_a,stored_b=sa,sb
        transformation=descriptor
    if require_bounded and transformation['status']=='uncertified': raise ValueError(transformation['reason'])
    binding={'schema':'dynacomp.engineering-import/1',**{k:encoded(v) for k,v in zip(('A','B','C','E','F'),(a,b,c,e,f))},
             'descriptor':descriptor,'standard_A':encoded(sa),'standard_B':encoded(sb),
             'output_divisor':str(divisor),'sampling_interval':str(interval) if interval is not None else None,'transformation':transformation}
    spec={'schema':'dynacomp.workload/1','name':path.stem,
          'model':{'A':encoded(stored_a),'B':encoded(stored_b),'C':encoded(tuple(tuple(v/divisor for v in row) for row in c))},
          'feedthrough':encoded(tuple(tuple(v/divisor for v in row) for row in f)),
          'output_input_semantics':'right-continuous-terminal-hold','import_certificate':binding,
          'initial_center':[0]*n,'initial_radius':[1]*n,'input_radii':[[1]*len(b[0]) for _ in range(20)],
          'tolerance':'1/100','min_rank':1,'max_rank':min(16,n),
          'provenance':{'mat_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'sampling_interval':binding['sampling_interval'],
                        'output_divisor':str(divisor),'scope':'exact stored coefficients; bounded ZOH sampling when replay reports bounded; sampling/physical errors otherwise unbounded; physical calibration always unbounded'}}
    replay_import(spec)
    return spec


def simulate_plan(folder,method,initial,inputs):
    """Execute a saved, independently checked plan on customer-supplied arrays."""
    folder=Path(folder)
    plan=read_json(folder/'plan.json')
    if plan.get('schema')!='dynacomp.pilot-plan/1' or digest(plan['workload'])!=plan['workload_sha256']:
        raise ValueError('Invalid saved workload plan')
    result=plan['methods'][method]
    if not result['accepted']:
        raise ValueError('Requested method did not meet the workload target')
    packet=_load_packet(folder,result)
    model,domain,tol,_,_=workload(plan['workload'])
    source=packet['source'];center,radius,radii,means=domain
    if (digest(source['model'])!=digest(model.spec()) or q(packet['tolerance'])!=tol or
        tuple(map(q,source['initial_center']))!=center or tuple(map(q,source['initial_radius']))!=radius or
        tuple(tuple(map(q,row)) for row in source['input_radii'])!=radii or
        tuple(tuple(map(q,row)) for row in source['input_centers'])!=means):
        raise ValueError('Saved packet differs from workload contract')
    from .engineering import output_contract, EngineeringRunner
    contract=output_contract(plan['workload'],model,domain,packet)
    upper=max(q(a)+q(b) for row,extra in zip(packet['combined_output_upper'],contract['additional_output_upper']) for a,b in zip(row,extra))
    if result.get('engineering',contract)!=contract or q(result.get('engineering_upper',upper))!=upper or upper>tol:
        raise ValueError('Engineering contract or tolerance does not replay')
    return EngineeringRunner(PreparedSequential(packet),contract).simulate(initial,inputs)


def main(argv=None):
    parser = argparse.ArgumentParser(description='Compile and evaluate a matched-accuracy workload')
    sub = parser.add_subparsers(dest='command',required=True)
    for name in ('compile','evaluate'):
        p = sub.add_parser(name); p.add_argument('spec'); p.add_argument('--output',required=True)
        if name=='evaluate':
            p.add_argument('--queries',type=int,default=256);p.add_argument('--trials',type=int,default=7)
            p.add_argument('--batches',type=int,default=1000)
    p = sub.add_parser('run');p.add_argument('folder');p.add_argument('--queries',type=int,default=256)
    p.add_argument('--trials',type=int,default=7);p.add_argument('--batches',type=int,default=1000)
    p = sub.add_parser('simulate');p.add_argument('folder');p.add_argument('--method',choices=['observable','balanced'],default='observable')
    p.add_argument('--initial',required=True);p.add_argument('--inputs',required=True);p.add_argument('--output',required=True)
    p = sub.add_parser('import-mat');p.add_argument('mat');p.add_argument('--dt');p.add_argument('--normalize-output',action='store_true');p.add_argument('--require-bounded',action='store_true');p.add_argument('--output',required=True)
    args = parser.parse_args(argv)
    if args.command=='simulate':
        import numpy as np
        result=simulate_plan(args.folder,args.method,np.load(args.initial,allow_pickle=False),np.load(args.inputs,allow_pickle=False))
        np.save(args.output,result);print(args.output);return
    if args.command=='import-mat':
        spec = import_mat(args.mat,args.dt,args.normalize_output,args.require_bounded)
        Path(args.output).write_text(json.dumps(spec,indent=2)+'\n');print(args.output);return
    if args.command in ('compile','evaluate'):
        plan = compile_workload(read_json(args.spec),args.output)
        if args.command=='compile':
            print(json.dumps({name:{'accepted':r['accepted'],'rank':r.get('rank')} for name,r in plan['methods'].items()}));return
        folder = args.output
    else:
        folder = args.folder
    result = evaluate_plan(folder,args.queries,args.trials,args.batches)
    print(json.dumps({name:{k:r[k] for k in ('rank','combined_upper','warm_speedup','break_even_batches','total_speedup')} for name,r in result['methods'].items()},indent=2))


if __name__=='__main__':
    main()
