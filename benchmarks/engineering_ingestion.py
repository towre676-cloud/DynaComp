"""Pinned public descriptor trial and retained sampled-source regressions.

Usage: PYTHONPATH=src python benchmarks/engineering_ingestion.py --models PATH
Acquire PATH with git clone https://github.com/SLICOT/Benchmark-ModelReduction
and checkout 98c1747030c0b1bf48d9353b2dcc9c4583867ca1.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
import numpy as np
from scipy.io import loadmat,savemat
from scipy.linalg import solve,solve_banded
from fractions import Fraction as Q
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from dynacomp.pilot import import_mat,compile_workload,evaluate_plan,read_json,digest,workload
from dynacomp.engineering import replay_import


def save(folder,name,value):
    folder.mkdir(parents=True,exist_ok=True)
    with gzip.GzipFile(filename=str(folder/name),mode='wb',mtime=0) as stream:
        stream.write(json.dumps(value,separators=(',',':')).encode())


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--models',type=Path)
    parser.add_argument('--output',type=Path,default=ROOT/'receipts/engineering_ingestion')
    args=parser.parse_args()
    if args.models is None:
        import subprocess,tempfile
        acquisition=tempfile.TemporaryDirectory(prefix='dynacomp-public-models-')
        args.models=Path(acquisition.name)/'models'
        subprocess.run(['git','clone','https://github.com/SLICOT/Benchmark-ModelReduction',str(args.models)],check=True)
        subprocess.run(['git','-C',str(args.models),'checkout','98c1747030c0b1bf48d9353b2dcc9c4583867ca1'],check=True)
    report={'schema':'dynacomp.engineering-benchmark/1','python':platform.python_version(),
        'upstream':'https://github.com/SLICOT/Benchmark-ModelReduction','commit':'98c1747030c0b1bf48d9353b2dcc9c4583867ca1','cases':{},'regressions':{}}
    # heat-disc is a genuinely public nonidentity, nonsingular descriptor.
    path=next(args.models.rglob('heat-disc.mat'));folder=args.output/'heat-disc'
    started=time.perf_counter();spec=import_mat(path,require_bounded=True);import_seconds=time.perf_counter()-started
    started=time.perf_counter();replay_import(spec);replay_seconds=time.perf_counter()-started
    data=loadmat(path);mass=data['E'];mass=mass.toarray() if hasattr(mass,'toarray') else mass
    a=data['A'];a=a.toarray() if hasattr(a,'toarray') else a
    b=data['B'];b=b.toarray() if hasattr(b,'toarray') else b
    started=time.perf_counter();reference=solve(mass,np.hstack((a,b)));solve_seconds=time.perf_counter()-started
    stored=np.asarray([[float(Q(v)) for v in row] for row in spec['model']['A']]);difference=float(abs(stored-reference[:,:len(a)]).max())
    bands=np.zeros((3,len(a)));bands[0,1:]=np.diag(mass,1);bands[1]=np.diag(mass);bands[2,:-1]=np.diag(mass,-1)
    independent=solve_banded((1,1),bands,np.hstack((a,b)))
    independent_difference=float(abs(reference-independent).max())
    # Keep the physical output units; the supplied boxes and tolerance are
    # experimental contracts, not calibrated operating data.
    save(folder,'workload.json.gz',spec)
    plan=compile_workload(spec,folder)
    result=evaluate_plan(folder,queries=256,trials=7,batches=1000) if any(r['accepted'] for r in plan['methods'].values()) else None
    report['cases']['heat-disc']={'mat_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'workload_sha256':digest(spec),
        'import_including_certification_seconds':import_seconds,'independent_replay_seconds':replay_seconds,
        'uncertified_reference_solve_seconds':solve_seconds,'diagnostic_max_coefficient_difference':difference,'independent_banded_vs_dense_difference':independent_difference,
        'descriptor':spec['import_certificate']['descriptor'],'methods':plan['methods'],'evaluation':result,
        'scope':'public discrete descriptor; physical output units; experimental [-1,1] boxes and twenty transitions'}
    # Explicitly synthetic public-derived stress model, not a physical benchmark.
    import tempfile
    with tempfile.TemporaryDirectory() as temporary:
        stress_path=Path(temporary)/'feedthrough-stress.mat'
        savemat(stress_path,{'A':np.diag([-1.,-2.,-3.,-4.]),'E':mass[:4,:4],
            'B':np.ones((4,1)),'C':np.array([[1.,0.,0.,0.]]),'D':np.array([[.25]])})
        started=time.perf_counter();stress=import_mat(stress_path,'1/10',require_bounded=True);stress_import=time.perf_counter()-started
        stress['max_rank']=4;sf=args.output/'feedthrough-stress';save(sf,'workload.json.gz',stress)
        sp=compile_workload(stress,sf);sr=evaluate_plan(sf,256,7,1000) if any(v['accepted'] for v in sp['methods'].values()) else None
        report['cases']['feedthrough-stress']={'scope':'synthetic dynamics and feedthrough using the public heat-disc principal mass block; not an intact public physical model',
            'import_including_certification_seconds':stress_import,'sampling':stress['import_certificate']['transformation'],'evaluation':sr,'methods':sp['methods']}
    for name in ('heat-cont','pde','build'):
        print('Regression',name,flush=True)
        old=read_json(ROOT/'receipts/customer_pilot'/name/'plan.json')
        original=read_json(ROOT/'receipts/customer_pilot'/name/'results.json')
        folder=args.output/'regression'/name
        new=compile_workload(old['workload'],folder)
        compared={}
        for method,previous in old['methods'].items():
            current=new['methods'][method]
            equal=all(current.get(k)==previous.get(k) for k in ('accepted','rank','combined_upper','model_upper','execution_upper'))
            if not equal:raise ValueError('Regression changed '+name+' '+method)
            compared[method]={'unchanged_mathematical_acceptance':equal,'rank':current.get('rank'),'combined_upper':current.get('combined_upper')}
        result=evaluate_plan(folder,256,7,1000)
        report['regressions'][name]={'workload_sha256':old['workload_sha256'],'original_results_sha256':hashlib.sha256((ROOT/'receipts/customer_pilot'/name/'results.json').read_bytes()).hexdigest(),
            'acceptance':compared,'original_warm_speedups':{k:v['warm_speedup'] for k,v in original['methods'].items()},'new_evaluation':result}
    # Singular public descriptor refusal is retained, not silently solved.
    path=next(args.models.rglob('peec.mat'));started=time.perf_counter()
    try:import_mat(path)
    except ValueError as error:report['cases']['peec']={'refused':True,'reason':str(error),'seconds':time.perf_counter()-started,'mat_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    else:raise ValueError('Unsupported public descriptor unexpectedly accepted')
    for path in args.output.rglob('plan.json'):
        save(path.parent,'plan.json.gz',json.loads(path.read_text()));path.unlink()
    (args.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'descriptor':report['cases']['heat-disc']['descriptor'],'regressions':{n:r['acceptance'] for n,r in report['regressions'].items()}},indent=2))

if __name__=='__main__':main()
