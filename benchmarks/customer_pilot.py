"""Reproduce all pinned external workload evaluations without downloading data."""
import argparse
import gzip
import hashlib
import html
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from dynacomp.pilot import compile_workload,evaluate_plan,render_report


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--case',choices=['all','heat-cont','pde','build'],default='all')
    parser.add_argument('--queries',type=int,default=256)
    parser.add_argument('--trials',type=int,default=7)
    parser.add_argument('--batches',type=int,default=1000)
    parser.add_argument('--output',type=Path,default=ROOT/'receipts/customer_pilot')
    args=parser.parse_args()
    manifest=json.loads((ROOT/'benchmarks/slicot_manifest.json').read_text())
    names=['heat-cont','pde','build'] if args.case=='all' else [args.case]
    reports=[]
    for name in names:
        with gzip.open(ROOT/f'receipts/external_slicot/{name}_rank8.json.gz','rt') as stream:
            retained=json.load(stream)
        digest=hashlib.sha256(json.dumps(retained['model'],sort_keys=True,separators=(',',':')).encode()).hexdigest()
        if digest!=manifest['sampled_model_sha256'][name]:
            raise ValueError('Pinned external sampled source changed')
        spec={'schema':'dynacomp.workload/1','name':f'SLICOT {name}',
              'model':retained['model'],'tolerance':retained['tolerance'],
              **{key:retained[key] for key in ('initial_center','initial_radius','input_radii','input_centers')},
              'min_rank':1,'max_rank':16,
              'provenance':{'upstream':manifest['upstream'],'commit':manifest['commit'],
                            'mat_sha256':manifest['sha256'][name+'.mat'],
                            'sampled_model_sha256':digest,
                            'scope':'retained sampled normalized-output model; experimental boxes, not measured customer trajectories'}}
        folder=args.output/name;folder.mkdir(parents=True,exist_ok=True)
        (folder/'workload.json').write_text(json.dumps(spec,indent=2)+'\n')
        print(f'Compiling {name}',flush=True)
        compile_workload(spec,folder)
        report=evaluate_plan(folder,args.queries,args.trials,args.batches)
        # Retained dense rational sources are large; compressed JSON remains
        # directly readable by the installed pilot commands.
        for filename in ('plan.json','workload.json'):
            path=folder/filename
            with gzip.GzipFile(filename=str(path)+'.gz',mode='wb',mtime=0) as stream:
                stream.write(path.read_bytes())
            path.unlink()
        reports.append(report)
        print(json.dumps({'case':name,'methods':{key:{k:r[k] for k in ('rank','combined_upper','warm_speedup','break_even_batches','total_speedup')} for key,r in report['methods'].items()}}),flush=True)
    (args.output/'summary.json').write_text(json.dumps({'schema':'dynacomp.customer-pilot/1','cases':reports},indent=2)+'\n')
    frames=''.join('<section><iframe title="'+html.escape(report['name'],quote=True)+'" srcdoc="'+
                   html.escape(render_report(report),quote=True)+'"></iframe></section>' for report in reports)
    (args.output/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8">'
                                        '<meta name="viewport" content="width=device-width"><title>DynaComp external workload evidence</title>'
                                        '<style>body{margin:0;background:#101923;color:#e5eef6;font:18px system-ui}header{padding:32px;max-width:1100px;margin:auto}'
                                        'iframe{border:0;width:100%;height:1100px}section{border-top:2px solid #39546c}</style>'
                                        '<header><h1>DynaComp: executable external workloads</h1><p>Three pinned sampled sources. '
                                        'One declared error ceiling. Actual native rounding is included. Setup costs and failed cases remain visible.</p></header>'+frames+'</html>')


if __name__=='__main__':main()
