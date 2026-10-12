#!/usr/bin/env python3
"""Publish a source-bound browser engine and deterministic teaching scenarios."""
import hashlib, json, subprocess, sys, tomllib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from dynacomp.explorer import scenario

def export():
    examples=[]
    for key in ('warmth','supplies'):
        request={'example':key,'start':20,'strength':2,'steps':16,'output':'whole','tolerance':0}
        result=scenario(request)
        examples.append({'id':key,'request':request,'result':result})
    files={p.relative_to(ROOT/'src').as_posix():p.read_text() for p in sorted((ROOT/'src/dynacomp').glob('*.py'))}
    hashes={name:hashlib.sha256(raw.encode()).hexdigest() for name,raw in files.items()}
    web={p.name:p.read_text() for p in sorted((ROOT/'web/explorer').iterdir()) if p.is_file() and p.suffix in ('.html','.css','.js','.mjs','.svg','.webmanifest')}
    body=json.dumps({'schema':'dynacomp.explorer-payload/1','version':tomllib.loads((ROOT/'pyproject.toml').read_text())['project']['version'],'engine_files':files,'source_hashes':hashes,'web_files':web,'examples':examples,
                     'monograph':'https://github.com/towre676-cloud/DynaComp/blob/main/docs/EXPLORER_AND_LIVE_APPS.md'},sort_keys=True,separators=(',',':'))
    envelope={'schema':'dynacomp.explorer-live/1','source_parent':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'payload_json':body,'payload_sha256':hashlib.sha256(body.encode()).hexdigest()}
    (ROOT/'data').mkdir(exist_ok=True)
    (ROOT/'data/explorer_live.json').write_text(json.dumps(envelope,sort_keys=True,indent=2)+'\n')
    (ROOT/'receipts/explorer').mkdir(exist_ok=True)
    for item in examples:(ROOT/'receipts/explorer'/f"{item['id']}.json").write_text(json.dumps(item,sort_keys=True,indent=2)+'\n')
    return {'files':len(files),'examples':len(examples),'payload_sha256':envelope['payload_sha256']}
if __name__=='__main__':print(json.dumps(export()))
