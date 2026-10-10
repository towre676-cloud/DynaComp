"""Regenerate deterministic thermal workbench packets and report."""
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from dynacomp.thermal import workbench


if __name__=='__main__':
    approximate,family,report=workbench()
    for name,data in [('thermal_uniform_plan',approximate),('thermal_robust_family',family),('thermal_workbench',report)]:
        (ROOT/'receipts'/f'{name}.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({key:report[key] for key in ('full_states','approximate_rank','uniform_worst_reduction_error','robust_exact_family_rank')}))
