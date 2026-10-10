"""Regenerate uniform adaptive, coefficient-family and memory demonstrations."""
import json
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from dynacomp.cli import compile_spec, replay
from dynacomp import LinearModel, adaptive_coordinates, optimal_coordinates, closure_conditions
from dynacomp.memory import simulate_memory


def main():
    summary={}
    for name in ('adaptive_driver','parameter_family','input_memory'):
        spec=json.loads((ROOT/'examples'/f'{name}.json').read_text())
        packet=compile_spec(spec)
        assert replay(packet)
        (ROOT/'receipts'/f'{name}.json').write_text(json.dumps(packet,indent=2)+'\n')
        if name=='adaptive_driver':
            summary[name]={'selected':packet['candidate']['selected'],
                           'worst_bounds':[step['worst_bound'] for step in packet['trace']],
                           'status':packet['status']}
        elif name=='parameter_family':
            conditions=closure_conditions([item['A'] for item in spec['coefficients']],[[1,0,0]])
            summary[name]={'robust_rank':len(packet['E']),'one_state_closure_conditions':conditions}
        else:
            ys=simulate_memory(packet,[2,5],[[1]]*6)
            summary[name]={'order':packet['order'],'output_coefficients':packet['output_coefficients'],
                          'input_kernels':packet['input_kernels'],'outputs':[str(v) for v in ys]}
    # Directional diffusive transport: upstream target with attenuated drivers.
    n=8
    A=[['0']*n for _ in range(n)]
    for i in range(n):
        A[i][i]='1/2'
        if i+1<n: A[i][i+1]='1/10'
    model=LinearModel.from_rows(A,[[0]]*(n-1)+[[1]],[[1]+[0]*(n-1)])
    curve=[]
    for tolerance in ('1/10','1/100','1/1000','0'):
        result=adaptive_coordinates(model,[0]*n,[1]*n,[[1]]*8,tolerance)
        assert replay(result)
        curve.append({'tolerance':tolerance,'rank':len(result['candidate']['selected']),
                      'worst_bound':result['trace'][-1]['worst_bound'],
                      'selected':result['candidate']['selected'],'status':result['status']})
        (ROOT/'receipts'/f'transport_rank{len(result["candidate"]["selected"])}.json').write_text(json.dumps(result,indent=2)+'\n')
    summary['transport_tolerance_curve']=curve
    competing=LinearModel.from_rows([[0,3,1],[0,0,0],[0,0,0]],[[0],[0],[0]],[[1,0,0]])
    options={'coordinate_costs':[1,5,1]}
    greedy=adaptive_coordinates(competing,[0]*3,[1]*3,[[0]],1,**options)
    optimal=optimal_coordinates(competing,[0]*3,[1]*3,[[0]],1,**options)
    summary['least_cost_search']={'greedy_selected':greedy['candidate']['selected'],
                                'greedy_cost':greedy['cost'],'optimal_selected':optimal['candidate']['selected'],
                                'optimal_cost':optimal['cost'],'status':optimal['status']}
    (ROOT/'receipts'/'least_cost_search.json').write_text(json.dumps(optimal,indent=2)+'\n')
    (ROOT/'receipts'/'adaptive_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
