"""Reproduce GF, algebraic tails and rounded execution of an infinite-source reduction."""
import json
from pathlib import Path
from fractions import Fraction as Q
from dynacomp.cli import compile_spec, replay
from dynacomp.machine import compose_execution, Binary64Runtime
from dynacomp.binary64 import bits, dot_bits, from_bits

ROOT=Path(__file__).resolve().parents[1]


def main():
    destination=ROOT/'receipts/machine_generating'
    destination.mkdir(parents=True,exist_ok=True)
    summary={}
    for name in ('gf_fibonacci','gf_cancel','gf_delayed','gf_cauchy','binomial_sqrt','binary64_quarter','binary64_execution','machine_composition'):
        packet=compile_spec(json.loads((ROOT/'examples'/f'{name}.json').read_text()))
        assert replay(packet)
        (destination/f'{name}.json').write_text(json.dumps(packet,indent=2)+'\n')
    source=json.loads((ROOT/'receipts/infinite_spectral/heat_projected_3_twenty_samples.json').read_text())
    packet=compose_execution(source,'1/100')
    assert replay(packet)
    runtime=Binary64Runtime(packet)
    # A common uniform initial field and bounded alternating heater; the
    # finite source domain also covers independent modal initial conditions.
    n=len(packet['execution']['E'][0])
    outputs=runtime.simulate([1.0]*n,[[float(t%2)] for t in range(20)])
    summary['heat_64_to_3']={'source_model_upper':source['combined_upper'],
        'execution_upper':packet['execution']['worst_execution_upper'],
        'combined_upper':packet['combined_upper'],'combined_upper_approx':float(Q(packet['combined_upper'])),
        'execution_upper_approx':float(Q(packet['execution']['worst_execution_upper'])),
        'tolerance':packet['tolerance'],'status':packet['status'],
        'last_output':outputs[-1][0],'horizon':20,'runtime':'integer products and exact accumulation; one final RN-even per row'}
    (destination/'heat_machine.json').write_text(json.dumps(packet,indent=2)+'\n')
    summary['cancellation']={'sequential':(1e16+1.0)-1e16,
        'exact_dot_round_once':from_bits(dot_bits(list(map(bits,[1e16,1.0,-1e16])),[bits(1.0)]*3)),
        'claim':'arithmetic behavior; no runtime speedup benchmark'}
    (destination/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
