"""Reproduce analytic infinite-source receipts; no dataset download required."""
import json
from pathlib import Path
from decimal import Decimal, localcontext
from dynacomp import (compile_spectral, replay_spectral, compress_spectral_modes,
    replay_spectral_subset, compose_spectral, replay_spectral_composition,
    ramanujan_value, replay_ramanujan_value, simulate_spectral, SparseLinearModel)
from dynacomp.rational import q

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'receipts' / 'infinite_spectral'


def save(name, packet, reader):
    reader(packet)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / (name+'.json')).write_text(json.dumps(packet, indent=2)+'\n')
    key = next(k for k in ('combined_all_time_upper','combined_upper','all_time_upper') if k in packet)
    return {'rank': packet.get('rank', packet.get('modes', packet.get('finite', {}).get('rank'))), 'status': packet['status'],
            'upper': packet[key], 'upper_approx': float(q(packet[key]))}


def main():
    summary = {'schema': 'dynacomp.infinite-spectral-demo/1',
               'scope': 'analytic infinite sources, not external measured workloads; exact model arithmetic, runtime roundoff excluded'}
    r = compile_spectral({'kind':'ramanujan_log','start':2}, 32, '1/50')
    summary['ramanujan_32'] = save('ramanujan_32', r, replay_spectral)
    one = compress_spectral_modes(r, [0], '1/50')
    summary['ramanujan_one_forever'] = save('ramanujan_one_forever', one, replay_spectral_subset)
    values = simulate_spectral(one, [0]*10)
    summary['ramanujan_one_forever']['constant_output_approx'] = float(values[0][0])
    heat = compile_spectral({'kind':'heat_rod'}, 64, '1/50')
    summary['heat_64'] = save('heat_64', heat, replay_spectral)
    subset = compress_spectral_modes(heat, range(11), '1/50')
    summary['heat_eleven_forever'] = save('heat_eleven_forever', subset, replay_spectral_subset)
    cold = compile_spectral({'kind':'heat_rod','initial_amplitude':0}, 32, '1/50')
    one_heat = compress_spectral_modes(cold, [0], '1/50')
    summary['heat_one_from_rest_forever'] = save('heat_one_from_rest_forever', one_heat, replay_spectral_subset)
    with localcontext() as ctx:
        ctx.prec = 80
        # Independent Machin arctangent series (Decimal rather than rational).
        def atan(x):
            return sum((-1)**j*x**(2*j+1)/Decimal(2*j+1) for j in range(120))
        pi = 16*atan(Decimal(1)/5)-4*atan(Decimal(1)/239)
        actual = pi*pi/12-8/(pi*pi)
        summary['heat_one_from_rest_forever']['constant_input_limiting_error_approx'] = float(actual)
    try:
        from dynacomp import StableProposalContext
        context = StableProposalContext(SparseLinearModel.from_spec(heat['model']), 20)
        for rank in (1,2,3):
            candidate = context.propose(rank)
            packet = compose_spectral(heat, *(candidate[k] for k in ('E','G','H','D')), 20, '1/50', require_contraction=True)
            summary[f'heat_projected_{rank}_twenty_samples'] = save(f'heat_projected_{rank}_twenty_samples', packet, replay_spectral_composition)
    except ImportError:
        summary['optional_projection'] = 'Install numerical extras to reproduce general-direction comparison.'
    cusp = ramanujan_value('1/1000000')
    replay_ramanujan_value(cusp)
    (OUT / 'reciprocal_cusp.json').write_text(json.dumps(cusp, indent=2)+'\n')
    summary['reciprocal_cusp'] = {'branch':cusp['branch'], 'product_terms':cusp['product_terms'],
                                'interval_width_approx':float(q(cusp['interval'][1])-q(cusp['interval'][0]))}
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k: {a:b for a,b in v.items() if a != 'upper'} if isinstance(v, dict) else v for k,v in summary.items()}, indent=2))


if __name__ == '__main__':
    main()
