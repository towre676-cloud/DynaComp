"""Reproduce the robust thermal plan and its extremal mode/input experiment."""
import json
import sys
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from dynacomp.thermal import plate_model
from dynacomp.robust import switching_enrichment,replay_switching


def main():
    modes={name:plate_model(12,k).spec() for name,k in [('low','1/16'),('high','1/8')]}
    packet=switching_enrichment(modes,modes['low']['C'],[20]*144,[2]*144,[[1]]*20,'1/10',input_centers=[[1]]*20,max_rank=12)
    assert replay_switching(packet)
    (ROOT/'receipts'/'thermal_switching_plan.json').write_text(json.dumps(packet,indent=2)+'\n')
    cert=packet['certificate']
    summary={'full_states':144,'compiled_states':packet['rank'],'horizon':20,'conductance_interval':['1/16','1/8'],
             'sharp_worst_error':cert['worst_bound'],'approximately':float(Fraction(cert['worst_bound'])),
             'symbolic_states':cert['symbolic_states'],'endpoint_words':cert['endpoint_words_at_horizon'],
             'attaining_mode_word':cert['attaining_witness']['modes'],'status':packet['status'],
             'scope':'all independent temperature/heater boxes and arbitrary time-varying conductance in the declared interval; rational model reduction only'}
    (ROOT/'receipts'/'thermal_switching_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
