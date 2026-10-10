"""Compare automatic compilation with SciPy sparse and hand-built aggregation.

Synthetic fixed-routing networks, not an industrial dataset. Timings include
all-time outputs and per-batch initial encoding. Generation/JSON I/O excluded;
construction, identity replay and numerical preparation reported separately.
"""
import argparse
import json
import os
import platform
import sys
import time
from pathlib import Path
from statistics import median
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"src"))
from dynacomp.network import compile_network,replay_network,Network
from make_demo import model


def measure(cells,queries,steps,trials,heterogeneous=False):
    import numpy as np
    import scipy
    from scipy.sparse import csr_matrix
    source=model(cells,heterogeneous=heterogeneous)
    started=time.perf_counter()
    packet=compile_network(source,max_groups=4*cells)
    replay_network(packet)
    compile_seconds=time.perf_counter()-started
    # Prepare a strong sparse baseline and target matrix, not Python loops.
    started=time.perf_counter()
    net=Network.load(source)
    from dynacomp.network_numeric import NumpyEvolution
    full_executor=NumpyEvolution(net,packet,"full")
    full_setup=time.perf_counter()-started
    started=time.perf_counter()
    reduced_executor=NumpyEvolution(net,packet,"reduced")
    reduced_setup=time.perf_counter()-started
    rng=np.random.default_rng(20261010+cells)
    initial=rng.uniform(5,15,(net.op.n,queries))
    inputs=rng.uniform(0,cells,(steps,4,queries))
    def full():
        return full_executor.simulate(initial,inputs)
    def reduced():
        return reduced_executor.simulate(initial,inputs)
    methods=[("full_sparse",full),("compiled",reduced)]
    if not heterogeneous:
        # Independent four-compartment formula: obvious optimal manual baseline.
        hand_g=np.zeros((4,4))
        hand_d=np.vstack((np.eye(4),np.diag([.02,.03,.04,.05])))
        for r in range(4):
            transfer=(3+r)/100
            hand_g[r,r]=1-transfer-(2+r)/100
            hand_g[(r+1)%4,r]=transfer
        def hand():
            z=initial.reshape(4,cells,queries).sum(axis=1)
            output=np.empty((steps+1,8,queries))
            output[0]=hand_d@z
            for t,u in enumerate(inputs):
                z=hand_g@z+u
                output[t+1]=hand_d@z
            return output
        methods.append(("hand_aggregate",hand))
    for _,function in methods:
        function()
    timings={name:[] for name,_ in methods}
    outputs={}
    for trial in range(trials):
        order=methods[trial%len(methods):]+methods[:trial%len(methods)]
        for name,function in order:
            started=time.perf_counter()
            outputs[name]=function()
            timings[name].append(time.perf_counter()-started)
    times={name:median(values) for name,values in timings.items()}
    error=float(np.max(np.abs(outputs["full_sparse"]-outputs["compiled"])))
    scale=float(np.max(np.abs(outputs["full_sparse"])))
    assert np.allclose(outputs["full_sparse"],outputs["compiled"],rtol=1e-11,atol=1e-9)
    if "hand_aggregate" in outputs:
        assert np.allclose(outputs["full_sparse"],outputs["hand_aggregate"],rtol=1e-11,atol=1e-9)
    f,r=times["full_sparse"],times["compiled"]
    return {"case":"unique per-cell dispatch rates" if heterogeneous else "heterogeneous internal routing; common regional exit fractions",
            "states":net.op.n,"compiled_states":packet["rank"],"nonzeros":full_executor.a.nnz,"queries":queries,"steps":steps,
            "compile_and_replay_seconds":compile_seconds,"full_numerical_setup_seconds":full_setup,
            "compiled_numerical_setup_seconds":reduced_setup,"median_batch_seconds":times,"raw_timings_seconds":timings,
            "execution_speedup_over_sparse":f/r,
            "first_batch_speedup_including_setup_and_compilation":(f+full_setup)/(r+compile_seconds+reduced_setup),
            "break_even_batches":max(0,compile_seconds+reduced_setup-full_setup)/(f-r) if f>r else None,
            "break_even_queries_estimate":max(0,compile_seconds+reduced_setup-full_setup)/((f-r)/queries) if f>r else None,
            "max_absolute_output_difference_all_times":error,"max_difference_scaled_to_largest_output":error/max(1,scale),
            "python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells",type=int,nargs="+",default=[32,512,2048])
    parser.add_argument("--queries",type=int,default=128)
    parser.add_argument("--steps",type=int,default=96)
    parser.add_argument("--trials",type=int,default=5)
    parser.add_argument("--output",type=Path,default=ROOT/"receipts"/"network_benchmark.json")
    args=parser.parse_args()
    if min(args.cells)<2 or min(args.queries,args.steps,args.trials)<1:
        parser.error("Positive budgets and at least two cells required")
    cases=[]
    for cells in args.cells:
        cases.append(measure(cells,args.queries,args.steps,args.trials))
        print(json.dumps(cases[-1]),flush=True)
    cases.append(measure(8,args.queries,args.steps,args.trials,heterogeneous=True))
    result={"schema":"dynacomp.network-benchmark/1","cases":cases,"platform":platform.platform(),"cpu_count":os.cpu_count(),
            "scope":"synthetic linear routing; output equality is exact algebraically but timings use binary64; automatic grouping provides no new mathematical advantage over a correctly hand-built aggregate",
            "excluded_costs":"synthetic generation, JSON I/O, scenario parsing and dependency installation; no compiled plan caching benefit assumed",
            "measurement":"five warmed rotated trials by default; initial encoding and every-time stock/dispatch outputs included; no demand caps or congestion"}
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
