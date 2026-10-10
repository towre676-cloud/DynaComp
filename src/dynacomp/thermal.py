"""Discrete heated-plate scenario workbench built on sparse uniform compilation.

Temperatures are excess above a fixed ambient; units are model-defined. Heat
stencil coefficients are per step. No continuum discretization error or real
installation validation is claimed by the algebraic output certificates.
"""
import argparse
import json
from pathlib import Path
from .rational import q,matrix,mv,addv,encoded
from .sparse import SparseOperator
from .uniform import SparseLinearModel,uniform_enrichment,replay_uniform,output_box_range
from .sparse_family import compile_sparse_family,replay_sparse_family


def plate_model(grid=12,conductance="1/8",loss="1/100"):
    if type(grid) is not int or grid<2 or grid%2:
        raise ValueError("Positive even grid size at least two required")
    k,l=q(conductance),q(loss)
    if k<0 or l<0 or 4*k+l>1:
        raise ValueError("Nonnegative conductance/loss with 4*k+loss <= 1 required")
    n=grid*grid;entries=[]
    for r in range(grid):
        for c in range(grid):
            near=[(r+dr,c+dc) for dr,dc in ((1,0),(-1,0),(0,1),(0,-1)) if 0<=r+dr<grid and 0<=c+dc<grid]
            i=r*grid+c
            entries.append((i,i,1-k*len(near)-l))
            entries.extend((i,rr*grid+cc,k) for rr,cc in near)
    middle={grid//2-1,grid//2};corners={0,grid-1}
    b=matrix([[q(int(i//grid in middle and i%grid in middle))/4] for i in range(n)])
    c=matrix([[q(int(i//grid in corners and i%grid in corners))/4 for i in range(n)]])
    return SparseLinearModel(SparseOperator.from_entries(n,entries),b,c)


def _simulate(model,packet,x0,inputs):
    e,g,h,d=map(matrix,(packet[k] for k in ('E','G','H','D')))
    x=tuple(q(v) for v in x0);z=mv(e,x)
    full,reduced=[mv(model.C,x)],[mv(d,z)]
    for u in inputs:
        x=addv(model.operator.step(x),mv(model.B,u))
        z=addv(mv(g,z),mv(h,u))
        full.append(mv(model.C,x));reduced.append(mv(d,z))
    return full,reduced


def workbench(grid=12,steps=20,tolerance="1/10",max_rank=None,threshold=21):
    if type(steps) is not int or steps<0:
        raise ValueError("Nonnegative integer horizon required")
    model=plate_model(grid)
    n=model.operator.n
    center,radius=[20]*n,[2]*n
    inputs,means=[[1]]*steps,[[1]]*steps
    approximate=uniform_enrichment(model,center,radius,inputs,tolerance,input_centers=means,max_rank=max_rank,strategy='prefix')
    assert replay_uniform(approximate)
    ranges=output_box_range(model,center,radius,inputs,means)
    classifications=[]
    for lo,hi in zip(ranges['lower'],ranges['upper']):
        classifications.append('all_below_threshold' if q(hi[0])<q(threshold) else 'all_at_or_above_threshold' if q(lo[0])>=q(threshold) else 'threshold_crossing_possible')
    scenarios=[]
    for name,initial,word in [('nominal',center,means),('heater_outage',center,[[0]]*steps),('upper_box_corner',[22]*n,[[2]]*steps)]:
        full,reduced=_simulate(model,approximate,initial,word)
        errors=[abs(a[0]-b[0]) for a,b in zip(full,reduced)]
        assert all(error<=q(bound) for error,bound in zip(errors,approximate['certificate']['bound']))
        scenarios.append({'name':name,'full_corner_mean':encoded(full),'reduced_corner_mean':encoded(reduced),
                          'max_actual_reduction_error':str(max(errors))})
    modes={}
    for name,conductance in [('low','1/16'),('high','1/8')]:
        mode=plate_model(grid,conductance)
        modes[name]=mode.spec()
    family=compile_sparse_family(modes,model.C,n)
    assert replay_sparse_family(family)
    # Test an interior, time-varying conductance word in exact arithmetic.
    e,d=matrix(family['E']),matrix(family['D'])
    x=tuple(q((7*i)%11)/10 for i in range(n));z=mv(e,x)
    low,high=plate_model(grid,'1/16'),plate_model(grid,'1/8')
    for t in range(steps):
        lam=q((t*7)%13)/12
        u=(q(int(t%3==0)),)
        xlow,xhigh=low.operator.step(x),high.operator.step(x)
        x=addv(tuple((1-lam)*a+lam*b for a,b in zip(xlow,xhigh)),mv(model.B,u))
        gl,gh=matrix(family['G']['low']),matrix(family['G']['high'])
        hl,hh=matrix(family['H']['low']),matrix(family['H']['high'])
        z=addv(tuple((1-lam)*a+lam*b for a,b in zip(mv(gl,z),mv(gh,z))),tuple((1-lam)*a+lam*b for a,b in zip(mv(hl,u),mv(hh,u))))
        assert mv(model.C,x)==mv(d,z)
    summary={'schema':'dynacomp.thermal-workbench/1','grid':grid,'full_states':n,'horizon':steps,
             'approximate_rank':approximate['rank'],'approximate_status':approximate['status'],
             'uniform_worst_reduction_error':approximate['certificate']['worst_bound'],
             'temperature_semantics':'mean of four corners, excess above ambient; per-step model units',
             'fixed_nominal_conductance':'1/8','threshold':str(q(threshold)),
             'sharp_fixed_model_temperature_ranges':ranges,'threshold_classification':classifications,'scenarios':scenarios,
             'robust_exact_family_rank':family['rank'],'robust_conductance_interval':['1/16','1/8'],
             'robust_exact_switching_checks':steps,
             'scope':'approximation/range guarantee fixes conductance at 1/8; separate exact family covers any time-varying conductance in [1/16,1/8]; no continuum or physical installation validation'}
    return approximate,family,summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grid',type=int,default=12)
    parser.add_argument('--steps',type=int,default=20)
    parser.add_argument('--tolerance',default='1/10')
    parser.add_argument('--max-rank',type=int)
    parser.add_argument('--threshold',default='21')
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    try:
        approximate,family,report=workbench(args.grid,args.steps,args.tolerance,args.max_rank,args.threshold)
        args.output_dir.mkdir(parents=True,exist_ok=True)
        for filename,packet in [('approximate_plan.json',approximate),('robust_family.json',family),('report.json',report)]:
            (args.output_dir/filename).write_text(json.dumps(packet,indent=2)+'\n')
        print(json.dumps({key:report[key] for key in ('full_states','approximate_rank','approximate_status','uniform_worst_reduction_error','robust_exact_family_rank')}))
    except (ValueError,KeyError,TypeError,OverflowError,OSError) as exc:
        parser.exit(2,f'dynacomp-thermal: {exc}\n')


if __name__=='__main__':
    main()
