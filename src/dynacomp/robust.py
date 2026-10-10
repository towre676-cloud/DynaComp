"""Sharp box errors under arbitrary time-varying convex mode combinations.

Commuting full/projected operators compress endpoint words into a count lattice.
A signed-support dynamic program retains order-dependent forcing while avoiding
word enumeration. Noncommuting families use budgeted endpoint enumeration.
Both routes exploit multi-affinity: a worst case occurs at mode vertices.
"""
from itertools import product
from .uniform import SparseLinearModel,_domain,_project,_sign
from .rational import q,matrix,mul,mv,sub,encoded,addv
from .krylov import Echelon


def _models(modes,C):
    if not isinstance(modes,dict) or not modes or any(not isinstance(k,str) or not k for k in modes):
        raise ValueError('Nonempty named modes required')
    if any('C' in spec and matrix(spec['C'])!=matrix(C) for spec in modes.values()):
        raise ValueError('Mode targets must match shared target')
    models={name:SparseLinearModel.from_spec({**spec,'C':C}) for name,spec in modes.items()}
    first=next(iter(models.values()))
    if any(m.operator.n!=first.operator.n or len(m.B[0])!=len(first.B[0]) for m in models.values()):
        raise ValueError('Modes must share state/input dimensions')
    return models


def _sparse_product(a,b):
    rows={}
    for i,j,v in b.entries:rows.setdefault(i,[]).append((j,v))
    out={}
    for i,k,v in a.entries:
        for j,w in rows.get(k,[]):out[i,j]=out.get((i,j),q(0))+v*w
    return {key:value for key,value in out.items() if value}


def _commute(models,gs):
    names=list(models)
    for i,s in enumerate(names):
        for t in names[i+1:]:
            if _sparse_product(models[s].operator,models[t].operator)!=_sparse_product(models[t].operator,models[s].operator) or mul(gs[s],gs[t])!=mul(gs[t],gs[s]):
                return False
    return True


def _candidate(models,rows):
    first=next(iter(models.values()))
    v,e,_,_,d=_project(first,rows)
    gs={s:mul(m.operator.pullback(e),v) for s,m in models.items()}
    hs={s:mul(e,m.B) for s,m in models.items()}
    return v,e,gs,hs,d


def _dot(a,b):return sum(x*y for x,y in zip(a,b))

def _sharp_for_word(models,candidate,domain,word,output):
    v,e,gs,hs,d=candidate
    center,radius,us,means=domain
    first=next(iter(models.values()))
    a,b=(first.C[output],),(d[output],)
    kernels=[None]*len(word)
    for j in range(len(word)-1,-1,-1):
        s=word[j]
        kernels[j]=sub(mul(a,models[s].B),mul(b,hs[s]))[0]
        a=models[s].operator.pullback(a);b=mul(b,gs[s])
    k=sub(a,mul(b,e))[0]
    mean=_dot(k,center)+sum(_dot(row,means[j]) for j,row in enumerate(kernels))
    width=sum(abs(x)*r for x,r in zip(k,radius))+sum(sum(abs(x)*r for x,r in zip(row,us[j])) for j,row in enumerate(kernels))
    return abs(mean)+width,k,kernels,(-1 if mean<0 else 1)


def _witness(models,candidate,domain,time,output,word,orientation,bound):
    _,k,kernels,_=_sharp_for_word(models,candidate,domain,word,output)
    center,radius,us,means=domain
    x0=[str(c+orientation*_sign(a)*r) for c,a,r in zip(center,k,radius)]
    inputs=[]
    for j,(mean,radii) in enumerate(zip(means,us)):
        row=kernels[j] if j<time else (q(0),)*len(mean)
        inputs.append([str(c+orientation*_sign(a)*r) for c,a,r in zip(mean,row,radii)])
    names=list(models)
    return {'time':time,'output':output,'modes':list(word)+[names[0]]*(len(us)-time),
            'x0':x0,'inputs':inputs,'signed_error':str(orientation*bound)}


def switching_box_bound(models,candidate,domain,max_symbolic_states=10000):
    """Exact worst output errors for all endpoint words and box choices.

    The same suprema cover independently chosen convex weights at every time.
    Resource exhaustion raises rather than emitting an accepted partial bound.
    """
    v,e,gs,hs,d=candidate
    first=next(iter(models.values()));names=list(models);p=len(first.C)
    center,radius,us,means=domain;horizon=len(us)
    if type(max_symbolic_states) is not int or max_symbolic_states<1:
        raise ValueError('Positive symbolic-state budget required')
    exact=mul(d,e)==first.C and all(m.operator.pullback(e)==mul(gs[s],e) and mul(e,m.B)==hs[s] for s,m in models.items())
    if exact:
        envelopes=[[q(0)]*p for _ in range(horizon+1)]
        engine='exact_common_closure';states=1
        worst=(q(0),0,0,(),1)
    elif _commute(models,gs):
        engine='commuting_count_lattice'
        root=(0,)*len(names)
        # a,b are row blocks C A_suffix and D G_suffix. Commutation makes
        # them independent of suffix order; input contributions still aren't.
        cache={root:(first.C,d)};layers=[[root]]
        transitions={}
        for depth in range(horizon):
            layer=[]
            for count in layers[-1]:
                a,b=cache[count]
                for index,s in enumerate(names):
                    after=list(count);after[index]+=1;after=tuple(after)
                    if after not in cache:
                        if len(cache)>=max_symbolic_states:raise OverflowError('Switching symbolic-state budget exhausted')
                        cache[after]=(models[s].operator.pullback(a),mul(b,gs[s]));layer.append(after)
                    transitions[count,s]=sub(mul(a,models[s].B),mul(b,hs[s]))
            layers.append(layer)
        initial={count:sub(a,mul(b,e)) for count,(a,b) in cache.items()}
        states=len(cache);envelopes=[];worst=None
        for time in range(horizon+1):
            row=[]
            for output in range(p):
                best=None
                for orientation in (1,-1):
                    scores={root:q(0)};paths={root:()}
                    for depth in range(time):
                        j=time-1-depth
                        next_scores={};next_paths={}
                        for count in layers[depth]:
                            for index,s in enumerate(names):
                                after=list(count);after[index]+=1;after=tuple(after)
                                kernel=transitions[count,s][output]
                                weight=orientation*_dot(kernel,means[j])+sum(abs(a)*r for a,r in zip(kernel,us[j]))
                                score=scores[count]+weight
                                if after not in next_scores or score>next_scores[after]:
                                    next_scores[after]=score;next_paths[after]=paths[count]+(s,)
                        scores,paths=next_scores,next_paths
                    for count in layers[time]:
                        k=initial[count][output]
                        score=scores[count]+orientation*_dot(k,center)+sum(abs(a)*r for a,r in zip(k,radius))
                        if best is None or score>best[0]:best=(score,tuple(reversed(paths[count])),orientation)
                row.append(best[0])
                if worst is None or best[0]>worst[0]:worst=(best[0],time,output,best[1],best[2])
            envelopes.append(row)
    else:
        engine='endpoint_word_enumeration'
        states=sum(len(names)**depth for depth in range(horizon+1))
        if states>max_symbolic_states:raise OverflowError('Noncommuting endpoint-word budget exhausted')
        envelopes=[];worst=None
        for time in range(horizon+1):
            row=[q(0)]*p
            for word in product(names,repeat=time):
                for output in range(p):
                    bound,_,_,orientation=_sharp_for_word(models,candidate,domain,word,output)
                    row[output]=max(row[output],bound)
                    if worst is None or bound>worst[0]:worst=(bound,time,output,word,orientation)
            envelopes.append(row)
    bound,time,output,word,orientation=worst
    return {'schema':'dynacomp.switching-box-bound/1','engine':engine,'symbolic_states':states,
            'endpoint_words_at_horizon':str(len(names)**horizon),
            'output_bounds':encoded(envelopes),'bound':[str(max(row)) for row in envelopes],'worst_bound':str(bound),
            'attaining_witness':_witness(models,candidate,domain,time,output,word,orientation,bound),
            'coverage':'all independent initial/input box choices and all independently chosen convex mode weights at each time',
            'claim':'sharp componentwise/horizon output-error suprema for this switching family; rounding excluded'}


def switching_enrichment(modes,C,initial_center,initial_radius,input_radii,tolerance,input_centers=None,max_rank=None,max_symbolic_states=10000,seed_rows=None):
    models=_models(modes,C);first=next(iter(models.values()));n=first.operator.n
    domain=_domain(first,initial_center,initial_radius,input_radii,input_centers)
    tol=q(tolerance);budget=n if max_rank is None else max_rank
    if tol<0 or type(budget) is not int or not 1<=budget<=n:
        raise ValueError('Nonnegative tolerance and rank budget in [1,n] required')
    span=Echelon()
    for row in first.C:span.add(row)
    if seed_rows is not None:
        for row in matrix(seed_rows):
            if len(row)!=n:raise ValueError('Wrong seed width')
            span.add(row)
    if not len(span) or len(span)>budget:raise ValueError('Nonzero target and sufficient starting budget required')
    history=[]
    while True:
        rows=matrix([row for _,row in span.rows]);candidate=_candidate(models,rows)
        cert=switching_box_bound(models,candidate,domain,max_symbolic_states)
        history.append({'rank':len(span),'worst_bound':cert['worst_bound'],'engine':cert['engine'],'symbolic_states':cert['symbolic_states']})
        if q(cert['worst_bound'])<=tol:status='tolerance_met';break
        if len(span)>=budget:status='rank_budget_exhausted';break
        # Enrich in a general direction from a missing mode pullback. This
        # eventually reaches common exact closure, without assuming monotonicity.
        chosen=None
        for s,model in models.items():
            for i,row in enumerate(model.operator.pullback(rows)):
                if any(span.reduce(row)):
                    chosen=row;history[-1]['added_direction']={'mode':s,'row':i};break
            if chosen is not None:break
        if chosen is None:raise ArithmeticError('Common exact closure should have zero error')
        assert span.add(chosen)
    v,e,gs,hs,d=candidate;center,radius,us,means=domain
    return {'schema':'dynacomp.switching-enrichment/1','modes':{s:m.spec() for s,m in models.items()},'C':encoded(first.C),
            'V':encoded(v),'E':encoded(e),'G':{s:encoded(g) for s,g in gs.items()},'H':{s:encoded(h) for s,h in hs.items()},'D':encoded(d),
            'initial_center':[str(a) for a in center],'initial_radius':[str(a) for a in radius],
            'input_radii':encoded(us) if us else [],'input_centers':encoded(means) if means else [],
            'tolerance':str(tol),'rank':len(e),'max_rank':budget,'max_symbolic_states':max_symbolic_states,
            'status':status,'certificate':cert,'history':history,
            'claim':'general-direction projection accepted uniformly under time-varying convex mode uncertainty; rank search is heuristic'}


def replay_switching(packet):
    if packet.get('schema')!='dynacomp.switching-enrichment/1':raise ValueError('Unknown switching schema')
    models=_models(packet['modes'],packet['C']);first=next(iter(models.values()))
    e=matrix(packet['E']);v,_,gs,hs,d=_candidate(models,e)
    if encoded(v)!=packet['V'] or encoded(d)!=packet['D'] or {s:encoded(g) for s,g in gs.items()}!=packet['G'] or {s:encoded(h) for s,h in hs.items()}!=packet['H']:
        raise ValueError('Projection identities fail')
    if type(packet['rank']) is not int or type(packet['max_rank']) is not int or not 1<=packet['rank']<=packet['max_rank']<=first.operator.n or len(e)!=packet['rank']:
        raise ValueError('Invalid rank or budget')
    domain=_domain(first,packet['initial_center'],packet['initial_radius'],packet['input_radii'],packet['input_centers'])
    fresh=switching_box_bound(models,(v,e,gs,hs,d),domain,packet['max_symbolic_states'])
    if fresh!=packet['certificate']:raise ValueError('Bound or extremal witness differs')
    tol=q(packet['tolerance']);expected='tolerance_met' if q(fresh['worst_bound'])<=tol else 'rank_budget_exhausted'
    if tol<0 or packet['status']!=expected or (expected=='rank_budget_exhausted' and packet['rank']!=packet['max_rank']):raise ValueError('False acceptance status')
    witness=fresh['attaining_witness'];x=tuple(q(v) for v in witness['x0']);z=mv(e,x)
    for s,u in zip(witness['modes'][:witness['time']],witness['inputs'][:witness['time']]):
        x=addv(models[s].operator.step(x),mv(models[s].B,u));z=addv(mv(gs[s],z),mv(hs[s],u))
    error=mv(first.C,x)[witness['output']]-mv(d,z)[witness['output']]
    if error!=q(witness['signed_error']) or abs(error)!=q(fresh['worst_bound']):raise ValueError('Witness fails to attain the sharp bound')
    return True
