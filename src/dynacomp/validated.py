"""Rigorous finite-horizon box brackets using bounded-precision integer rows.

This accepts a supplied encoder and reduced dynamics; EV=I is not required.
Exact inequalities concern those supplied matrices, not runtime roundoff.
Outward radii track every coefficient conversion and integer dot rounding.
"""
from fractions import Fraction
from .rational import q,matrix,encoded,mul,transpose,identity,sub
from .uniform import as_sparse,_domain,_sign,SparseLinearModel
from .stability import positive_definite


def _nearest(x,scale):
    y=x*scale
    return (2*y.numerator+y.denominator)//(2*y.denominator)


def _quantize(rows,scale):
    centers=[];radii=[]
    for row in rows:
        c=[_nearest(v,scale) for v in row]
        centers.append(c);radii.append([int(Fraction(a,scale)!=v) for a,v in zip(c,row)])
    return centers,radii


def _columns(rows,scale):
    c,r=_quantize(rows,scale)
    return [[(i,a,b) for i,(a,b) in enumerate(zip(col,rad)) if a or b] for col,rad in zip(zip(*c),zip(*r))]


def _advance(row,columns,scale):
    x,r=row;centers=[];radii=[]
    for col in columns:
        dot=spread=0
        for i,a,b in col:
            dot+=x[i]*a
            spread+=abs(x[i])*b+r[i]*(abs(a)+b)
        center=(2*dot+scale)//(2*scale)
        # |exact product-center/scale| <= spread/scale^2 plus
        # the exact distance of the rounded integer dot from its midpoint.
        radius=(spread+abs(dot-center*scale)+scale-1)//scale
        centers.append(center);radii.append(radius)
    return centers,radii


def _subtract(a,b):
    return [x-y for x,y in zip(a[0],b[0])],[x+y for x,y in zip(a[1],b[1])]


class BoxKernelCache:
    """Reuse source output pullbacks and input kernels across candidate models."""
    def __init__(self,model,horizon):
        if type(horizon) is not int or horizon<0:raise ValueError('Nonnegative integer horizon required')
        self.model=as_sparse(model);self.horizon=horizon;self._cache={}
        self.source_passes=0

    def source(self,bits):
        if bits not in self._cache:
            scale=1<<bits;model=self.model
            columns=[[] for _ in range(model.operator.n)]
            for i,j,v in model.operator.entries:
                c=_nearest(v,scale);r=int(Fraction(c,scale)!=v)
                if c or r:columns[j].append((i,c,r))
            bcols=_columns(model.B,scale);c,r=_quantize(model.C,scale)
            rows=list(zip(c,r));states=[];forcing=[]
            for t in range(self.horizon+1):
                states.append(rows);forcing.append([_advance(row,bcols,scale) for row in rows])
                if t<self.horizon:rows=[_advance(row,columns,scale) for row in rows]
            self._cache[bits]=(states,forcing);self.source_passes+=1
        return self._cache[bits]


def _dimensions(model,E,G,H,D):
    e,g,h,d=map(matrix,(E,G,H,D));n=model.operator.n;r=len(g);m=len(model.B[0]);p=len(model.C)
    if len(e)!=r or len(e[0])!=n or len(g[0])!=r or len(h)!=r or len(h[0])!=m or len(d)!=p or len(d[0])!=r:
        raise ValueError('Expected E[r,n], G[r,r], H[r,m], D[p,r]')
    return e,g,h,d


def _bracket(cache,e,g,h,d,domain,bits):
    scale=1<<bits;states,forcing=cache.source(bits)
    ec,gc,hc=_columns(e,scale),_columns(g,scale),_columns(h,scale)
    dc,dr=_quantize(d,scale);rows=list(zip(dc,dr));ks=[];js=[]
    for t in range(cache.horizon+1):
        ks.append([_subtract(full,_advance(reduced,ec,scale)) for full,reduced in zip(states[t],rows)])
        js.append([_subtract(full,_advance(reduced,hc,scale)) for full,reduced in zip(forcing[t],rows)])
        if t<cache.horizon:rows=[_advance(row,gc,scale) for row in rows]
    center,radius,us,means=domain;bounds=[];lowers=[];best=None
    for t in range(cache.horizon+1):
        upper_row=[];lower_row=[]
        for output in range(len(d)):
            x,r=ks[t][output]
            mean=sum(a*c for a,c in zip(x,center));mean_error=sum(a*abs(c) for a,c in zip(r,center))
            width_hi=sum((abs(a)+b)*v for a,b,v in zip(x,r,radius))
            width_signed_lo=sum((abs(a)-b)*v for a,b,v in zip(x,r,radius))
            for j in range(t):
                k,s=js[t-1-j][output]
                mean+=sum(a*c for a,c in zip(k,means[j]));mean_error+=sum(a*abs(c) for a,c in zip(s,means[j]))
                width_hi+=sum((abs(a)+b)*v for a,b,v in zip(k,s,us[j]))
                width_signed_lo+=sum((abs(a)-b)*v for a,b,v in zip(k,s,us[j]))
            upper=(abs(mean)+mean_error+width_hi)/scale
            # This lower endpoint is attained by the sign corner of the
            # midpoint kernels. Ambiguous coefficient signs subtract radii.
            lower=max(q(0),(abs(mean)-mean_error+width_signed_lo)/scale)
            upper_row.append(upper);lower_row.append(lower)
            if best is None or lower>best[0]:best=(lower,t,output,1 if mean>=0 else -1)
        bounds.append(upper_row);lowers.append(lower_row)
    low,t,output,orientation=best;x,r=ks[t][output]
    x0=[str(c+orientation*_sign(a)*v) for c,a,v in zip(center,x,radius)]
    inputs=[]
    for j in range(cache.horizon):
        row=js[t-1-j][output][0] if j<t else [0]*len(us[j])
        inputs.append([str(c+orientation*_sign(a)*v) for c,a,v in zip(means[j],row,us[j])])
    return {'schema':'dynacomp.dyadic-box-bracket/1','bits':bits,'output_upper':encoded(bounds),'output_attained_lower':encoded(lowers),
            'worst_upper':str(max(max(row) for row in bounds)),'worst_attained_lower':str(low),
            'witness':{'time':t,'output':output,'orientation':orientation,'x0':x0,'inputs':inputs,'absolute_error_at_least':str(low)},
            'claim':'encloses the sharp worst error; lower endpoint reached or exceeded by the supplied corner; finite horizon; evolution rounding excluded'}


def euclidean_contraction(G):
    g=matrix(G)
    if len(g[0])!=len(g):raise ValueError('Square reduced operator required')
    return positive_definite(sub(identity(len(g)),mul(transpose(g),g)))


def validated_box(model,E,G,H,D,initial_center,initial_radius,input_radii,tolerance,input_centers=None,
                  precisions=(64,96,128,192),slack='1/1000000',cache=None,require_contraction=False):
    model=as_sparse(model);e,g,h,d=_dimensions(model,E,G,H,D)
    domain=_domain(model,initial_center,initial_radius,input_radii,input_centers);tol=q(tolerance);gap=q(slack)
    if tol<0 or gap<0:raise ValueError('Nonnegative tolerance and slack required')
    precisions=tuple(precisions)
    if not precisions or any(type(b) is not int or not 8<=b<=4096 for b in precisions):raise ValueError('Precision bits must be integers in [8,4096]')
    if type(require_contraction) is not bool:raise ValueError('Boolean contraction requirement expected')
    if cache is None:cache=BoxKernelCache(model,len(domain[2]))
    if cache.model!=model or cache.horizon!=len(domain[2]):raise ValueError('Cache source or horizon differs')
    stable=euclidean_contraction(g) if require_contraction else None
    history=[]
    for bits in precisions:
        certificate=_bracket(cache,e,g,h,d,domain,bits)
        lo,hi=q(certificate['worst_attained_lower']),q(certificate['worst_upper'])
        history.append({'bits':bits,'lower':str(lo),'upper':str(hi)})
        if hi<=tol:status='tolerance_met' if stable is not False else 'stability_not_certified';break
        if lo>tol:status='refuted';break
        status='undecided'
        if hi-lo<=gap:break
    center,radius,us,means=domain
    return {'schema':'dynacomp.validated-box/1','model':model.spec(),'E':encoded(e),'G':encoded(g),'H':encoded(h),'D':encoded(d),
            'initial_center':[str(v) for v in center],'initial_radius':[str(v) for v in radius],
            'input_radii':encoded(us) if us else [],'input_centers':encoded(means) if means else [],'tolerance':str(tol),
            'require_contraction':require_contraction,'reduced_contraction_certified':stable,'status':status,'rank':len(g),
            'certificate':certificate,'precision_history':history,
            'claim':'supplied encoder and reduced model accepted over stated boxes; no EV identity or physical-model validation implied'}


def replay_validated(packet):
    if packet.get('schema')!='dynacomp.validated-box/1':raise ValueError('Unknown validated-box schema')
    model=SparseLinearModel.from_spec(packet['model'])
    fresh=validated_box(model,packet['E'],packet['G'],packet['H'],packet['D'],packet['initial_center'],packet['initial_radius'],
                        packet['input_radii'],packet['tolerance'],input_centers=packet['input_centers'],
                        precisions=(packet['certificate']['bits'],),require_contraction=packet['require_contraction'])
    for key in ('certificate','status','rank','reduced_contraction_certified'):
        if fresh[key]!=packet[key]:raise ValueError('Validated certificate or acceptance differs: '+key)
    return True
