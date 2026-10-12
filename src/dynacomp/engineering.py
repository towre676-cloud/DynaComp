"""Exact descriptor elimination and replayable ZOH/coefficient/output bounds.

Engineering feedthrough F is distinct from the reduction's output matrix D.
Outputs use u[t] at t<T and the last held input at the terminal sample.
"""
from fractions import Fraction as Q
from .rational import matrix, mul, inverse, identity, encoded, norm_inf, q
from .machine import _up, U, ETA, LIMIT
from .native import sequential_budget


def standardize(a,b,e):
    n=len(a)
    if e == identity(n): return a,b
    if all(not e[i][j] for i in range(n) for j in range(n) if i!=j):
        if any(not e[i][i] for i in range(n)): raise ValueError('Singular descriptor mass E')
        return tuple(tuple(v/e[i][i] for v in row) for i,row in enumerate(a)),tuple(tuple(v/e[i][i] for v in row) for i,row in enumerate(b))
    if n>32: raise ValueError('General descriptor exact elimination supports at most 32 states; diagonal E has no such limit')
    inv=inverse(e)
    return mul(inv,a),mul(inv,b)


def descriptor_transform(a,b,e,candidate=None):
    """Exact solve, or strict row dominance plus an exact residual bound.

    ||E^-1||inf <= 1/min_i(|Eii|-sum_j!=i |Eij|) by the
    maximum-component argument. This proves nonsingularity independently
    of the numerical solve used only to propose stored coefficients.
    """
    n=len(a)
    diagonal=all(not e[i][j] for i in range(n) for j in range(n) if i!=j)
    if n<=32 or diagonal:
        sa,sb=standardize(a,b,e)
        if candidate is not None and candidate!=(sa,sb): raise ValueError('Exact descriptor coefficients do not replay')
        return sa,sb,{'status':'exact','A_error_inf':'0','B_error_inf':'0'}
    margin=min(abs(e[i][i])-sum(abs(v) for j,v in enumerate(row) if i!=j) for i,row in enumerate(e))
    if margin<=0: raise ValueError('Unsupported descriptor: above 32 states require diagonal E or strict row diagonal dominance; singular systems are not approximated')
    if n>256: raise ValueError('Residual-certified descriptor solve limited to 256 states')
    if candidate is None:
        import numpy as np
        from scipy.linalg import solve
        stored=solve(np.asarray(e,dtype=float),np.asarray([ar+br for ar,br in zip(a,b)],dtype=float))
        if not np.isfinite(stored).all(): raise ValueError('Descriptor solve produced nonfinite values')
        candidate=(matrix([[str(Q(float(v))) for v in row[:n]] for row in stored]),
                   matrix([[str(Q(float(v))) for v in row[n:]] for row in stored]))
    sa,sb=candidate
    if len(sa)!=n or len(sa[0])!=n or len(sb)!=n or len(sb[0])!=len(b[0]): raise ValueError('Descriptor candidate dimensions differ')
    def residual(stored,target):
        # Sparse source rows avoid cubic work for banded engineering masses.
        columns=tuple(zip(*stored))
        return max(sum(abs(sum(v*col[j] for j,v in enumerate(row) if v)-wanted) for col,wanted in zip(columns,tr)) for row,tr in zip(e,target))
    da,db=residual(sa,a)/margin,residual(sb,b)/margin
    return sa,sb,{'status':'bounded','A_error_inf':str(_up(da)),'B_error_inf':str(_up(db)),
                  'method':'strict row diagonal dominance and exact rational solve residuals','mass_dominance_margin':str(margin)}


def enclosure(a,b,dt,stored_a,stored_b):
    """Taylor polynomial in Q; infinity-norm tail from a geometric majorant.

    For rho=||M||inf, k+2>rho bounds the tail after k by
    rho^(k+1)/(k+1)! /(1-rho/(k+2)). No floating proof arithmetic.
    """
    n,m=len(a),len(b[0]); size=n+m
    if size>12: return {'status':'uncertified','reason':'verified Taylor enclosure limited to augmented dimension 12'}
    z=matrix([[dt*v for v in a[i]+b[i]] if i<n else [0]*size for i in range(size)])
    rho=norm_inf(z)
    if rho>8: return {'status':'uncertified','reason':'verified Taylor enclosure requires augmented infinity norm <= 8'}
    term=identity(size); total=term; scalar=Q(1)
    for k in range(1,161):
        term=tuple(tuple(v/k for v in row) for row in mul(term,z))
        total=tuple(tuple(x+y for x,y in zip(row,other)) for row,other in zip(total,term))
        scalar=scalar*rho/k
        if k+2>rho:
            tail=scalar*rho/(k+1)/(1-rho/Q(k+2))
            if tail<=Q(1,1<<100): break
    else: return {'status':'uncertified','reason':'Taylor enclosure work limit'}
    error_a=norm_inf(tuple(tuple(abs(total[i][j]-stored_a[i][j]) for j in range(n)) for i in range(n)))+tail
    error_b=norm_inf(tuple(tuple(abs(total[i][n+j]-stored_b[i][j]) for j in range(m)) for i in range(n)))+tail
    return {'status':'bounded','terms':k,'tail_inf':str(_up(tail)),
            'A_error_inf':str(_up(error_a)),'B_error_inf':str(_up(error_b)),
            'method':'exact rational Taylor polynomial and geometric tail; includes stored coefficient discrepancy'}


def replay_import(spec):
    binding=spec.get('import_certificate')
    if binding is None: return None
    if binding.get('schema')!='dynacomp.engineering-import/1': raise ValueError('Unknown import certificate')
    a,b,c,e,f=map(matrix,(binding[k] for k in ('A','B','C','E','F')))
    n,m,p=len(a),len(b[0]),len(c)
    if len(a[0])!=n or len(b)!=n or len(c[0])!=n or len(e)!=n or len(e[0])!=n or len(f)!=p or len(f[0])!=m: raise ValueError('Import certificate dimensions differ')
    stored_a,stored_b=map(matrix,(spec['model'][k] for k in ('A','B')))
    candidate=tuple(map(matrix,(binding[k] for k in ('standard_A','standard_B')))) if 'standard_A' in binding else None
    sa,sb,descriptor=descriptor_transform(a,b,e,candidate)
    if descriptor!=binding.get('descriptor',descriptor): raise ValueError('Descriptor elimination bound does not replay')
    divisor=q(binding['output_divisor'])
    if divisor<=0: raise ValueError('Positive output divisor required')
    expected_c=tuple(tuple(v/divisor for v in row) for row in c)
    expected_f=tuple(tuple(v/divisor for v in row) for row in f)
    if matrix(spec['model']['C'])!=expected_c or matrix(spec['feedthrough'])!=expected_f: raise ValueError('Original outputs or units changed')
    dt=binding['sampling_interval']
    if dt is None:
        if (stored_a,stored_b)!=(sa,sb): raise ValueError('Descriptor source identities do not replay')
        certificate=descriptor
    else:
        dt=q(dt)
        if dt<=0: raise ValueError('Positive sampling interval required')
        certificate=enclosure(sa,sb,dt,stored_a,stored_b)
        if descriptor['status']!='exact':
            certificate={'status':'uncertified','reason':'continuous sampling after residual-bounded descriptor elimination is not yet enclosed'}
    if certificate!=binding['transformation']: raise ValueError('Import transformation enclosure does not replay')
    return certificate


def output_contract(spec,model,domain,packet=None):
    """Uniform imported-source + feedthrough execution budgets, time matched."""
    from .binary64 import round_bits, exact, classify
    cert=replay_import(spec)
    n,m,p=model.operator.n,len(model.B[0]),len(model.C)
    f=matrix(spec.get('feedthrough',[[0]*m for _ in range(p)]))
    if len(f)!=p or len(f[0])!=m: raise ValueError('Feedthrough must have shape [outputs,inputs]')
    if spec.get('output_input_semantics','right-continuous-terminal-hold')!='right-continuous-terminal-hold': raise ValueError('Unsupported output input semantics')
    center,radius,radii,means=domain
    if any(v for row in f for v in row) and not radii: raise ValueError('Feedthrough requires at least one held input')
    words=[[round_bits(v) for v in row] for row in f]
    if any(classify(v)=='infinity' for row in words for v in row): raise ValueError('Feedthrough loading overflows')
    loaded=tuple(tuple(exact(v) for v in row) for row in words)
    x=max(abs(c)+r for c,r in zip(center,radius)); error=Q(0)
    a=[[Q(0)]*n for _ in range(n)]
    for i,j,v in model.operator.entries:a[i][j]=v
    anorm=norm_inf(a); bnorm=norm_inf(model.B)
    da=q(cert.get('A_error_inf','0')) if cert else Q(0)
    db=q(cert.get('B_error_inf','0')) if cert else Q(0)
    bounded=cert is None or cert['status']!='uncertified'
    rows=[]; imports=[]
    for t in range(len(radii)+1):
        ub=tuple(abs(c)+r for c,r in zip(means[min(t,len(radii)-1)],radii[min(t,len(radii)-1)])) if radii else (Q(0),)*m
        out=[]; imp=[]
        for j,(row,hat,crow) in enumerate(zip(f,loaded,model.C)):
            imported=_up(sum(abs(v) for v in crow)*error); imp.append(str(imported))
            # Exact zero feedthrough preserves the original engine and budget.
            extra=Q(0)
            if any(row) and packet is not None:
                magnitude=sum(abs(v)*u for v,u in zip(hat,ub))
                ferr=sum(abs(v-w)*u for v,w,u in zip(row,hat,ub))+sequential_budget(magnitude,m)
                base=sum(abs(v)*x for v in crow)+q(packet['combined_output_upper'][t][j])
                total=base+magnitude+ferr
                add=U*total+ETA
                if total+add>LIMIT: raise ValueError('Cannot exclude feedthrough addition overflow')
                extra=_up(ferr+add)
            out.append(str(_up(imported+extra)))
        rows.append(out); imports.append(imp)
        if t<len(radii):
            u=max(abs(c)+r for c,r in zip(means[t],radii[t]))
            error=_up((anorm+da)*error+da*x+db*u)
            x=_up(anorm*x+bnorm*u)
    return {'schema':'dynacomp.engineering-output/1','feedthrough_words':words,'additional_output_upper':rows,
            'import_output_upper':imports,'engineering_import_bounded':bounded,
            'continuous_sampling_bounded':bounded if cert and spec['import_certificate']['sampling_interval'] is not None else None,
            'reference':('continuous ZOH samples' if spec['import_certificate']['sampling_interval'] is not None else 'original discrete descriptor') if cert and bounded else 'exact stored discrete coefficients',
            'semantics':'right-continuous-terminal-hold'}


class EngineeringRunner:
    """Checked base runtime plus separate sequential feedthrough dot and add."""
    def __init__(self,runner,contract,certified=True):
        import numpy as np
        from .binary64 import from_bits
        self.runner=runner; self.f=np.ascontiguousarray([[from_bits(v) for v in row] for row in contract['feedthrough_words']],dtype=np.float64); self.certified=certified
        self.f.flags.writeable=False
    def simulate(self,initial,inputs):
        import numpy as np
        output=self.runner.simulate(initial,inputs)
        if not np.any(self.f): return output
        u=np.asarray(inputs,dtype=np.float64)
        if self.certified:
            import ctypes
            from .native import _library
            u=np.ascontiguousarray(u)
            pointer=ctypes.POINTER(ctypes.c_double)
            status=_library().feedthrough(self.f.shape[1],self.f.shape[0],len(u),u.shape[2],
                *(v.ctypes.data_as(pointer) for v in (self.f,u,output)))
            if status: raise RuntimeError('Native feedthrough refused floating environment')
        else:
            for t in range(len(output)):
                output[t]=output[t]+self.f@u[min(t,len(u)-1)]
        if not np.isfinite(output).all(): raise OverflowError('Nonfinite engineering output')
        return output
