"""Optional reusable switching scenario executor with one-time exact replay.

Two endpoint modes are interpolated per scenario/time. The certified model is
interpreted with that scalar parameter; QR and binary64 rounding are excluded
from the rational reduction-error bound. Input-box membership remains exact.
"""
import numpy as np
from scipy.sparse import csr_matrix
from fractions import Fraction
from .rational import q
from .robust import replay_switching,_models


class PreparedSwitching:
    """Validate once, then run many batches with a shared or varying initial field.

    The stored plan is privately copied so later caller mutations cannot change
    the validation domain or numerical model. No global cache or credentials.
    """
    def __init__(self,packet,engine='reduced'):
        import copy
        if engine not in ('full','reduced'):raise ValueError('Engine must be full or reduced')
        self.packet=copy.deepcopy(packet)
        if self.packet.get('schema')!='dynacomp.switching-enrichment/1':raise ValueError('Unknown switching plan')
        if engine=='reduced':replay_switching(self.packet)
        else:
            from .uniform import _domain
            source=_models(self.packet['modes'],self.packet['C'])
            _domain(next(iter(source.values())),self.packet['initial_center'],self.packet['initial_radius'],self.packet['input_radii'],self.packet['input_centers'])
        models=_models(self.packet['modes'],self.packet['C'])
        if len(models)!=2:raise ValueError('Prepared interpolation executor requires two endpoint modes')
        self.names=tuple(models)
        first=next(iter(models.values()))
        self.n,self.m,self.p=first.operator.n,len(first.B[0]),len(first.C)
        self.horizon=len(self.packet['input_radii'])
        self.engine=engine
        self.basis=None
        if engine=='reduced':
            e=np.asarray([[float(q(v)) for v in row] for row in self.packet['E']])
            self.basis=np.linalg.qr(e.T,mode='reduced')[0]
        self.a,self.b=[],[]
        for model in models.values():
            rows,cols,values=zip(*model.operator.entries) if model.operator.entries else ([],[],[])
            a=csr_matrix((np.asarray(values,dtype=float),(rows,cols)),shape=(self.n,self.n))
            b=np.asarray(model.B,dtype=float)
            self.a.append(a if self.basis is None else self.basis.T@(a@self.basis))
            self.b.append(b if self.basis is None else self.basis.T@b)
        c=np.asarray(first.C,dtype=float)
        self.c=c if self.basis is None else c@self.basis
        self.initial_intervals=[(q(c)-q(r),q(c)+q(r)) for c,r in zip(self.packet['initial_center'],self.packet['initial_radius'])]
        self.input_intervals=[[(q(c)-q(r),q(c)+q(r)) for c,r in zip(centers,radii)] for centers,radii in zip(self.packet['input_centers'],self.packet['input_radii'])]
        lower,upper=zip(*(self._float_endpoints(lo,hi) for lo,hi in self.initial_intervals))
        self.initial_lower=np.asarray(lower)[:,None];self.initial_upper=np.asarray(upper)[:,None]
        input_pairs=[[self._float_endpoints(lo,hi) for lo,hi in intervals] for intervals in self.input_intervals]
        self.input_lower=np.asarray([[lo for lo,hi in row] for row in input_pairs]).reshape(self.horizon,self.m,1)
        self.input_upper=np.asarray([[hi for lo,hi in row] for row in input_pairs]).reshape(self.horizon,self.m,1)
        self.cache_stats={'packet_replays':int(engine=='reduced'),'batches':0,'initial_fields_encoded':0}

    @staticmethod
    def _float_endpoints(low,high):
        # The smallest/largest binary64 inside each rational endpoint gives
        # exact box membership for finite binary64 values, without per-value
        # Fraction conversion. Endpoints themselves are compared exactly once.
        def rounded(value,up):
            try:result=float(value)
            except OverflowError:return np.inf if value>0 else -np.inf
            actual=Fraction(result)
            if (up and actual<value) or (not up and actual>value):
                result=float(np.nextafter(result,np.inf if up else -np.inf))
            return result
        return rounded(low,True),rounded(high,False)

    def simulate(self,initial,inputs,mixing,check_domain=True):
        """Return all-time outputs, shape [steps+1, outputs, scenarios].

        initial: [nodes], [nodes,1] or [nodes,scenarios]. inputs: [steps,channels,
        scenarios]. mixing: [steps] or [steps,scenarios], in [0,1], weighting
        the second named mode. Shared fields are projected only once. Disabling
        domain checks is for already-validated workloads, without an input claim.
        """
        initial,inputs,mixing=(np.asarray(v,dtype=float) for v in (initial,inputs,mixing))
        if initial.ndim==1:initial=initial[:,None]
        if inputs.ndim!=3 or inputs.shape[0]!=self.horizon or inputs.shape[1]!=self.m or inputs.shape[2]<1:
            raise ValueError('Wrong input batch dimensions or certified horizon')
        queries=inputs.shape[2]
        if initial.ndim!=2 or initial.shape[0]!=self.n or initial.shape[1] not in (1,queries):raise ValueError('Wrong initial field dimensions')
        if mixing.ndim==1:mixing=mixing[:,None]
        if mixing.shape not in ((self.horizon,1),(self.horizon,queries)) or any(not np.isfinite(v).all() for v in (initial,inputs,mixing)) or (mixing<0).any() or (mixing>1).any():raise ValueError('Finite matching mode interpolation in [0,1] required')
        if check_domain:
            if (initial<self.initial_lower).any() or (initial>self.initial_upper).any():
                raise ValueError('Initial field outside certified box')
            if (inputs<self.input_lower).any() or (inputs>self.input_upper).any():
                raise ValueError('Input outside certified box')
        state=initial.copy() if self.basis is None else self.basis.T@initial
        if state.shape[1]==1 and queries!=1:state=np.repeat(state,queries,axis=1)
        output=np.empty((self.horizon+1,self.p,queries));output[0]=self.c@state
        for t,u in enumerate(inputs):
            alpha=mixing[t]
            state=(1-alpha)*(self.a[0]@state+self.b[0]@u)+alpha*(self.a[1]@state+self.b[1]@u)
            output[t+1]=self.c@state
        if not np.isfinite(output).all():raise OverflowError('Nonfinite numerical output')
        self.cache_stats['batches']+=1
        if self.basis is not None:self.cache_stats['initial_fields_encoded']+=initial.shape[1]
        return output
