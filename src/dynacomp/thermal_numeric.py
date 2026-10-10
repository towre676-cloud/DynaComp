"""Optional batched thermal execution in a numerically orthogonal basis.

QR changes coordinates within the rationally certified subspace. Binary64 QR
and evolution introduce rounding not covered by the exact algebraic bound.
"""
import numpy as np
from scipy.sparse import csr_matrix
from .rational import q


class ThermalEvolution:
    """Executor for a validated model/plan pair; columns are scenarios."""
    def __init__(self,model,packet,engine='reduced'):
        if engine not in ('full','reduced'):
            raise ValueError("Engine must be full or reduced")
        self.n=model.operator.n;self.m=len(model.B[0]);self.p=len(model.C)
        rows,cols,values=zip(*model.operator.entries) if model.operator.entries else ([],[],[])
        a=csr_matrix((np.asarray(values,dtype=float),(rows,cols)),shape=(self.n,self.n))
        b=np.asarray(model.B,dtype=float);c=np.asarray(model.C,dtype=float)
        if engine=='full':
            self.q=None;self.a,self.b,self.c=a,b,c
        else:
            e=np.asarray([[float(q(v)) for v in row] for row in packet['E']])
            self.q=np.linalg.qr(e.T,mode='reduced')[0]
            self.a=self.q.T@(a@self.q)
            self.b=self.q.T@b;self.c=c@self.q

    def simulate(self,initial,inputs):
        initial,inputs=np.asarray(initial,dtype=float),np.asarray(inputs,dtype=float)
        if initial.ndim!=2 or initial.shape[0]!=self.n or initial.shape[1]<1 or inputs.ndim!=3 or inputs.shape[1:]!=(self.m,initial.shape[1]) or not np.isfinite(initial).all() or not np.isfinite(inputs).all():
            raise ValueError("Finite initial [nodes,scenarios] and input [steps,channels,scenarios] arrays required")
        state=initial.copy() if self.q is None else self.q.T@initial
        output=np.empty((len(inputs)+1,self.p,initial.shape[1]))
        output[0]=self.c@state
        for t,u in enumerate(inputs):
            state=self.a@state+self.b@u
            output[t+1]=self.c@state
        if not np.isfinite(output).all():
            raise OverflowError("Nonfinite numerical outputs")
        return output


def run_thermal_batch(packet,initial,inputs,engine='reduced'):
    """Public convenience path validates the packet and its box domain."""
    from .uniform import SparseLinearModel,replay_uniform
    replay_uniform(packet)
    model=SparseLinearModel.from_spec(packet['model'])
    initial,inputs=np.asarray(initial,dtype=float),np.asarray(inputs,dtype=float)
    runner=ThermalEvolution(model,packet,engine)
    # Validate shapes before broadcasting domain bounds.
    if initial.ndim!=2 or initial.shape[0]!=model.operator.n or inputs.ndim!=3 or inputs.shape[1:]!=(len(model.B[0]),initial.shape[1]) or inputs.shape[0]!=len(packet['input_radii']):
        raise ValueError("Batch dimensions or horizon differ from certified domain")
    if not np.isfinite(initial).all() or not np.isfinite(inputs).all():
        raise ValueError("Finite scenario values required")
    from fractions import Fraction
    # Interpret binary64 values as exact rationals for box membership, so a
    # rounded decimal endpoint cannot accidentally admit an outside value.
    for i,(c,r) in enumerate(zip(packet['initial_center'],packet['initial_radius'])):
        low,high=q(c)-q(r),q(c)+q(r)
        if any(not low<=Fraction(float(value))<=high for value in initial[i]):
            raise ValueError("Scenario lies outside the certified initial box")
    for t,(centers,radii) in enumerate(zip(packet['input_centers'],packet['input_radii'])):
        for k,(c,r) in enumerate(zip(centers,radii)):
            low,high=q(c)-q(r),q(c)+q(r)
            if any(not low<=Fraction(float(value))<=high for value in inputs[t,k]):
                raise ValueError("Scenario lies outside a certified input box")
    return runner.simulate(initial,inputs)
