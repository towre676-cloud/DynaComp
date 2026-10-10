"""Optional vectorized batch execution shared by the application and benchmark."""
import numpy as np
from scipy.sparse import csr_matrix
from .rational import q


class NumpyEvolution:
    """Prepare an already validated Network/plan pair; columns are scenarios.

    Public JSON entry points validate the exact plan before constructing this
    executor. Numeric traces use ordinary binary64 arithmetic, not intervals.
    """
    def __init__(self, net, packet, engine):
        self.n, self.inputs, self.outputs = net.op.n, len(net.channels), 2*len(net.names)
        self.engine = engine
        if engine == "full":
            rows,cols,values=zip(*net.op.entries) if net.op.entries else ([],[],[])
            self.a=csr_matrix((np.asarray(values,dtype=float),(rows,cols)),shape=(net.op.n,net.op.n))
            self.b=csr_matrix(np.asarray(net.B,dtype=float))
            self.c=csr_matrix(([1.]*net.op.n+list(map(float,net.exits)),
                (list(net.regions)+[len(net.names)+r for r in net.regions],list(range(net.op.n))*2)),shape=(self.outputs,net.op.n))
            self.e=None
        elif engine == "reduced":
            self.e=csr_matrix((np.ones(net.op.n),(packet["classes"],np.arange(net.op.n))),shape=(packet["rank"],net.op.n))
            self.a,self.b,self.c=(np.asarray([[float(q(v)) for v in row] for row in packet[k]]) for k in ("G","H","D"))
            if packet["rank"]>=32:
                self.a,self.b,self.c=map(csr_matrix,(self.a,self.b,self.c))
        else:
            raise ValueError("Engine must be reduced or full")

    def simulate(self, initial, inputs):
        initial=np.asarray(initial,dtype=float)
        inputs=np.asarray(inputs,dtype=float)
        if initial.ndim != 2 or initial.shape[0] != self.n or initial.shape[1] < 1 or inputs.ndim != 3 or inputs.shape[1:] != (self.inputs,initial.shape[1]):
            raise ValueError("Expected initial (nodes, scenarios), input (steps, channels, scenarios)")
        if not np.isfinite(initial).all() or not np.isfinite(inputs).all() or (initial<0).any() or (inputs<0).any():
            raise ValueError("Finite nonnegative numerical stocks and inputs required")
        state=initial.copy() if self.e is None else self.e@initial
        output=np.empty((len(inputs)+1,self.outputs,initial.shape[1]))
        output[0]=self.c@state
        for t,u in enumerate(inputs):
            state=self.a@state+self.b@u
            output[t+1]=self.c@state
        if not np.isfinite(output).all():
            raise OverflowError("Nonfinite numerical output; rescale quantities or use exact arithmetic")
        return output
