"""Reusable numerical execution of the actual validated encoder/dynamics."""
import copy
import numpy as np
from .rational import q
from .uniform import SparseLinearModel,_domain
from .validated import replay_validated
from .prepared import PreparedSwitching


class PreparedValidated:
    def __init__(self,packet,engine='reduced',allow_unaccepted=False):
        if engine not in ('full','reduced'):raise ValueError('Full or reduced engine required')
        self.packet=copy.deepcopy(packet)
        if self.packet.get('schema')!='dynacomp.validated-box/1':raise ValueError('Unknown validated plan')
        if engine=='reduced':
            replay_validated(self.packet)
            if self.packet['status']!='tolerance_met' and not allow_unaccepted:raise ValueError('Reduced plan has not met tolerance')
        model=SparseLinearModel.from_spec(self.packet['model'])
        _domain(model,self.packet['initial_center'],self.packet['initial_radius'],self.packet['input_radii'],self.packet['input_centers'])
        self.n,self.m,self.p=model.operator.n,len(model.B[0]),len(model.C);self.horizon=len(self.packet['input_radii'])
        self.encoder=None
        if engine=='full':
            # Sampled external matrices are dense; dense BLAS is the baseline.
            self.a=np.zeros((self.n,self.n))
            for i,j,v in model.operator.entries:self.a[i,j]=float(v)
            self.b=np.asarray(model.B,dtype=float);self.c=np.asarray(model.C,dtype=float)
        else:
            def array(name):return np.asarray([[float(q(v)) for v in row] for row in self.packet[name]])
            self.encoder=array('E');self.a=array('G');self.b=array('H');self.c=array('D')
        values=(self.a,self.b,self.c) if self.encoder is None else (self.a,self.b,self.c,self.encoder)
        if any(not np.isfinite(v).all() for v in values):raise OverflowError('Nonfinite numerical matrices')
        pairs=[PreparedSwitching._float_endpoints(q(c)-q(r),q(c)+q(r)) for c,r in zip(self.packet['initial_center'],self.packet['initial_radius'])]
        self.initial_lower=np.asarray([lo for lo,hi in pairs])[:,None];self.initial_upper=np.asarray([hi for lo,hi in pairs])[:,None]
        pairs=[[PreparedSwitching._float_endpoints(q(c)-q(r),q(c)+q(r)) for c,r in zip(center,radius)] for center,radius in zip(self.packet['input_centers'],self.packet['input_radii'])]
        self.input_lower=np.asarray([[lo for lo,hi in row] for row in pairs]).reshape(self.horizon,self.m,1)
        self.input_upper=np.asarray([[hi for lo,hi in row] for row in pairs]).reshape(self.horizon,self.m,1)
        self.cache_stats={'packet_replays':int(engine=='reduced'),'batches':0,'initial_fields_encoded':0}

    def simulate(self,initial,inputs):
        initial,inputs=np.asarray(initial,dtype=float),np.asarray(inputs,dtype=float)
        if initial.ndim==1:initial=initial[:,None]
        if inputs.ndim!=3 or inputs.shape[:2]!=(self.horizon,self.m) or inputs.shape[2]<1:raise ValueError('Wrong input dimensions or certified horizon')
        count=inputs.shape[2]
        if initial.ndim!=2 or initial.shape[0]!=self.n or initial.shape[1] not in (1,count):raise ValueError('Wrong initial dimensions')
        if not np.isfinite(initial).all() or not np.isfinite(inputs).all():raise ValueError('Finite batch required')
        if (initial<self.initial_lower).any() or (initial>self.initial_upper).any() or (inputs<self.input_lower).any() or (inputs>self.input_upper).any():raise ValueError('Batch outside certified boxes')
        state=initial.copy() if self.encoder is None else self.encoder@initial
        if state.shape[1]==1 and count!=1:state=np.repeat(state,count,axis=1)
        output=np.empty((self.horizon+1,self.p,count));output[0]=self.c@state
        for t,u in enumerate(inputs):state=self.a@state+self.b@u;output[t+1]=self.c@state
        if not np.isfinite(output).all():raise OverflowError('Nonfinite numerical output')
        self.cache_stats['batches']+=1
        if self.encoder is not None:self.cache_stats['initial_fields_encoded']+=initial.shape[1]
        return output
