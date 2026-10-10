"""Optional Lyapunov-coordinate discovery with exact reduced contraction.

The full metric and basis are floating proposals. Only supplied reduced
matrices passing I-G^T G > 0 exactly claim reduced stability. Output accuracy
needs separate validated_box acceptance. No full-system stability proof here.
"""
from fractions import Fraction
from .uniform import as_sparse
from .rational import encoded
from .validated import euclidean_contraction


class StableProposalContext:
    """Prepare one numerical metric and extend a nested proposal basis."""
    def __init__(self,model,horizon):
        import numpy as np
        from scipy.linalg import solve_discrete_lyapunov,cholesky,solve_triangular
        self.model=as_sparse(model);self.horizon=horizon;n=self.model.operator.n
        if type(horizon) is not int or horizon<0:raise ValueError('Nonnegative integer horizon required')
        a=np.zeros((n,n))
        for i,j,v in self.model.operator.entries:a[i,j]=float(v)
        b=np.asarray(self.model.B,dtype=float);c=np.asarray(self.model.C,dtype=float)
        if not np.isfinite(a).all() or not np.isfinite(b).all() or not np.isfinite(c).all():raise ValueError('Source not representable in finite binary64')
        p=solve_discrete_lyapunov(a.T,np.eye(n));l=cholesky((p+p.T)/2,lower=False)
        li=solve_triangular(l,np.eye(n));self.at=l@a@li;self.bt=l@b;self.ct=c@li
        _,singular,vh=np.linalg.svd(self.ct,full_matrices=False)
        count=int(np.sum(singular>max(self.ct.shape)*np.finfo(float).eps*(singular[0] if len(singular) else 0)))
        if not count:raise ValueError('Nonzero numerical target required')
        self.basis=vh[:count].T;self.target_rank=count;self.pulls=[c]
        for _ in range(horizon):self.pulls.append(self.pulls[-1]@a)
        self.a,self.b,self.c,self.l,self.li=a,b,c,l,li
        self.condition=float(np.linalg.cond(p));self.metric_preparations=1

    def propose(self,rank):
        import numpy as np
        if type(rank) is not int or not self.target_rank<=rank<=self.model.operator.n:raise ValueError('Adequate rank budget in [1,n] required')
        a,b,c,l,li=self.a,self.b,self.c,self.l,self.li
        basis=self.basis
        while basis.shape[1]<rank:
            g=basis.T@self.at@basis;h=basis.T@self.bt;d=self.ct@basis;e=basis.T@l;dg=d
            residuals=[];forcing=[]
            for ca in self.pulls:
                residuals.append(ca-dg@e);forcing.append(ca@b-dg@h);dg=dg@g
            scores=[np.sum(abs(v),axis=1)+sum((np.sum(abs(f),axis=1) for f in forcing[:t]),start=np.zeros(len(c))) for t,v in enumerate(residuals)]
            t,k=max(((t,k) for t in range(len(scores)) for k in range(len(c))),key=lambda pair:scores[pair[0]][pair[1]])
            row=residuals[t][k]@li
            v=row-basis@(basis.T@row);v-=basis@(basis.T@v)
            threshold=1e-12*max(1,float(np.linalg.norm(row)))
            if np.linalg.norm(v)<=threshold:
                found=False
                for ca in self.pulls:
                    for row in ca@li:
                        v=row-basis@(basis.T@row);v-=basis@(basis.T@v)
                        if np.linalg.norm(v)>1e-12*max(1,float(np.linalg.norm(row))):found=True;break
                    if found:break
                if not found:break
            basis=np.column_stack([basis,v/np.linalg.norm(v)])
        self.basis=basis;basis=basis[:,:rank]
        g=basis.T@self.at@basis;h=basis.T@self.bt;d=self.ct@basis;e=basis.T@l
        if not all(np.isfinite(x).all() for x in (e,g,h,d)):raise OverflowError('Nonfinite proposal')
        def bind(x):return tuple(tuple(Fraction(float(v)) for v in row) for row in x)
        e,g,h,d=map(bind,(e,g,h,d));stable=euclidean_contraction(g)
        return {'E':encoded(e),'G':encoded(g),'H':encoded(h),'D':encoded(d),'rank':len(g),
                'reduced_contraction_certified':stable,'proposal':'Lyapunov-coordinate output-visible enrichment; binary64 matrices bound as exact dyadics',
                'full_metric_condition_estimate':self.condition,
                'claim':'reduced stability only if exact contraction accepted; accuracy requires independent box certificate'}


def stable_observable_proposal(model,horizon,rank,context=None):
    model=as_sparse(model)
    if type(rank) is not int or not 1<=rank<=model.operator.n:raise ValueError('Rank in [1,n] required')
    if context is None:context=StableProposalContext(model,horizon)
    if context.model!=model or context.horizon!=horizon:raise ValueError('Proposal context source or horizon differs')
    return context.propose(rank)
