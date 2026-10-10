"""Sharp finite-horizon box errors and general-direction enrichment.

Every output error is a known linear function of the initial state and the
independent input boxes. Its exact range follows from absolute coefficient
sums, with a retained attaining witness. Sparse target pullbacks are cached;
no full-state matrix powers or full-state reconstruction envelope is needed.
"""
from dataclasses import dataclass
from .rational import q, matrix, mul, sub, mv, addv, transpose, inverse, identity, encoded
from .sparse import SparseOperator
from .krylov import Echelon


@dataclass(frozen=True)
class SparseLinearModel:
    operator: SparseOperator
    B: tuple
    C: tuple

    @classmethod
    def from_spec(cls, spec):
        op = SparseOperator.from_entries(spec["n"], spec["entries"])
        b,c = matrix(spec["B"]),matrix(spec["C"])
        if len(b)!=op.n or len(c[0])!=op.n:
            raise ValueError("Expected B[n,m] and C[p,n]")
        return cls(op,b,c)

    def spec(self):
        return {"n":self.operator.n,"entries":[[i,j,str(v)] for i,j,v in self.operator.entries],
                "B":encoded(self.B),"C":encoded(self.C)}


def as_sparse(model):
    if isinstance(model,SparseLinearModel):
        return model
    n=len(model.A)
    return SparseLinearModel(SparseOperator.from_entries(n,[(i,j,v) for i,row in enumerate(model.A) for j,v in enumerate(row) if v]),model.B,model.C)


def _domain(model,center,radius,input_radii,input_centers):
    center,radius=tuple(q(v) for v in center),tuple(q(v) for v in radius)
    if len(center)!=model.operator.n or len(radius)!=model.operator.n or any(v<0 for v in radius):
        raise ValueError("Initial box needs n coordinates and nonnegative radii")
    us=tuple(tuple(q(v) for v in row) for row in input_radii)
    m=len(model.B[0])
    if any(len(row)!=m or any(v<0 for v in row) for row in us):
        raise ValueError("Input radii must have m nonnegative coordinates")
    centers=tuple(tuple(q(v) for v in row) for row in input_centers) if input_centers is not None else tuple((q(0),)*m for _ in us)
    if len(centers)!=len(us) or any(len(row)!=m for row in centers):
        raise ValueError("One matching input center per time required")
    return center,radius,us,centers


def _pullbacks(model,horizon):
    rows=[model.C]
    for _ in range(horizon):
        rows.append(model.operator.pullback(rows[-1]))
    return rows


def _project(model,rows):
    e=matrix(rows)
    v=mul(transpose(e),inverse(mul(e,transpose(e))))
    g=mul(model.operator.pullback(e),v)
    return v,e,g,mul(e,model.B),mul(model.C,v)


def _kernels(model,e,g,h,d,pullbacks):
    dg=d
    state,forcing=[],[]
    for ca in pullbacks:
        state.append(sub(ca,mul(dg,e)))
        forcing.append(sub(mul(ca,model.B),mul(dg,h)))
        dg=mul(dg,g)
    return state,forcing


def _sign(x):
    return int(x>0)-int(x<0)


def _bound(model,candidate,domain,pullbacks):
    v,e,g,h,d=candidate
    center,radius,us,centers=domain
    state,forcing=_kernels(model,e,g,h,d,pullbacks)
    envelopes,central=[],[]
    for t,kernel in enumerate(state):
        mean=list(mv(kernel,center))
        uncertain=[sum(abs(a)*r for a,r in zip(row,radius)) for row in kernel]
        for j in range(t):
            lag=forcing[t-1-j]
            mean=list(addv(mean,mv(lag,centers[j])))
            uncertain=[a+sum(abs(v)*r for v,r in zip(row,us[j])) for a,row in zip(uncertain,lag)]
        central.append(mean)
        envelopes.append(tuple(abs(a)+b for a,b in zip(mean,uncertain)))
    t,output=max(((t,k) for t in range(len(envelopes)) for k in range(len(model.C))),key=lambda pair:envelopes[pair[0]][pair[1]])
    orientation=-1 if central[t][output]<0 else 1
    x0=tuple(c+orientation*_sign(a)*r for c,a,r in zip(center,state[t][output],radius))
    word=[]
    for j,(u,c) in enumerate(zip(us,centers)):
        coefficients=forcing[t-1-j][output] if j<t else (q(0),)*len(u)
        word.append(tuple(v+orientation*_sign(a)*r for v,a,r in zip(c,coefficients,u)))
    signed=central[t][output]+orientation*(envelopes[t][output]-abs(central[t][output]))
    cert={"schema":"dynacomp.sharp-box-bound/1","arithmetic":"rational",
          "coverage":"all initial states and independently chosen inputs in the declared componentwise boxes",
          "bound":[str(max(row)) for row in envelopes],"output_bounds":encoded(envelopes),
          "worst_bound":str(envelopes[t][output]),
          "attaining_witness":{"time":t,"output":output,"x0":[str(v) for v in x0],"inputs":encoded(word) if word else [],"signed_error":str(signed)},
          "claim":"sharp componentwise and worst-horizon output-error suprema for the fixed full and reduced models; excludes floating-point rounding"}
    return cert,state,forcing


def sharp_box_bound(model,V,E,G,H,initial_center,initial_radius,input_radii,input_centers=None):
    model=as_sparse(model)
    v,e,g,h=map(matrix,(V,E,G,H))
    n,r=model.operator.n,len(g)
    if len(v)!=n or len(v[0])!=r or len(e)!=r or len(e[0])!=n or len(g[0])!=r or len(h)!=r or len(h[0])!=len(model.B[0]) or mul(e,v)!=identity(r):
        raise ValueError("Expected V[n,r], E[r,n], G[r,r], H[r,m] and EV=I")
    domain=_domain(model,initial_center,initial_radius,input_radii,input_centers)
    return _bound(model,(v,e,g,h,mul(model.C,v)),domain,_pullbacks(model,len(domain[2])))[0]


def uniform_enrichment(model,initial_center,initial_radius,input_radii,tolerance,input_centers=None,max_rank=None,seed_rows=None,strategy="worst_kernel"):
    """Grow general rational directions with sharp box-uniform acceptance.

    Worst-kernel enrichment adds a missing output sensitivity at the worst
    time. If that row is dependent, it falls back to the first missing target
    pullback. Prefix strategy always chooses the latter. No monotonic error
    or globally optimal rank claim is made. Finite-horizon closure suffices.
    """
    model=as_sparse(model)
    n=model.operator.n
    domain=_domain(model,initial_center,initial_radius,input_radii,input_centers)
    tol=q(tolerance)
    budget=n if max_rank is None else max_rank
    if tol<0 or type(budget) is not int or not 1<=budget<=n or strategy not in ("worst_kernel","prefix"):
        raise ValueError("Nonnegative tolerance, rank budget in [1,n] and valid strategy required")
    span=Echelon()
    for row in model.C:
        span.add(row)
    if seed_rows is not None:
        for row in matrix(seed_rows):
            if len(row)!=n:
                raise ValueError("Seed directions need n coordinates")
            span.add(row)
    if not len(span) or len(span)>budget:
        raise ValueError("Nonzero targets and sufficient initial rank budget required")
    pullbacks=_pullbacks(model,len(domain[2]))
    history=[]
    while True:
        candidate=_project(model,[row for _,row in span.rows])
        cert,state,forcing=_bound(model,candidate,domain,pullbacks)
        history.append({"rank":len(span),"worst_bound":cert["worst_bound"],"worst_time":cert["attaining_witness"]["time"],"worst_output":cert["attaining_witness"]["output"]})
        if q(cert["worst_bound"])<=tol:
            status="tolerance_met"
            break
        if len(span)>=budget:
            status="rank_budget_exhausted"
            break
        witness=cert["attaining_witness"]
        t,k=witness["time"],witness["output"]
        chosen=None
        if strategy=="worst_kernel" and any(span.reduce(state[t][k])):
            chosen=state[t][k]
            label={"kind":"worst_output_kernel","time":t,"output":k}
        if chosen is None:
            for depth,ca in enumerate(pullbacks):
                for output,row in enumerate(ca):
                    if any(span.reduce(row)):
                        chosen=row
                        label={"kind":"target_pullback","time":depth,"output":output}
                        break
                if chosen is not None:
                    break
        if chosen is None:
            raise ArithmeticError("Finite-horizon target span should have zero kernels")
        history[-1]["added_direction"]=label
        assert span.add(chosen)
    v,e,g,h,d=candidate
    center,radius,us,centers=domain
    return {"schema":"dynacomp.uniform-enrichment/1","model":model.spec(),"status":status,"tolerance":str(tol),"max_rank":budget,"rank":len(e),"strategy":strategy,
            "initial_center":[str(x) for x in center],"initial_radius":[str(x) for x in radius],
            "input_radii":encoded(us) if us else [],"input_centers":encoded(centers) if centers else [],
            "V":encoded(v),"E":encoded(e),"G":encoded(g),"H":encoded(h),"D":encoded(d),
            "certificate":cert,"history":history,
            "claim":"general-direction projection with sharp fixed-model box-uniform output-error acceptance; enrichment is heuristic"}


def replay_uniform(packet):
    """Recompute projection and exact box certificate, without enrichment."""
    if packet.get("schema")!="dynacomp.uniform-enrichment/1":
        raise ValueError("Unknown uniform enrichment schema")
    model=SparseLinearModel.from_spec(packet["model"])
    v,e,g,h,d=map(matrix,(packet[k] for k in ("V","E","G","H","D")))
    rank=packet["rank"]
    if type(rank) is not int or type(packet["max_rank"]) is not int or not 1<=rank<=packet["max_rank"]<=model.operator.n or len(e)!=rank or len(e[0])!=model.operator.n:
        raise ValueError("Invalid dimensions or budget")
    if _project(model,e)!=(v,e,g,h,d):
        raise ValueError("Projection identities fail")
    certificate=sharp_box_bound(model,v,e,g,h,packet["initial_center"],packet["initial_radius"],packet["input_radii"],packet["input_centers"])
    if certificate!=packet["certificate"]:
        raise ValueError("Stored bound or attaining witness does not match")
    tol=q(packet["tolerance"])
    expected="tolerance_met" if q(certificate["worst_bound"])<=tol else "rank_budget_exhausted"
    if tol<0 or packet["status"]!=expected or (expected=="rank_budget_exhausted" and rank!=packet["max_rank"]):
        raise ValueError("Invalid acceptance status")
    # Also execute the worst-case witness as a separate acceptance check.
    witness=certificate["attaining_witness"]
    x=tuple(q(v) for v in witness["x0"])
    z=mv(e,x)
    for u in witness["inputs"][:witness["time"]]:
        x=addv(model.operator.step(x),mv(model.B,u))
        z=addv(mv(g,z),mv(h,u))
    actual=mv(model.C,x)[witness["output"]]-mv(d,z)[witness["output"]]
    if actual!=q(witness["signed_error"]) or abs(actual)!=q(certificate["worst_bound"]):
        raise ValueError("Attaining witness does not attain the bound")
    return True


def output_box_range(model,initial_center,initial_radius,input_radii,input_centers=None):
    """Sharp full-model output intervals at each time for independent boxes.

    These are ranges of the output itself, distinct from reduction error.
    Different output/time extrema may require different admissible inputs.
    """
    model=as_sparse(model)
    center,radius,us,centers=_domain(model,initial_center,initial_radius,input_radii,input_centers)
    ca=_pullbacks(model,len(us))
    cb=[mul(row,model.B) for row in ca]
    means,widths=[],[]
    for t,row in enumerate(ca):
        mean=mv(row,center)
        width=tuple(sum(abs(v)*r for v,r in zip(output,radius)) for output in row)
        for j in range(t):
            kernel=cb[t-1-j]
            mean=addv(mean,mv(kernel,centers[j]))
            width=addv(width,tuple(sum(abs(v)*r for v,r in zip(output,us[j])) for output in kernel))
        means.append(mean);widths.append(width)
    return {"schema":"dynacomp.output-box-range/1","center":encoded(means),"radius":encoded(widths),
            "lower":encoded(tuple(tuple(a-b for a,b in zip(mean,width)) for mean,width in zip(means,widths))),
            "upper":encoded(tuple(tuple(a+b for a,b in zip(mean,width)) for mean,width in zip(means,widths))),
            "coverage":"sharp individual time/output intervals for a fixed model and independently chosen state/input boxes; not a simultaneous trajectory enclosure with independently attainable endpoints"}
