"""Common sparse observable closure, including time-varying affine families."""
from .uniform import SparseLinearModel
from .krylov import krylov
from .rational import matrix,rref,pivot_coordinates,mul,encoded,row_basis


def compile_sparse_family(modes,C,max_rank=128):
    if not isinstance(modes,dict) or not modes or any(not isinstance(name,str) or not name for name in modes):
        raise ValueError("Nonempty named mode dictionary required")
    models={name:SparseLinearModel.from_spec({**spec,"C":C}) for name,spec in modes.items()}
    first=next(iter(models.values()))
    n,m=first.operator.n,len(first.B[0])
    if any(model.operator.n!=n or len(model.B[0])!=m for model in models.values()) or type(max_rank) is not int or not 1<=max_rank<=n:
        raise ValueError("Matching mode/input dimensions and rank budget in [1,n] required")
    rows,labels=krylov([(i,row) for i,row in enumerate(first.C)],
                      {name:(lambda v,op=model.operator:op.pullback((v,))[0]) for name,model in models.items()},limit=max_rank)
    if not rows:
        raise ValueError("Nonzero target required")
    e,pivots=rref(rows)
    return {"schema":"dynacomp.sparse-family/1","modes":{name:model.spec() for name,model in models.items()},
            "C":encoded(first.C),"E":encoded(e),"D":encoded([pivot_coordinates(row,e,pivots) for row in first.C]),
            "G":{name:encoded([pivot_coordinates(row,e,pivots) for row in model.operator.pullback(e)]) for name,model in models.items()},
            "H":{name:encoded(mul(e,model.B)) for name,model in models.items()},
            "experiments":[{"output":i,"word":list(word)} for i,word in labels],
            "rank":len(e),"max_rank":max_rank,
            "coverage":"all initial states and input words under arbitrary time-varying linear combinations of the supplied mode/input pairs; output preservation does not establish stability"}


def replay_sparse_family(packet):
    if packet.get("schema")!="dynacomp.sparse-family/1":
        raise ValueError("Unknown sparse family schema")
    modes=packet["modes"]
    if not isinstance(modes,dict) or not modes or set(packet["G"])!=set(modes) or set(packet["H"])!=set(modes):
        raise ValueError("Mode counts disagree")
    c,e,d=map(matrix,(packet[k] for k in ("C","E","D")))
    rank=packet["rank"]
    n=len(e[0])
    if type(rank) is not int or type(packet["max_rank"]) is not int or not 1<=rank<=packet["max_rank"]<=n or len(e)!=rank or len(row_basis(e))!=rank or mul(d,e)!=c:
        raise ValueError("Invalid rank, encoder or target identity")
    channels=None
    for name,spec in modes.items():
        model=SparseLinearModel.from_spec(spec)
        if model.operator.n!=n or model.C!=c or (channels is not None and len(model.B[0])!=channels):
            raise ValueError("Mode dimensions or targets disagree")
        channels=len(model.B[0])
        if model.operator.pullback(e)!=mul(matrix(packet["G"][name]),e) or mul(e,model.B)!=matrix(packet["H"][name]):
            raise ValueError("Sparse family transition or input identity fails")
    return True
