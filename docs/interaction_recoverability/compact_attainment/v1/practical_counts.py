"""Evaluate public work/precision formulas at two fixed counts; not a benchmark."""
import json
import math
from pathlib import Path
import mpmath as mp
from estimator import schedule


def bell(r, m, x):
    if r == m == 0:
        return mp.mpf(1)
    if r == 0 or m == 0:
        return mp.mpf(0)
    return mp.fsum(math.comb(r-1,j-1)*x[j]*bell(r-j,m-1,x)
                   for j in range(1,r-m+2))


def compose(r, outer, inner):
    return mp.fsum(outer[m]*bell(r,m,inner) for m in range(1,r+1))


def count(N):
    s=schedule(N); K,M,R=s.depth,s.lattice,s.radius
    lam=mp.mpf(271)/290
    H=lambda B:{1:mp.mpf('1.1')+mp.mpf('.2')*B,2:mp.mpf('.2'),
                3:mp.mpf('.3'),4:mp.mpf('1.5'),5:mp.mpf('10.5')}
    a={1:mp.mpf(1),2:mp.mpf('.2'),3:mp.mpf('.42'),4:mp.mpf('2.22'),5:mp.mpf('17.328')}
    d={0:1,1:1,2:1,3:2,4:3}
    b={r:compose(r,a,H(mp.mpf(7)/4)) for r in range(1,5)}
    D={1:mp.mpf(1)}
    for r in range(2,5):
        D[r]=mp.fsum(b[m]*bell(r,m,D) for m in range(2,r+1))/(1-lam)
    f={r:compose(r,b,D) for r in range(1,5)}
    p={r:compose(r,{m:d[m-1] for m in range(1,5)},f) for r in range(1,5)}
    G={0:mp.mpf(1)}
    G.update({r:compose(r,d,a) for r in range(1,5)})
    rho={r:mp.fsum(math.comb(r,l)*G[l]*a[r-l+1] for l in range(r+1)) for r in range(5)}
    def A(B):
        hs=H(B); C={0:rho[0]}
        C.update({r:compose(r,rho,hs) for r in range(1,5)})
        return {r:4*mp.fsum(math.comb(r,l)*hs[r-l+1]*C[l] for l in range(r+1)) for r in range(5)}
    AB,AR=A(R+M+4),A(R)
    L={r:4*d[r]+AB[r] for r in range(5)}
    P={r:(2*M+1)*L[r] for r in range(5)}
    W={r:M*L[r] for r in range(5)}
    U={r:AR[r]+W[r] for r in range(5)}
    C4c=36*K*P[4]+18*K*mp.fsum(math.comb(4,l)*P[4-l]*p[l] for l in range(1,5))
    C4o=U[4]+W[4]+mp.fsum(math.comb(4,l)*(U[4-l]+W[4-l])*d[l-1] for l in range(1,5))
    C4=C4c+C4o; S=1+2*R
    mesh=min(mp.mpf(1),(mp.mpf(4320)/(8*N**4*C4*S))**mp.mpf('.25'))
    knots=mp.ceil(36*(N//5)/(1-lam)+36*K+2*(N//5)+4)
    panel_bound=int(knots+mp.ceil(S/mesh)+2)
    # Explicit conservative coefficient/data sensitivity constants in algorithm.md.
    LN=1+6900*N**3*(K+1)
    AT=(mp.mpf(LN)**s.iterations-1)/(LN-1)
    Lc=2+2*H(R+M+4)[1]*rho[1]
    Jc=2*N**2*(2*R*(2*H(R)[1]*rho[1]+2*M*Lc)+36*K*(2*M+1)*Lc)
    JY=2*R*(U[0]+W[0])+36*K*P[0]
    dc=1/(16*N*(1+Jc)); Cdata=10350*(K+1)
    tau=min(dc/(4*AT*Cdata),1/(32*N*(1+JY)))
    return {'N':N,'C4_bound_diagnostic':str(C4),'gauss_mesh':str(mesh),
            'public_panel_overcount':panel_bound,'log2_iteration_amplification':str(mp.log(AT,2)),
            'sufficient_input_fraction_bits_for_allocated_tau':int(mp.ceil(-mp.log(tau,2))),
            'interpretation':'Formula counts only, no panels, samples or complete estimator evaluated.'}


if __name__=='__main__':
    mp.mp.dps=40
    path=Path(__file__).with_name('practical_counts.json')
    if path.exists():
        raise FileExistsError('Preserve the existing count diagnostic')
    result={'status':'passed','scope':'Fixed formula evaluation, not certified interval rounding of constants.',
            'counts':[count(N) for N in (1000,1000000)]}
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
