"""Exact directed-interval enclosures of the two radial Robin branches.

For eta in [0,1], solve
  (1-eta) J0(sqrt(x)) - eta sqrt(x) J1(sqrt(x)) = 0.
The branches run alpha->0 and zeta->beta as eta runs 0->1.
Only Python integers and Fraction enter any certification decision.
"""
from fractions import Fraction as F
from pathlib import Path
from functools import lru_cache
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'shape_optimization'))
from calibration_interval import I,Q,j0,hj1

DEN=10**80
AL=578318596294678452117599575845580703507144180642368558708712371445606430488554437
BL=1468197064212389325721977776863010698968158698215497813417558917182975765058406195
YL=3047126234366208639907816317502275584842083370855398235487607831463637089182269021
ALPHA=I(F(AL,DEN),F(AL+1,DEN))
BETA=I(F(BL,DEN),F(BL+1,DEN))
ZETA=I(F(YL,DEN),F(YL+1,DEN))

def audit_endpoints():
    assert j0(I(F(AL,DEN))).lo>0 and j0(I(F(AL+1,DEN))).hi<0
    assert hj1(I(F(BL,DEN))).lo>0 and hj1(I(F(BL+1,DEN))).hi<0
    assert j0(I(F(YL,DEN))).hi<0 and j0(I(F(YL+1,DEN))).lo>0

def j1div(x,K=80):
    # sqrt(x)J1(sqrt(x))/x, continuously extended by 1/2 at zero.
    term=I(F(1,2));total=term
    for k in range(1,K+1):
        term=term*(-x)/(4*k*(k+1));total=total+term
    nxt=term*(-x)/(4*(K+1)*(K+2))
    err=I(F(nxt.absupper(),Q))/(1-x/(4*(K+2)*(K+3)))
    return total.expand(err)

def robin(x,eta):
    return (1-eta)*j0(x)-eta*hj1(x)

def robin_derivative(x,eta):
    return -(1-eta)*j1div(x)/2-eta*j0(x)/2

@lru_cache(maxsize=20000)
def root_at(eta,branch):
    eta=F(eta);assert 0<=eta<=1 and branch in (0,1)
    if eta==0:return ALPHA if branch==0 else ZETA
    if eta==1:return I(0) if branch==0 else BETA
    # The endpoints below enclose the whole branch, including the small
    # uncertainty in the fixed Dirichlet/Neumann Bessel roots.
    lo=F(0) if branch==0 else F(BL,DEN)
    hi=F(AL+1,DEN) if branch==0 else F(YL+1,DEN)
    sign=-1 if branch==0 else 1
    target=F(1,10**65)
    for iteration in range(300):
        if hi-lo<target:return I(lo,hi)
        mid=(lo+hi)/2
        fm=robin(I(mid),eta)
        # Ordinary monotone bisection keeps an exact root bracket.
        if sign*fm.lo>0 and sign*fm.hi>0:hi=mid
        elif sign*fm.lo<0 and sign*fm.hi<0:lo=mid
        else:
            # At this precision a midpoint whose sign is unresolved is
            # treated only by interval Newton, never as an exact root.
            pass
        derivative=robin_derivative(I(lo,hi),eta)
        if derivative.lo>0 or derivative.hi<0:
            newton=I(mid)-fm/derivative
            nl=max(lo,F(newton.lo,Q));nh=min(hi,F(newton.hi,Q))
            assert nl<=nh
            if nh-nl<hi-lo:lo,hi=nl,nh
    raise AssertionError(('root enclosure did not converge',eta,branch,lo,hi))

@lru_cache(maxsize=20000)
def root_boxes(eta_lo,eta_hi):
    """Return (lambda_interval, mu_interval), using strict monotonicity."""
    eta_lo=F(eta_lo);eta_hi=F(eta_hi)
    assert 0<=eta_lo<=eta_hi<=1
    lows=[root_at(eta_hi,j) for j in (0,1)]
    highs=[root_at(eta_lo,j) for j in (0,1)]
    return tuple(I(min(l.lo,h.lo),max(l.hi,h.hi),True) for l,h in zip(lows,highs))

if __name__=='__main__':
    import time
    start=time.time();audit_endpoints()
    for e in [F(0),F(1,10000),F(1,10),F(1,2),F(9,10),F(9999,10000),F(1)]:
        print(e,[r.approx() for r in root_boxes(e,e)],'seconds',time.time()-start,flush=True)
