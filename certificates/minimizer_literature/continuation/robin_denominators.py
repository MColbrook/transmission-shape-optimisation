"""Exact global bounds for E(x)=J0(sqrt(x))+sqrt(x)J1(sqrt(x))."""
from fractions import Fraction as F
from math import factorial,comb
from pathlib import Path
import json
from robin_roots import AL,BL,YL,DEN,audit_endpoints

K=32
P=[F((1-2*k)*(-1)**k,4**k*factorial(k)**2) for k in range(K+1)]
TAIL=F(2*(K+1)*9**(K+1),factorial(K+1)**2)/(1-F(9,(K+1)*(K+2)))

def bernstein(p,lo,hi):
    n=len(p)-1
    a=[sum(p[j]*comb(j,k)*lo**(j-k)*(hi-lo)**k for j in range(k,n+1)) for k in range(n+1)]
    b=[sum(a[k]*F(comb(j,k),comb(n,k)) for k in range(j+1)) for j in range(n+1)]
    return min(b)-TAIL,max(b)+TAIL

def certificate():
    audit_endpoints();out=[]
    for branch,interval in [(0,(F(0),F(AL+1,DEN))),(1,(F(BL,DEN),F(YL+1,DEN)))]:
        todo=[interval]
        while todo:
            lo,hi=todo.pop();lower,upper=bernstein(P,lo,hi)
            if (branch==0 and lower>F(1,2)) or (branch==1 and upper<F(-1,4)):
                out.append(dict(branch=branch,lo=str(lo),hi=str(hi),lower=str(lower),upper=str(upper)))
            else:
                assert hi-lo>F(1,10**8)
                mid=(lo+hi)/2;todo.extend([(lo,mid),(mid,hi)])
    return out

if __name__=='__main__':
    r=certificate();out={'status':'exact_E_bounds_passed','cells':len(r),'tail':str(TAIL),'partition':r}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print('certified cells',len(r));print('E(lambda)>1/2; E(mu)<-1/4')
