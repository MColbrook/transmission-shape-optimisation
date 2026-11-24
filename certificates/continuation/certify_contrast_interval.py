"""Exact interval certificate for a fixed smooth competitor on an n interval.

The polynomial trial bounds are inherited from the separately verified
near-optimal-domain certificate. The radial disk zero-free contour is checked
using directed fixed-point intervals and a convergent entire power series.
No floating point comparison is used in any certification assertion.
"""
from fractions import Fraction as F
from math import factorial, comb
from pathlib import Path
from functools import lru_cache
import sys, json, time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'shape_optimization'))
from calibration_interval import I,Q

NLO=F(259,50)
NHI=F(5281,1000)
NC=F(263447,50000)
K=32

class C:
    __slots__=('r','i')
    def __init__(self,r=0,i=0):
        if isinstance(r,C): self.r,self.i=r.r,r.i
        else: self.r,self.i=I(r),I(i)
    def __add__(self,b):
        b=C(b);return C(self.r+b.r,self.i+b.i)
    __radd__=__add__
    def __neg__(self):return C(-self.r,-self.i)
    def __sub__(self,b):return self+-C(b)
    def __rsub__(self,b):return C(b)+-self
    def __mul__(self,b):
        b=C(b);return C(self.r*b.r-self.i*b.i,self.r*b.i+self.i*b.r)
    __rmul__=__mul__
    def __truediv__(self,b):
        if isinstance(b,C):return self*C(b.r,-b.i)/(b.r*b.r+b.i*b.i)
        return C(self.r/b,self.i/b)
    def expand(self,e):return C(self.r.expand(e),self.i.expand(e))
    def nonzero(self):return self.r.lo>0 or self.r.hi<0 or self.i.lo>0 or self.i.hi<0
    def absupper(self):return max(self.r.absupper(),self.i.absupper())
    def approx(self):return [self.r.approx(),self.i.approx()]

def radius(n):return F(519631,100000)-F(5,4)*(n-NC)+(n-NC)**2/2
def radius_derivative(n):return n-NC-F(5,4)

def realpoly(p,x):
    out=I(0)
    for c in reversed(p):out=out*x+c
    return out

COEFF=[[F((-1)**k*comb(k,j)**2,4**k*(k+1)*factorial(k)**2)
        for j in range(k+1)] for k in range(K+1)]

@lru_cache(maxsize=20000)
def coefficients(lo,hi,derivative=0):
    n=I(lo,hi)
    polys=COEFF if not derivative else [[j*p[j] for j in range(1,len(p))] for p in COEFF]
    return tuple(realpoly(p,n) if p else I(0) for p in polys)

def horner(p,z):
    out=C(0)
    for c in reversed(p):out=out*z+c
    return out

def tail(order):
    # For 1<=n<=6 and |z|<=8, coefficients are majorized by 28^k/(k!)^2.
    # Either first partial derivative is majorized by k times this sequence.
    first=F((K+1)**order*28**(K+1),factorial(K+1)**2)
    ratio=F(28,(K+2)**2)*F(K+2,K+1)**order
    assert ratio<1
    return first/(1-ratio)

TAIL=tail(0);DTAIL=tail(1)

def circle(s):return C((1-s*s)/(1+s*s),2*s/(1+s*s))
def circle_box(lo,hi):
    # Both monotonicities are exact on [0,1].
    return C(I((1-hi*hi)/(1+hi*hi),(1-lo*lo)/(1+lo*lo)),
             I(2*lo/(1+lo*lo),2*hi/(1+hi*hi)))

def contour_cell(nlo,nhi,slo,shi):
    nm=(nlo+nhi)/2;sm=(slo+shi)/2
    dn=(nhi-nlo)/2;ds=(shi-slo)/2
    ni=I(nlo,nhi);si=I(slo,shi)
    rr=I(radius(nhi),radius(nlo))
    cc=circle_box(slo,shi);zz=cc*rr
    center=horner(coefficients(nm,nm),circle(sm)*radius(nm)).expand(TAIL)
    coeff=coefficients(nlo,nhi)
    dz=horner(tuple(k*coeff[k] for k in range(1,K+1)),zz).expand(DTAIL)
    dnpartial=horner(coefficients(nlo,nhi,1),zz).expand(DTAIL)
    compn=dnpartial+dz*cc*radius_derivative(ni)
    cs=C(-4*si,2*(1-si*si))/((1+si*si)**2)
    comps=dz*rr*cs
    variation_n=compn*I(-dn,dn)
    variation_s=comps*I(-ds,ds)
    result=center+variation_n+variation_s
    if result.nonzero():return True,result,None
    # Splitting choice affects only efficiency, never the verified bounds.
    split_n=variation_n.absupper()>=variation_s.absupper()
    return False,result,split_n

def padd(a,b):
    out=[F(0)]*max(len(a),len(b))
    for j,c in enumerate(a):out[j]+=c
    for j,c in enumerate(b):out[j]+=c
    return out
def pmul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for j,c in enumerate(a):
        for k,d in enumerate(b):out[j+k]+=c*d
    return out
def scale(a,s):return [s*c for c in a]
def affine_poly(p,lo,hi):
    # p(lo+(hi-lo)t), ordinary power coefficients in t.
    out=[F(0)]*len(p)
    for j,c in enumerate(p):
        for k in range(j+1):out[k]+=c*comb(j,k)*lo**(j-k)*(hi-lo)**k
    return out
def bernstein_bounds(p,lo,hi):
    a=affine_poly(p,lo,hi);degree=len(a)-1
    b=[sum(a[k]*F(comb(j,k),comb(degree,k)) for k in range(j+1)) for j in range(degree+1)]
    return min(b),max(b)

def real_trial_certificate():
    # Compact exact moment bounds already established in the parent certificate.
    Ahi=F(6923409,10**6);Blo=F(2941129141,10**10);Chi=F(7471409,400000000)
    t=[radius(F(0)),radius_derivative(F(0)),F(1,2)]
    p=padd([Ahi],padd(scale(pmul([1,1],t),-Blo),scale(pmul([0,1],pmul(t,t)),Chi)))
    todo=[(NLO,NHI)];parts=[]
    while todo:
        lo,hi=todo.pop();lower,upper=bernstein_bounds(p,lo,hi)
        if upper<0:parts.append([str(lo),str(hi),str(upper)])
        else:
            assert hi-lo>F(1,10**10)
            mid=(lo+hi)/2;todo.extend([(lo,mid),(mid,hi)])
    return parts

def normalized_angular(m,n,z):
    # D_m divided by its nonzero negative leading factor. D_m has no zero
    # on (0,z] if this is positive and D_m' changes sign at most once.
    # 2(m+1) int_0^1 r [J_m(sqrt(z)r)/(sqrt(z)/2)^m]
    #                     [J_m(sqrt(nz)r)/(sqrt(nz)/2)^m] dr.
    # Multiplying by (m!)^2 makes its value at z=0 equal to one.
    coeff=[]
    for k in range(K+1):
        coeff.append(sum(F((-1)**k*factorial(m)**2*(m+1),
                          4**k*(k+m+1)*factorial(j)*factorial(j+m)*factorial(k-j)*factorial(k-j+m))*n**j
                         for j in range(k+1)))
    out=I(0)
    for c in reversed(coeff):out=out*z+c
    # These coefficients in modulus are no larger than the m=0 ones.
    return out.expand(TAIL)

def angular_certificate():
    todo=[(NLO,NHI)];parts=[]
    while todo:
        lo,hi=todo.pop();ni=I(lo,hi);ti=I(radius(hi),radius(lo))
        vals=[normalized_angular(m,ni,ti) for m in (1,2)]
        if all(v.lo>0 for v in vals):parts.append([str(lo),str(hi),*[str(F(v.lo,Q)) for v in vals]])
        else:
            assert hi-lo>F(1,10**10)
            mid=(lo+hi)/2;todo.extend([(lo,mid),(mid,hi)])
    return parts

def bessel_zero_certificate():
    # J0(sqrt(z)) brackets the first two zeros. Root indexing uses the
    # elementary separation argument already given in the manuscript.
    def j0partial(z,degree):
        term=F(1);ans=term
        for k in range(1,degree+1):term*=-z/F(4*k*k);ans+=term
        return ans
    alo,ahi=F(5783,1000),F(723,125)
    blo,bhi=F(30471,1000),F(3809,125)
    assert j0partial(alo,21)>0 and j0partial(ahi,20)<0
    assert j0partial(blo,20)<0 and j0partial(bhi,21)>0
    assert NLO<blo/ahi and bhi/alo<NHI
    assert radius(NLO)<alo
    assert NHI*radius(NLO)<29<blo
    # Certify J3(sqrt(z))>0 on 0<z<=31. Its normalized series has value1
    # at zero; the omitted series is bounded by the larger TAIL above.
    p=[F((-1)**k*factorial(3),4**k*factorial(k)*factorial(k+3)) for k in range(K+1)]
    todo=[(F(0),F(31))];parts=[]
    while todo:
        lo,hi=todo.pop();lower,upper=bernstein_bounds(p,lo,hi)
        if lower>TAIL:parts.append([str(lo),str(hi),str(lower-TAIL)])
        else:
            assert hi-lo>F(1,10**8)
            mid=(lo+hi)/2;todo.extend([(lo,mid),(mid,hi)])
    # On this real series the term majorant (31/4)^k/(k!)^2 is smaller
    # than the 28^k/(k!)^2 majorant defining TAIL.
    return {'alpha':[str(alo),str(ahi)],'zeta':[str(blo),str(bhi)],
            'J3_positive_through':'31','J3_partition':parts}

def main():
    start=time.time()
    assert 1<NLO<NHI<6 and 0<radius(NHI)<radius(NLO)<F(16,3)
    # Each complex rectangle zz lies in |zz|<sqrt(2)*16/3<8.
    assert 2*F(16,3)**2<64
    trial=real_trial_certificate();print('Real trial:',len(trial),'cells',flush=True)
    angular=angular_certificate();print('Angular determinants:',len(angular),'cells',flush=True)
    bessel=bessel_zero_certificate()
    pointwise_comparison_n=F(131,25)
    pointwise_comparison_t=radius(pointwise_comparison_n)
    assert pointwise_comparison_n*pointwise_comparison_t**2<F(719,5)
    todo=[(NLO,NHI,F(0),F(1),0)];parts=[];nodes=0;maxdepth=0
    while todo:
        lo,hi,sl,sh,depth=todo.pop();nodes+=1
        good,val,split_n=contour_cell(lo,hi,sl,sh)
        if good:
            parts.append([str(lo),str(hi),str(sl),str(sh)])
            maxdepth=max(maxdepth,depth)
        else:
            if nodes%2000==0:print('nodes',nodes,'pending',len(todo),'accepted',len(parts),'depth',depth,flush=True)
            assert depth<65,(lo,hi,sl,sh,val.approx())
            if split_n:
                mid=(lo+hi)/2;todo.extend([(lo,mid,sl,sh,depth+1),(mid,hi,sl,sh,depth+1)])
            else:
                mid=(sl+sh)/2;todo.extend([(lo,hi,sl,mid,depth+1),(lo,hi,mid,sh,depth+1)])
    result={'status':'exact_assertions_passed','n_interval':[str(NLO),str(NHI)],
            'comparison_curve':'519631/100000 - (5/4)(n-263447/50000) + (n-263447/50000)^2/2',
            'series_degree':K,'series_tail':str(TAIL),'derivative_tail':str(DTAIL),
            'bessel_zero_certificate':bessel,
            'pointwise_comparison_n':str(pointwise_comparison_n),
            'pointwise_comparison_t':str(pointwise_comparison_t),
            'real_trial_partition':trial,'angular_partition':angular,
            'contour_cells':len(parts),'contour_nodes':nodes,'max_depth':maxdepth,
            'contour_partition':parts,'seconds':time.time()-start}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if 'partition' not in k},indent=2),flush=True)

if __name__=='__main__':main()
