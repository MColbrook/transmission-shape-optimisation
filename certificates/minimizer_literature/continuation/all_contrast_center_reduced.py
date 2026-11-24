"""Exact second-order Taylor-model certificate for the calibration center."""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
from pathlib import Path
import json,time
from robin_roots import I,Q,ALPHA,BETA,ZETA,root_boxes,audit_endpoints,hj1
from robin_denominators import certificate as denominator_certificate
K=30;SMAX=F(9,16)
TAIL=F(9**(K+1),factorial(K+1)**2)/(1-F(9,(K+2)**2))
COEFF=[[F((-1)**k,4**k*factorial(k)*factorial(k+m)) for k in range(K+1)] for m in range(5)]
@lru_cache(maxsize=50000)
def bn_cache(m,lo,hi):
    assert 0<=lo<=hi<=36*Q
    x=I(F(lo+hi,2*Q));ans=I(0)
    for c in reversed(COEFF[m]):ans=ans*x+c
    return ans.expand(TAIL).expand(F(hi-lo,8*Q*factorial(m+1)))
def bn(m,x):return bn_cache(m,x.lo,x.hi)
class J:
    __slots__=('v','e','s','ee','es','ss')
    def __init__(self,v=0,e=0,s=0,ee=0,es=0,ss=0):
        if isinstance(v,J):
            for a in self.__slots__:setattr(self,a,getattr(v,a))
        else:
            for a,b in zip(self.__slots__,(v,e,s,ee,es,ss)):setattr(self,a,I(b))
    def __add__(self,b):
        b=J(b);return J(*[getattr(self,a)+getattr(b,a) for a in self.__slots__])
    __radd__=__add__
    def __neg__(self):return J(*[-getattr(self,a) for a in self.__slots__])
    def __sub__(self,b):return self+-J(b)
    def __rsub__(self,b):return J(b)+-self
    def __mul__(self,b):
        b=J(b);a=self
        return J(a.v*b.v,a.e*b.v+a.v*b.e,a.s*b.v+a.v*b.s,
                 a.ee*b.v+2*a.e*b.e+a.v*b.ee,
                 a.es*b.v+a.e*b.s+a.s*b.e+a.v*b.es,
                 a.ss*b.v+2*a.s*b.s+a.v*b.ss)
    __rmul__=__mul__
    def inv(self):
        a=self;q=1/a.v
        return J(q,-a.e*q*q,-a.s*q*q,2*a.e*a.e*q**3-a.ee*q*q,
                 2*a.e*a.s*q**3-a.es*q*q,2*a.s*a.s*q**3-a.ss*q*q)
    def __truediv__(self,b):return self*J(b).inv()
    def __rtruediv__(self,b):return J(b)*self.inv()
    def __pow__(self,k):
        assert isinstance(k,int) and k>=0
        if k==0:return J(1)
        if k==1:return self
        f=k*self.v**(k-1);g=k*(k-1)*self.v**(k-2)
        return J(self.v**k,f*self.e,f*self.s,g*self.e*self.e+f*self.ee,
                 g*self.e*self.s+f*self.es,g*self.s*self.s+f*self.ss)
def BN(m,x):
    x=J(x);v=bn(m,x.v);d=-bn(m+1,x.v)/4;dd=bn(m+2,x.v)/16
    return J(v,d*x.e,d*x.s,dd*x.e*x.e+d*x.ee,dd*x.e*x.s+d*x.es,dd*x.s*x.s+d*x.ss)
@lru_cache(maxsize=20000)
def data(el,eh):
    ei=I(el,eh);eta=J(ei,1,0);roots=root_boxes(el,eh);xs=[];es=[]
    for branch,x in enumerate(roots):
        b0=bn(0,x);b1=bn(1,x);b2=bn(2,x)
        E=b0+x*b1/2
        E=I(max(E.lo,Q//2),E.hi,True) if branch==0 else I(E.lo,min(E.hi,-Q//4),True)
        assert E.lo<=E.hi
        Rx=-(1-ei)*b1/4-ei*b0/2
        Ex=b0/2-b1/4;Exx=-b1/8+b2/16
        Rxx=(1-ei)*b2/16+ei*b1/8
        if branch==0 and el>=F(1,2):
            rat=b1/(2*b0);den=ei*ei+ei*(1-ei)*rat
            xp=-2/den
            ratx=(b1*b1-b0*b2)/(8*b0*b0)
            denp=2*ei+(1-2*ei)*rat+ei*(1-ei)*ratx*xp
            xpp=2*denp/(den*den)
        elif branch==1 or eh<=F(1,2):
            den=(1-ei)**2+ei*ei*x
            xp=-2*x/den
            denp=-2*(1-ei)+2*ei*x+ei*ei*xp
            xpp=-2*xp/den+2*x*denp/(den*den)
        else:
            xp=E/Rx;xpp=(2*Ex*xp-Rxx*xp*xp)/Rx
        xs.append(J(x,xp,0,xpp,0,0));es.append(J(E,Ex*xp,0,Exx*xp*xp+Ex*xpp,0,0))
    return eta,xs[0],xs[1],es[0],es[1]
def evaluate(el,eh,sl,sh):
    eta,x,y,Ex,Ey=data(el,eh);ss=J(I(sl,sh),0,1);gap=y-x
    ax=x*ss;ay=y*ss;aa=J(ALPHA)*ss
    X=[BN(j,ax) for j in range(3)];Y=[BN(j,ay) for j in range(3)]
    V=(X[0]/Ex-Y[0]/Ey)/gap
    V1=(-x*X[1]/Ex+y*Y[1]/Ey)/(4*gap)
    V2=(x*x*X[2]/Ex-y*y*Y[2]/Ey)/(16*gap)
    H=(X[1]/Ex-Y[1]/Ey)/gap
    ha=hj1(ALPHA);W=BN(0,aa)/ha;W1=-J(ALPHA)*BN(1,aa)/(4*ha)
    t0=-V/(2*V1);a0=-1-2*ss*V2/V1;T0=t0*a0-ss
    f0=-eta*W/(2*V1);w=-eta*W1/V1
    c=x+y;d=x*y;R=d*H/(4*V1)
    G=1+a0-a0*a0-c*ss-(2*ss-t0)*R+d*t0*t0+ss*w*w+f0*w+J(ALPHA)*f0*f0
    D=(-G*(2*ss+t0)+t0*(t0*R+a0+f0*w+2))/T0
    h=eta*(V*W1-W*V1)/(V*V1+2*ss*V*V2-2*ss*V1*V1)
    return [h,D,-T0,V,-V1]
def cell(el,eh,sl,sh):
    full=evaluate(el,eh,sl,sh);em=(el+eh)/2;sm=(sl+sh)/2;mid=evaluate(em,em,sm,sm)
    he=(eh-el)/2;hs=(sh-sl)/2;de=I(-he,he);ds=I(-hs,hs);de2=I(0,he*he);ds2=I(0,hs*hs)
    raw=[c.v+c.e*de+c.s*ds+r.ee*de2/2+r.es*de*ds+r.ss*ds2/2 for c,r in zip(mid,full)]
    aa=1-I(sl,sh)*raw[0]*raw[0]
    out=[aa,4*aa-raw[1]*raw[1],*raw[2:]]
    bad=next((i for i,v in enumerate(out) if v.lo<=0),0)
    indices=[0,1] if bad<2 else [bad]
    ee=sum((mid[i].e*de+full[i].ee*de2/2).absupper() for i in indices)
    es=sum((mid[i].s*ds+full[i].ss*ds2/2).absupper() for i in indices)
    return out,ee>=es

def main():
    start=time.time();audit_endpoints();denominator_certificate();assert ZETA.hi<36*Q and (BETA-ALPHA).lo>8*Q
    todo=[(F(0),F(1),F(0),SMAX,0)];good=[];nodes=0;lower=[None]*5;upper=[None]*5;maxdepth=0
    while todo:
        el,eh,sl,sh,dep=todo.pop();nodes+=1
        try:res,split=cell(el,eh,sl,sh)
        except ZeroDivisionError:res=None;split=(eh-el)*SMAX>sh-sl
        if res is not None and all(v.lo>0 for v in res):
            good.append([str(el),str(eh),str(sl),str(sh)])
            for i,v in enumerate(res):
                lower[i]=v.lo if lower[i] is None else min(lower[i],v.lo);upper[i]=v.hi if upper[i] is None else max(upper[i],v.hi)
            maxdepth=max(maxdepth,dep)
        else:
            assert dep<40,(el,eh,sl,sh,None if res is None else [v.approx() for v in res])
            if split:
                mid=(el+eh)/2;todo.extend([(el,mid,sl,sh,dep+1),(mid,eh,sl,sh,dep+1)])
            else:
                mid=(sl+sh)/2;todo.extend([(el,eh,sl,mid,dep+1),(el,eh,mid,sh,dep+1)])
        if nodes%1000==0:print('nodes',nodes,'accepted',len(good),'pending',len(todo),'eta',float(el),float(eh),'s',float(sl),float(sh),'seconds',time.time()-start,flush=True)
    out={'status':'exact_center_PSD_assertions_passed','eta_interval':['0','1'],'s_interval':['0',str(SMAX)],'method':'exact second-order Taylor model of h and D, followed by PSD inequalities','bessel_degree':K,'tail':str(TAIL),'cells':len(good),'nodes':nodes,'max_depth':maxdepth,'certified_bounds':[[str(F(a,Q)),str(F(b,Q))] for a,b in zip(lower,upper)],'bounds_approximate':[[a/Q,b/Q] for a,b in zip(lower,upper)],'quantities':['a','4a-D^2','-T0','V','-Vprime'],'partition':good,'seconds':time.time()-start}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='partition'},indent=2),flush=True)
if __name__=='__main__':main()
