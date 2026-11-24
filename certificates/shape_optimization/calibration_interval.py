"""Exact fixed-point interval certificate for the nonnegative calibration profile.

All numerical enclosures use directed integer arithmetic at a fixed common
scale. This file uses the Python standard library only. It certifies scalar
inequalities; the differential calibration identity and integration arguments
are separate mathematical steps.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json,time

DIGITS=100
Q=10**DIGITS

def ceildiv(a,b):return -((-a)//b)

class I:
 __slots__=('lo','hi')
 def __init__(self,x=0,y=None,raw=False):
  if raw:self.lo,self.hi=x,y;return
  if isinstance(x,I):self.lo,self.hi=x.lo,x.hi;return
  x=F(x);y=x if y is None else F(y)
  self.lo=x.numerator*Q//x.denominator;self.hi=ceildiv(y.numerator*Q,y.denominator)
 def __add__(self,b):
  b=I(b);return I(self.lo+b.lo,self.hi+b.hi,True)
 __radd__=__add__
 def __neg__(self):return I(-self.hi,-self.lo,True)
 def __sub__(self,b):return self+-I(b)
 def __rsub__(self,b):return I(b)+-self
 def __mul__(self,b):
  b=I(b);p=[self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi]
  return I(min(p)//Q,ceildiv(max(p),Q),True)
 __rmul__=__mul__
 def __truediv__(self,b):
  b=I(b)
  if b.lo<=0<=b.hi:raise ZeroDivisionError((b.lo,b.hi))
  lower=[];upper=[]
  for aa in (self.lo,self.hi):
   for bb in (b.lo,b.hi):
    numerator,denominator=(aa,bb) if bb>0 else (-aa,-bb)
    lower.append(numerator*Q//denominator);upper.append(ceildiv(numerator*Q,denominator))
  return I(min(lower),max(upper),True)
 def __rtruediv__(self,b):return I(b)/self
 def __pow__(self,k):
  assert isinstance(k,int) and k>=0
  if k==0:return I(1)
  if k==1:return self
  b=self**(k//2);out=b*b
  return out*self if k%2 else out
 def absupper(self):return max(abs(self.lo),abs(self.hi))
 def expand(self,e):
  if isinstance(e,I):e=e.absupper()
  elif isinstance(e,F):e=ceildiv(e.numerator*Q,e.denominator)
  else:e=I(e).absupper()
  return I(self.lo-e,self.hi+e,True)
 def approx(self):return [self.lo/Q,self.hi/Q]
 def __repr__(self):return str(self.approx())

ZERO=I(0)
def padd(*polys):
 out=[I(0)]*max(map(len,polys))
 for p in polys:
  for j,c in enumerate(p):out[j]=out[j]+c
 return out

def pscale(p,c):return [a*c for a in p]
def pshift(p,k):return [I(0)]*k+p

def pmul(a,b):
 out=[I(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  if x.lo==x.hi==0:continue
  for j,y in enumerate(b):
   if y.lo!=0 or y.hi!=0:out[i+j]=out[i+j]+x*y
 return out

def pder(p):return [j*p[j] for j in range(1,len(p))]
def peval(p,s):
 out=I(0)
 for c in reversed(p):out=out*s+c
 return out

def j0(z,K=80):
 term=I(1);total=term
 for k in range(1,K+1):term=term*(-z)/(4*k*k);total=total+term
 nextterm=term*(-z)/(4*(K+1)**2)
 err=I(F(nextterm.absupper(),Q))/(1-z/(4*(K+2)**2))
 return total.expand(err)

def hj1(z,K=80):
 term=I(1);total=term
 for k in range(1,K+1):term=term*(-z)/(4*k*(k+1));total=total+term
 nextterm=term*(-z)/(4*(K+1)*(K+2))
 err=I(F(nextterm.absupper(),Q))/(1-z/(4*(K+2)*(K+3)))
 return z/2*total.expand(err)

ROOT_DEN=10**80
XL=578318596294678452117599575845580703507144180642368558708712371445606430488554437
YL=3047126234366208639907816317502275584842083370855398235487607831463637089182269021


def setup(K=60):
 xl,xh=F(XL,ROOT_DEN),F(XL+1,ROOT_DEN)
 yl,yh=F(YL,ROOT_DEN),F(YL+1,ROOT_DEN)
 assert j0(I(xl)).lo>0 and j0(I(xh)).hi<0
 assert j0(I(yl)).hi<0 and j0(I(yh)).lo>0
 # Root indexing is proved in the companion note: J0>0 on[0,2], and
 # consecutive zeros above2 have spacing>=4*pi/sqrt17>2, so j03>6.
 assert F(5)<xl<xh<F(6) and F(30)<yl<yh<F(31)
 x=I(xl,xh);y=I(yl,yh);c=x+y;d=x*y
 hx=hj1(x);hy=hj1(y)
 assert hx.lo>Q and hy.hi< -Q # |sqrt(z)J1(sqrt(z))|>1.
 assert (y-x).lo>Q and y.hi<36*Q
 # Center: u(r)=V(r^2), V=(J0(sqrt(x)sqrt(s))/hx - .../hy)/(y-x).
 vx=[I(1)];vy=[I(1)]
 for k in range(1,K+1):
  vx.append(vx[-1]*(-x)/(4*k*k));vy.append(vy[-1]*(-y)/(4*k*k))
 V=[(a/hx-b/hy)/(y-x) for a,b in zip(vx,vy)]
 center=[V]
 for j in range(3):center.append(pder(center[-1]))
 sm=F(9,16)
 # Uniform tail for all derivatives of orders 0..3 on 0<=s<=9/16.
 ratio=9*sm*F((K+2)**3,(K+1)**3*(K+2)**2)
 tail_center=2*F(9**(K+1)*(K+1)**3,factorial(K+1)**2)*sm**(K-2)/(1-ratio)
 # Endpoint: s=1-r, u(s)=sum aj*s^j, normalized u'''(1)=-1.
 # [(1-s)^3 D^4 -2(1-s)^2 D^3 +(c(1-s)^3-(1-s))D^2
 #  -(c(1-s)^2+1)D+d(1-s)^3]u=0.
 rpoly=[I(1),I(-1)];r2=pmul(rpoly,rpoly);r3=pmul(r2,rpoly)
 ops={4:r3,3:pscale(r2,-2),2:padd(pscale(r3,c),pscale(rpoly,-1)),
      1:padd(pscale(r2,-c),[I(-1)]),0:pscale(r3,d)}
 coeff=[I(0),I(0),I(0),I(F(1,6))]
 for k in range(K-3):
  val=I(0)
  for der,poly in ops.items():
   for j,fac in enumerate(poly):
    index=k-j+der
    if index<der or index>=len(coeff):continue
    mult=1
    for h in range(der):mult*=index-h
    val=val+fac*mult*coeff[index]
  coeff.append(-val/((k+1)*(k+2)*(k+3)*(k+4)))
 coeff[4]=I(F(1,12)) # exact ODE consequence, also used in zero-factor audit.
 U=coeff[3:];P=[(k+3)*v for k,v in enumerate(U)];QQ=[(k+2)*v for k,v in enumerate(P)]
 LL=[(k+1)*v for k,v in enumerate(QQ)]
 KK=[]
 for k in range(len(U)+1):KK.append(((U[k] if k<len(U) else I(0))-(U[k-1] if k else I(0)))/(k+4))
 A=padd(pmul(rpoly,P),pscale(pshift(U,1),-1))
 Num=padd(pmul(r2,padd(pmul(P,LL),pscale(pmul(QQ,QQ),-1))),A,
          pscale(pshift(pmul(KK,A),4),-d),pscale(pshift(pmul(r2,pmul(U,U)),4),d))
 # These two coefficients vanish EXACTLY, not merely by interval inclusion:
 # U0=1/6,U1=1/12; P0=1/2,P1=1/3; Q0=L0=1,Q1=1,L1=2.
 # constant=(1/2-1)+1/2=0;
 # linear=(4/3-2)-2*(1/2-1)+(1/3-1/2-1/6)=0.
 # The following assertions only audit consistency of their enclosures.
 assert Num[0].lo<=0<=Num[0].hi and Num[1].lo<=0<=Num[1].hi
 W4=Num[2:]
 endpoint={'U':U,'P':P,'Q':QQ,'L':LL,'K':KK,'W4':W4}
 smax=F(1,4);M=K-3
 # All coefficient tails U,P,Q,L,K are bounded by E(s)=432(6s)^(M+1)
 # /((M+1)! (1-6smax/(M+2))). This follows from |u^(j)|<=2*6^j.
 Econst=F(432*6**(M+1),factorial(M+1))/(1-6*smax/F(M+2))
 Emax=Econst*smax**(M+1)
 assert Emax<F(1,10**20)
 # Exact coarse field bounds used in the W4 product remainder estimate.
 for name in ('U','P','Q','L','K'):
  bound=sum(F(q.absupper(),Q)*smax**j for j,q in enumerate(endpoint[name]))+Emax
  assert bound<10000
 assert d.hi<200*Q
 # Difference of W4 numerators <=1e8*E(s); after division by s^2:
 Werror=10**8*Econst*smax**(M-1)
 return c,d,center,tail_center,endpoint,Emax,Werror


def center_cell(interval,data):
 c,d,polys,tail,_,_,_=data;s=I(*interval)
 V,V1,V2,V3=[peval(p,s).expand(tail) for p in polys]
 if V.lo<=0 or V1.hi>=0:return None
 t0=-V/(2*V1);a=-1-2*s*V2/V1;T0=t0*a-s
 if T0.hi>=0:return None
 R=-c-4*(2*V2+s*V3)/V1
 G=1+a-a*a-c*s-(2*s-t0)*R+d*t0*t0
 D=(-G*(2*s+t0)+t0*(a+t0*R+2))/T0
 return D,T0,V,-V1


def endpoint_cell(interval,data):
 c,d,_,_,polys,tail,Werr=data;s=I(*interval);r=1-s
 U,P,QQ,K=[peval(polys[name],s).expand(tail) for name in ('U','P','Q','K')]
 W4=peval(polys['W4'],s).expand(Werr)
 Z=P*P-U*QQ
 if U.lo<=0 or P.lo<=0 or Z.lo<=0:return None
 D=((2*r*P+s*U)*W4+U*U-d*s**4*U*U*K-r*U*P*QQ-2*s*U*P*P)/(r*P*Z)
 return D,Z,U,P


def partition(kind,data,bound=F(3,2),min_bound=F(0),max_depth=22):
 end=F(9,16) if kind=='center' else F(1,4);fun=center_cell if kind=='center' else endpoint_cell
 pending=[(F(0),end,0)];accepted=[];failed=0
 while pending:
  a,b,depth=pending.pop()
  try:res=fun((a,b),data)
  except ZeroDivisionError:res=None
  if res is not None and res[0].lo>min_bound*Q and res[0].hi<bound*Q:
   accepted.append((a,b,res[0].lo,res[0].hi));continue
  if depth>=max_depth:
   print('FAILED',kind,float(a),float(b),None if res is None else res[0].approx(),flush=True)
   raise AssertionError('interval subdivision did not certify the cell')
  mid=(a+b)/2;pending.extend([(a,mid,depth+1),(mid,b,depth+1)])
 return {'kind':kind,'cells':len(accepted),'D_bounds_approximate':[min(v[2] for v in accepted)/Q,max(v[3] for v in accepted)/Q],
         'smallest_cell_approximate':float(min(b-a for a,b,_,_ in accepted)),
         'partition':[[str(a),str(b)] for a,b,_,_ in accepted]}

if __name__=='__main__':
 t=time.time();data=setup();print('setup',time.time()-t,'tails',float(data[3]),float(data[5]),float(data[6]),flush=True)
 rows=[]
 for kind in ['center','endpoint']:
  row=partition(kind,data);rows.append(row);print({k:v for k,v in row.items() if k!='partition'},'seconds',time.time()-t,flush=True)
 result={'status':'all_exact_interval_assertions_passed','digits':DIGITS,'polynomial_degree':60,
         'D_strict_bounds':['0','3/2'],'center_tail':str(data[3]),'endpoint_field_tail':str(data[5]),'endpoint_W4_tail':str(data[6]),
         'partitions':rows,'seconds':time.time()-t}
 Path(__file__).with_name('calibration_interval_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
 print('DONE',time.time()-t,flush=True)
