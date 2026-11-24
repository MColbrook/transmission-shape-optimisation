"""Exact rational certificate rho(Omega,N) < 5.1963, |Omega|=pi.

Python standard library only. Map and trial coefficients are printed integers;
no numerical data files, numerical eigensolvers, or numerical quadrature are used.
The inverse derivative is approximated by a fixed-denominator polynomial H and
its error is certified by the exact polynomial residual 1-gH.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json,time

# f(z)=S^(-1/2)(z+sum a_j z^(2j+3)), S=1+sum(2j+3)a_j^2.
MAP_DEN=10**8
MAP_NUM=[23661008,944390,841714,-267345,182938,-115928,80287,-57758,42974,-32778,25473,-19906]
# U=(1-r^2)^2 sum_m sum_k b_mk r^m cos(m theta) P_k^(4,m)(2r^2-1).
TRIAL_DEN=10**10
TRIAL_NUM=[10000000000,-714354563,-56400458,-2125776,1583967,
 3570773036,-529889309,23843272,-597204,292049,
 785588028,-133167479,9457497,-308100,59819,
 134049214,-29915606,2893523,-135368,8429,
 29725599,-4266203,416310,-17814,2314,
 -1030520,-1130923,153893,-11423,128,
 6047546,300537,6035,2833,353]
NR=5;NA=7


def mul_analytic(p,h):
 out={}
 for (a,b),v in p.items():
  for k,w in h.items():out[a+k,b]=out.get((a+k,b),F(0))+v*w
 return {key:v for key,v in out.items() if v}


def norm(p):
 """Exact integral of |p(z,zbar)|^2 over D, divided by pi."""
 groups={}
 for (a,b),v in p.items():groups.setdefault(a-b,[]).append((a,b,v))
 ans=F(0)
 for terms in groups.values():
  for i,(a,b,v) in enumerate(terms):
   ans+=v*v/F(a+b+1)
   for c,d,w in terms[i+1:]:ans+=2*v*w/F(a+d+1)
 return ans


def inner(p,q):
 ans=F(0)
 for (a,b),v in p.items():
  for (c,d),w in q.items():
   if a+c==b+d:ans+=v*w/F(a+c+1)
 return ans


def jacobi_shifted(k,m):
 """Integer coefficients of P_k^(4,m)(2s-1)."""
 out=[0]*(k+1)
 for j in range(k+1):
  v=comb(k+4,j)*comb(k+m,k-j)
  for ell in range(k-j+1):out[j+ell]+=v*comb(k-j,ell)*(-1)**(k-j-ell)
 return out


def trial():
 U={}
 for j in range(NA):
  m=2*j+1;p=[F(0)]*NR
  for k in range(NR):
   for ell,c in enumerate(jacobi_shifted(k,m)):p[ell]+=F(TRIAL_NUM[j*NR+k],TRIAL_DEN)*c
  q=[F(0)]*(NR+2)
  for ell,c in enumerate(p):
   q[ell]+=c;q[ell+1]-=2*c;q[ell+2]+=c
  for ell,c in enumerate(q):
   U[ell+m,ell]=c/2;U[ell,ell+m]=c/2
 lap={(a-1,b-1):4*a*b*v for (a,b),v in U.items() if a and b and v}
 return U,lap


def residual_inverse(ci,K=300,T=10**24):
 """Construct H(w), then bound its error by exact residual R=1-gH."""
 d=[T]
 for k in range(1,K+1):
  numerator=-sum(ci[j-1]*d[k-j] for j in range(1,min(k,len(ci))+1))
  q,r=divmod(numerator,MAP_DEN)
  d.append(q+(2*r>=MAP_DEN))
 residual=[]
 for k in range(K+len(ci)+1):
  v=MAP_DEN*d[k] if k<=K else 0
  v+=sum(ci[j-1]*d[k-j] for j in range(1,len(ci)+1) if 0<=k-j<=K)
  if k==0:v-=MAP_DEN*T
  residual.append(v)
 R=F(sum(abs(x) for x in residual),MAP_DEN*T)
 H={2*k:F(v,T) for k,v in enumerate(d) if v}
 return H,R


def j0partial(x,K):
 term=F(1);ans=term
 for k in range(1,K+1):term*=-x/F(4*k*k);ans+=term
 return ans


def main():
 start=time.time();a=[F(x,MAP_DEN) for x in MAP_NUM];orders=[2*j+3 for j in range(len(a))]
 S=1+sum(o*x*x for o,x in zip(orders,a));ci=[o*x for o,x in zip(orders,MAP_NUM)]
 eta=F(sum(abs(x) for x in ci),MAP_DEN)
 assert eta<1 # Re f'>0 proves injectivity on the closed disk and a smooth boundary.
 g={0:F(1)};g.update({2*(j+1):F(c,MAP_DEN) for j,c in enumerate(ci)})
 U,lap=trial()
 B=-inner(U,lap);C=norm(mul_analytic(U,g))/S;Adisk=norm(lap)
 H,R=residual_inverse(ci)
 err=R/(1-eta)
 Apoly=S*norm(mul_analytic(lap,H))
 Aerror=S*Adisk*(2*err/(1-eta)+err*err)
 Aupper=Apoly+Aerror
 tau=F(51963,10000);n0=F(263447,50000)
 slope=tau*tau*C-tau*B
 qupper=Aupper-(n0+1)*tau*B+n0*tau*tau*C
 assert slope<0
 assert qupper < -F(1,50000)
 # Compact rational bounds suitable for printing in the paper.
 Blo=F(2941129141,10**10);Chi=F(7471409,400000000);Ahi=F(6923409,10**6)
 assert B>Blo and C<Chi and Aupper<Ahi
 compact_q=Ahi-(n0+1)*tau*Blo+n0*tau*tau*Chi
 assert compact_q== -F(44004164878243713,2*10**21)
 assert compact_q < -F(1,50000)
 assert err<F(1,10**10)
 # To isolate the first two zeros, first bound their spacing.
 # J0(x)>0 on [0,2] by pairing its alternating series. Between zeros x>=2,
 # y=sqrt(x)J0(x) satisfies y''+(1+1/(4x^2))y=0. Wirtinger gives each
 # zero gap >=4*pi/sqrt(17)>2. Hence j03>6. A negative J0 at sqrt(bl)<6,
 # after the already located first zero, must therefore precede j02.
 # Odd/even partial sums enclose the tails at the following arguments.
 ah=F(5783185963,10**9);bl=F(30471262343,10**9)
 assert F(4)<ah<bl<F(36)
 assert j0partial(ah,20)<0
 assert j0partial(bl,20)<0
 assert bl/ah>n0
 result={
  'status':'all_exact_rational_assertions_passed','tau':str(tau),'n0':str(n0),
  'map_denominator':MAP_DEN,'map_numerators':MAP_NUM,
  'trial_denominator':TRIAL_DEN,'trial_numerators':TRIAL_NUM,
  'S':str(S),'eta':str(eta),'B_over_pi':str(B),'C_over_pi':str(C),
  'compact_B_lower_bound':str(Blo),'compact_C_upper_bound':str(Chi),'compact_I2_upper_bound':str(Ahi),
  'compact_q_upper_bound':str(compact_q),
  'I2_upper_over_pi_approximate':float(Aupper),'I2_error_bound_over_pi_approximate':float(Aerror),
  'q_upper_over_pi_approximate':float(qupper),'q_upper_over_pi_strict_rational_bound':'-1/50000',
  'reciprocal_residual_norm_approximate':float(R),'reciprocal_uniform_error_approximate':float(err),
  'seconds':time.time()-start
 }
 output=Path(__file__).with_name('near_optimal_domain_certificate.json');output.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))

if __name__=='__main__':main()
