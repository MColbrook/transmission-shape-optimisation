"""Efficient exact endpoint partition with a shorter, audited Taylor order.

Only the W recurrence is used, with its omitted coefficients bounded by
an explicit positive exponential majorant. The double factor in a is
removed from Z-H before multiplying by Z+H.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
from functools import lru_cache
import json,time,sys
import all_contrast_endpoint as ep

ep.M=38;ep.SMAX=F(1,4)
ep.EC=F(432*6**(ep.M+1),factorial(ep.M+1))/(1-6*ep.SMAX/F(ep.M+2))
ep.ERR=10**22*ep.EC*ep.SMAX**(ep.M-6)
assert ep.ERR<F(1,10**9)

I=ep.I;Q=ep.Q;padd=ep.padd;pscale=ep.pscale;pshift=ep.pshift;pder=ep.pder

def mul(a,b):
    # Only coefficients needed by the W recurrence are formed.
    out=[I(0)]*min(len(a)+len(b)-1,ep.M+1)
    for j,x in enumerate(a):
        if j>=len(out):break
        if not x.absupper():continue
        for k,y in enumerate(b[:len(out)-j]):
            if y.absupper():out[j+k]=out[j+k]+x*y
    return out

@lru_cache(maxsize=8192)
def data(lo,hi):
    x,y=ep.root_boxes(lo,hi);c=x+y;d=x*y;m=ep.M
    ua=[I(0),I(0),I(F(1,2)),I(0)];ub=[I(0),I(0),I(0),I(F(1,6))]
    r=[I(1),I(-1)];r2=mul(r,r);r3=mul(r2,r)
    ops={4:r3,3:pscale(r2,-2),2:padd(pscale(r3,c),pscale(r,-1)),
         1:padd(pscale(r2,-c),[I(-1)]),0:pscale(r3,d)}
    for k in range(m+1):
        for coeff in (ua,ub):
            v=I(0)
            for der,poly in ops.items():
                for j,fac in enumerate(poly):
                    index=k-j+der
                    if index<der or index>=len(coeff):continue
                    mult=1
                    for h in range(der):mult*=index-h
                    v=v+fac*mult*coeff[index]
            coeff.append(-v/((k+1)*(k+2)*(k+3)*(k+4)))
    phi=[I(0),I(1)]
    for k in range(m+2):
        phi.append(((k+1)**2*phi[k+1]-ep.ALPHA*phi[k]+(ep.ALPHA*phi[k-1] if k else 0))/((k+2)*(k+1)))
    def fields(coeff):
        U=coeff[2:m+3]
        return dict(U=U,
          P=[(k+2)*coeff[k+2] for k in range(m+1)],
          Q=[(k+1)*(k+2)*coeff[k+2] for k in range(m+1)],
          L=[(k+1)*(k+2)*(k+3)*coeff[k+3] for k in range(m+1)],
          K=[(U[k]-(U[k-1] if k else 0))/(k+3) for k in range(m+1)])
    fa,fb=fields(ua),fields(ub);FF=phi[1:m+2];J=[(k+1)*phi[k+1] for k in range(m+1)]
    Ua,Pa,Qa,Ka=[fa[k] for k in ('U','P','Q','K')]
    Ub,Pb,Qb,Kb=[fb[k] for k in ('U','P','Q','K')]
    M0=padd([I(1)],pscale(pshift(Kb,3),-d));M1=padd([I(-1)],pscale(pshift(Ka,3),-d))
    h2=padd(mul(Qa,Pa),pscale(mul(FF,J),-1))[1:]
    h1=padd(mul(Qa,Pb),mul(Qb,Pa))[1:];h0=mul(Qb,Pb)[1:]
    Rs=[padd(mul(r,h0),pscale(mul(Pb,Pb),2),pscale(mul(M0,Ub),-1)),
        padd(mul(r,h1),pscale(mul(Pa,Pb),4),pscale(mul(M0,Ua),-1),pscale(mul(M1,Ub),-1)),
        padd(mul(r,h2),pscale(mul(Pa,Pa),2),pscale(mul(M1,Ua),-1))]
    Ws=[]
    for j in range(3):
        w=[I(0) if j==0 else I(F(1,3))]
        for k in range(1,m-2):w.append((Rs[j][k]+(k+1)*w[k-1])/(k+3))
        Ws.append(w)
    return c,d,fa,fb,FF,J,Ws
ep.data=data

def positive_poly_add(*polys):
    out=[F(0)]*max(map(len,polys))
    for p in polys:
        for j,v in enumerate(p):out[j]+=v
    return out
def shift(p,k):return [F(0)]*k+p
def scale(p,c):return [c*v for v in p]

def audit_W_tail():
    # N=s^3 W is majorized by P6(s)e^(6s)+P12(s)e^(12s).
    r=[F(1),F(1)];r2=[F(1),F(2),F(1)]
    P6=positive_poly_add(scale(r2,9),scale(shift(r,1),9),scale(shift(r,1),36),
                        [F(0),F(0),F(18)],scale(shift(r2,2),54))
    P12=positive_poly_add(scale(r2,1296),scale(shift(r2,1),7776),
                         scale(shift(r2,4),64800),scale(shift(r,4),57600),
                         [F(0)]*5+[F(28800)])
    a=ep.SMAX;K=ep.M+1;nv=F(0);nd=F(0)
    def tail(z,k):
        assert k>z
        return z**k/F(factorial(k))/(1-z/F(k+1))
    for z,p in ((F(6),P6),(F(12),P12)):
        for j,v in enumerate(p):
            if not v:continue
            k=K-j;tv=tail(z*a,k)
            nv+=v*a**j*tv
            nd+=v*a**j*(z*a*tail(z*a,k-1)+j*tv)
    # W is retained through degree M-3. The value and derivative tails
    # are bounded by nv/a^3 and nd/a^4. Chart A/B divisions cost at most4.
    extra=4*nv/a**3+nd/a**4
    assert extra<ep.ERR/F(100)
    return extra

def W_value_tail():
    # The derivative-inclusive budget is also a valid (loose) bound for
    # W/s on chart B, by the preceding calculation.
    return audit_W_tail()

def chart_B_poly(sab,dab):
    ss=I(*sab);delta=I(*dab);eta=ss*delta
    elo=sab[0]*dab[0];ehi=sab[1]*dab[1]
    c,d,fa,fb,FF,J,Ws=data(elo,ehi);p={}
    for key in ('U','P','Q','K'):
        assert fb[key][0].lo==fb[key][0].hi==0
        p[key]=padd(pscale(fa[key],delta),fb[key][1:])
    p.update(F=FF,J=J)
    W=padd(Ws[0][1:],pscale(Ws[1],delta),pscale(pshift(Ws[2],1),delta**2))
    W[0]=(1+4*delta)/12;p['W']=W
    # Pointwise analytic errors before any rational operation.
    eraw=ep.EC*ep.SMAX**ep.M;ew=W_value_tail()
    BB=F(1)+max(sum(F(v.absupper(),Q)*ep.SMAX**j for j,v in enumerate(poly)) for poly in p.values())
    if BB>10**6:return None
    U,P,QQ,KK,Fp,JJ,WW=[p[k] for k in ('U','P','Q','K','F','J','W')]
    pm=ep.pmul;r=[I(1),I(-1)]
    Zpol=padd(pm(P,P),pscale(pm(U,QQ),-1))
    Hpol=pscale(padd(pm(Fp,P),pscale(pm(U,JJ),-1)),delta)
    Mpol=padd([I(1),-delta],pscale(pshift(KK,4),-d))
    Npol=padd(pscale(pm(r,pm(P,WW)),2),pm(Mpol,pm(U,U)),
              pm(r,pm(U,padd(pscale(pm(Fp,JJ),delta**2),pscale(pm(QQ,P),-1)))),
              pshift(pm(U,WW),1),pscale(pshift(pm(U,pm(P,P)),1),-2))
    ez=8*BB*eraw;eh=ez
    # N has coefficient-weight at most 208 and degree at most three;
    # telescoping products on fields bounded by 2*BB gives <2500 BB²
    # times the sum of the component errors. The factor below is larger.
    en=10000*BB**2*(eraw+ew)
    Uv=ep.peval(U,ss).expand(eraw);Pv=ep.peval(P,ss).expand(eraw)
    Fv=ep.peval(Fp,ss).expand(eraw);Zv=ep.peval(Zpol,ss).expand(ez)
    if min(Uv.lo,Pv.lo,Fv.lo,Zv.lo)<=0:return None
    diff=ep.peval(padd(Zpol,pscale(Hpol,-1)),ss).expand(ez+eh)
    plus=ep.peval(padd(Zpol,Hpol),ss).expand(ez+eh)
    aa=diff*plus/Zv**2
    DD=ep.peval(Npol,ss).expand(en)/((1-ss)*Pv*Zv)
    return aa,4*aa-DD**2,Zv

def chart_E(eab,sab):
    eta=I(*eab);ss=I(*sab);p=ep.original_polys(eta,data(*eab))
    ZZ=padd(ep.pmul(p['P'],p['P']),pscale(ep.pmul(p['U'],p['Q']),-1))
    HH=pscale(padd(ep.pmul(p['F'],p['P']),pscale(ep.pmul(p['U'],p['J']),-1)),eta)
    diff=padd(ZZ,pscale(HH,-1));plus=padd(ZZ,HH)
    for k in range(2):assert diff[k].lo<=0<=diff[k].hi
    ep.field_bound([p[k] for k in p]+[pder(p['P']),pder(p['W'])],ep.SMAX)
    U,P,QQ,FF,WW=[ep.evalp(p[k],ss) for k in ('U','P','Q','F','W')]
    PP=ep.evalp(pder(p['P']),ss);WP=ep.evalp(pder(p['W']),ss)
    Z=P**2-U*QQ
    if min(U.lo,P.lo,FF.lo,Z.lo)<=0:return None
    aa=ep.evalp(diff[2:],ss)*ep.evalp(plus,ss)/Z**2
    DD=2*WW/P**2-(U*P/Z)*(WW/P**2+ss*WP/P**2-2*ss*WW*PP/P**3)
    return aa,4*aa-DD**2,Z
ep.chart_E=chart_E

def certify(kind):
    if kind=='B':rect=(F(0),F(1,4),F(0),F(1));fun=chart_B_poly
    elif kind=='Esmall':rect=(F(1,16),F(1,4),F(0),F(1,4));fun=chart_E
    else:rect=(F(1,4),F(1),F(0),F(1,4));fun=chart_E
    total=(rect[1]-rect[0])*(rect[3]-rect[2]);covered=F(0)
    todo=[(*rect,0,0)];good=[];count=0;start=time.time()
    while todo:
        a,b,c,d,da,db=todo.pop();count+=1
        try:res=fun((a,b),(c,d))
        except (ZeroDivisionError,ep.PolynomialWidthError):res=None
        if res is not None and all(v.lo>0 for v in res):
            good.append([str(a),str(b),str(c),str(d)]);covered+=(b-a)*(d-c)
        else:
            assert da+db<30,(kind,a,b,c,d,None if res is None else [v.approx() for v in res])
            if da<=db:
                m=(a+b)/2;todo.extend([(a,m,c,d,da+1,db),(m,b,c,d,da+1,db)])
            else:
                m=(c+d)/2;todo.extend([(a,b,c,m,da,db+1),(a,b,m,d,da,db+1)])
        if count%200==0:
            print('PROGRESS',kind,count,len(good),len(todo),'covered',float(covered/total),'seconds',time.time()-start,flush=True)
    assert covered==total
    return dict(chart=kind,cells=len(good),tested=count,covered_fraction=str(covered/total),seconds=time.time()-start,
                partition=good,polynomial_degree=ep.M,common_error_bound=str(ep.ERR),W_tail_bound=str(audit_W_tail()))

if __name__=='__main__':
    ep.audit_endpoints();extra=audit_W_tail()
    print('ERRORS',float(ep.ERR),float(extra),flush=True)
    for kind in sys.argv[1:] or ['B','E']:
        row=certify(kind)
        Path(__file__).with_name('all_contrast_endpoint_fast_'+kind+'.json').write_text(json.dumps(row,indent=2)+'\n')
        print('PASSED',{k:v for k,v in row.items() if k!='partition'},flush=True)
