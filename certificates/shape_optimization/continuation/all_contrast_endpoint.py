"""Exact two-parameter endpoint certificates for n>N.

Three charts cover eta in [0,1], s=1-r in [0,1/4], with the singular
corner eta=s=0 resolved by s=eta*x and eta=s*delta.  Analytic identities
and tail majorants are documented in the companion proof note.
"""
from pathlib import Path
import sys,time,json
from fractions import Fraction as F
from math import factorial
from functools import lru_cache
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'minimizer_literature'/'continuation'))
from robin_roots import root_boxes,ALPHA,I,Q,j0,hj1,audit_endpoints
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from calibration_interval import padd,pmul,pscale,pshift,pder,peval

M=48
SMAX=F(1,4)
# The actual ubar and phi derivatives are <=2*6^j.  This field tail
# majorant covers U,P,Q,L,F,J,K.  The larger final margin also covers
# exact divisions and one derivative in the s=eta*x chart.
EC=F(432*6**(M+1),factorial(M+1))/(1-6*SMAX/F(M+2))
ERR=10**30*EC*SMAX**(M-6)
assert ERR<F(1,10**10)

def addmul(a,b,c):return padd(a,pscale(b,c))
def trunc(p,m=M+1):return p[:m]

@lru_cache(maxsize=8192)
def data(lo,hi):
    x,y=root_boxes(lo,hi);c=x+y;d=x*y
    # u coefficients are affine in eta: eta*ua+ub, at fixed c,d.
    ua=[I(0),I(0),I(F(1,2)),I(0)]
    ub=[I(0),I(0),I(0),I(F(1,6))]
    r=[I(1),I(-1)];r2=pmul(r,r);r3=pmul(r2,r)
    ops={4:r3,3:pscale(r2,-2),2:padd(pscale(r3,c),pscale(r,-1)),
         1:padd(pscale(r2,-c),[I(-1)]),0:pscale(r3,d)}
    for k in range(M+1):
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
    for k in range(M+2):
        phi.append(((k+1)**2*phi[k+1]-ALPHA*phi[k]+(ALPHA*phi[k-1] if k else 0))/((k+2)*(k+1)))
    def fields(coeff):
        U=coeff[2:M+3]
        P=[(k+2)*coeff[k+2] for k in range(M+1)]
        QQ=[(k+1)*(k+2)*coeff[k+2] for k in range(M+1)]
        L=[(k+1)*(k+2)*(k+3)*coeff[k+3] for k in range(M+1)]
        K=[(U[k]-(U[k-1] if k else 0))/(k+3) for k in range(M+1)]
        return dict(U=U,P=P,Q=QQ,L=L,K=K)
    fa,fb=fields(ua),fields(ub);FF=phi[1:M+2];J=[(k+1)*phi[k+1] for k in range(M+1)]
    Ua,Pa,Qa,La,Ka=[fa[k] for k in ('U','P','Q','L','K')]
    Ub,Pb,Qb,Lb,Kb=[fb[k] for k in ('U','P','Q','L','K')]
    Ra=padd(pmul(r,Pa),pscale(pshift(Ua,1),-1));Rb=padd(pmul(r,Pb),pscale(pshift(Ub,1),-1))
    N2=padd(pmul(r2,padd(pmul(J,J),pscale(pmul(Qa,Qa),-1))),
            pshift(padd(pmul(r2,pmul(Pa,La)),pscale(pmul(r,pmul(FF,J)),-1),pscale(Ra,-1)),1),
            pscale(pshift(pmul(r2,pmul(FF,FF)),2),ALPHA),
            pscale(pshift(padd(pmul(r2,pmul(Ua,Ua)),pscale(pmul(Ka,Ra),-1)),4),d))
    N1=padd(pscale(pmul(r2,pmul(Qa,Qb)),-2),
            pshift(padd(pmul(r2,padd(pmul(Pa,Lb),pmul(Pb,La))),Ra,pscale(Rb,-1)),1),
            pscale(pshift(padd(pscale(pmul(r2,pmul(Ua,Ub)),2),pscale(pmul(Ka,Rb),-1),pscale(pmul(Kb,Ra),-1)),4),d))
    N0=padd(pscale(pmul(r2,pmul(Qb,Qb)),-1),pshift(padd(pmul(r2,pmul(Pb,Lb)),Rb),1),
            pscale(pshift(padd(pmul(r2,pmul(Ub,Ub)),pscale(pmul(Kb,Rb),-1)),4),d))
    for nn in (N0,N1,N2):
        for k in range(3):assert nn[k].lo<=0<=nn[k].hi
    Ws=[N0[3:],N1[3:],N2[3:]]
    # A low-order recurrence reduces interval dependency dramatically.
    M0=padd([I(1)],pscale(pshift(Kb,3),-d))
    M1=padd([I(-1)],pscale(pshift(Ka,3),-d))
    h2=padd(pmul(Qa,Pa),pscale(pmul(FF,J),-1))[1:]
    h1=padd(pmul(Qa,Pb),pmul(Qb,Pa))[1:]
    h0=pmul(Qb,Pb)[1:]
    Rs=[padd(pmul(r,h0),pscale(pmul(Pb,Pb),2),pscale(pmul(M0,Ub),-1)),
        padd(pmul(r,h1),pscale(pmul(Pa,Pb),4),pscale(pmul(M0,Ua),-1),pscale(pmul(M1,Ub),-1)),
        padd(pmul(r,h2),pscale(pmul(Pa,Pa),2),pscale(pmul(M1,Ua),-1))]
    for j in range(3):
        Ws[j][0]=I(0) if j==0 else I(F(1,3))
        for k in range(1,M-2):Ws[j][k]=(Rs[j][k]+(k+1)*Ws[j][k-1])/(k+3)
    return c,d,fa,fb,FF,J,Ws

def original_polys(eta,dat):
    c,d,fa,fb,FF,J,Ws=dat
    fields={key:padd(pscale(fa[key],eta),fb[key]) for key in fa}
    fields.update(F=FF,J=J)
    W=padd(Ws[0],pscale(Ws[1],eta),pscale(Ws[2],eta**2))
    W[0]=eta*(eta+1)/3
    fields['W']=W
    return fields

class PolynomialWidthError(Exception):
    pass

def field_bound(polys,end):
    # Required hypotheses for the very conservative algebraic tail margin.
    for p in polys:
        if sum(F(v.absupper(),Q)*end**j for j,v in enumerate(p))>=10**6:
            raise PolynomialWidthError('coefficient enclosure requires subdivision')

def evalp(p,x):return peval(p,x).expand(ERR)

def chart_A(eab,xab):
    eta=I(*eab);xx=I(*xab);dat=data(*eab)
    orig=original_polys(eta,dat)
    p={}
    for key in ('U','P','Q','W'):
        first={'U':I(F(1,2)),'P':I(1),'Q':I(1),'W':(eta+1)/3}[key]
        p[key]=[first]+[v*eta**(j-1) for j,v in enumerate(orig[key]) if j]
    for key in ('F','J'):
        p[key]=[v*eta**j for j,v in enumerate(orig[key])]
    ZZ=padd(pmul(p['P'],p['P']),pscale(pmul(p['U'],p['Q']),-1))
    HH=padd(pmul(p['F'],p['P']),pscale(pmul(p['U'],p['J']),-1))
    anum=padd(pmul(ZZ,ZZ),pscale(pmul(HH,HH),-1))
    for k in range(2):assert anum[k].lo<=0<=anum[k].hi
    field_bound([p[k] for k in p]+[pder(p['P']),pder(p['W'])],F(1))
    U,P,QQ,FF,WW=[evalp(p[k],xx) for k in ('U','P','Q','F','W')]
    PP=evalp(pder(p['P']),xx);WP=evalp(pder(p['W']),xx)
    Z=P**2-U*QQ
    if min(U.lo,P.lo,FF.lo,Z.lo)<=0:return None
    a2=evalp(anum[2:],xx)/Z**2
    D1=2*WW/P**2-(U*P/Z)*(WW/P**2+xx*WP/P**2-2*xx*WW*PP/P**3)
    return a2,4*a2-D1**2,Z

def chart_E(eab,sab):
    eta=I(*eab);ss=I(*sab);dat=data(*eab)
    p=original_polys(eta,dat)
    ZZ=padd(pmul(p['P'],p['P']),pscale(pmul(p['U'],p['Q']),-1))
    HH=pscale(padd(pmul(p['F'],p['P']),pscale(pmul(p['U'],p['J']),-1)),eta)
    anum=padd(pmul(ZZ,ZZ),pscale(pmul(HH,HH),-1))
    for k in range(2):assert anum[k].lo<=0<=anum[k].hi
    field_bound([p[k] for k in p]+[pder(p['P']),pder(p['W'])],SMAX)
    U,P,QQ,FF,WW=[evalp(p[k],ss) for k in ('U','P','Q','F','W')]
    PP=evalp(pder(p['P']),ss);WP=evalp(pder(p['W']),ss)
    Z=P**2-U*QQ
    if min(U.lo,P.lo,FF.lo,Z.lo)<=0:return None
    a2=evalp(anum[2:],ss)/Z**2
    D1=2*WW/P**2-(U*P/Z)*(WW/P**2+ss*WP/P**2-2*ss*WW*PP/P**3)
    return a2,4*a2-D1**2,Z

def chart_B(sab,dab):
    ss=I(*sab);delta=I(*dab);eta=ss*delta
    elo=sab[0]*dab[0];ehi=sab[1]*dab[1]
    c,d,fa,fb,FF,J,Ws=data(elo,ehi)
    p={}
    for key in ('U','P','Q','K'):
        assert fb[key][0].lo==fb[key][0].hi==0
        p[key]=padd(pscale(fa[key],delta),fb[key][1:])
    p['L']=padd(fb['L'],pscale(pshift(fa['L'],1),delta))
    p.update(F=FF,J=J)
    assert Ws[0][0].lo==Ws[0][0].hi==0
    W=padd(Ws[0][1:],pscale(Ws[1],delta),pscale(pshift(Ws[2],1),delta**2))
    W[0]=(1+4*delta)/12;p['W']=W
    field_bound(list(p.values()),SMAX)
    U,P,QQ,FF,JJ,KK,WW=[evalp(p[k],ss) for k in ('U','P','Q','F','J','K','W')]
    Z=P**2-U*QQ;r=1-ss;HH=1-ss*delta-d*ss**4*KK
    if min(U.lo,P.lo,FF.lo,Z.lo)<=0:return None
    gg=delta*(FF*P-U*JJ)/Z;aa=1-gg**2
    DD=(2*WW+HH*U**2/(r*P)+U*(delta**2*FF*JJ-QQ*P)/P+ss*U*WW/(r*P)-2*ss*U*P/r)/Z
    return aa,4*aa-DD**2,Z

def certify(kind,maxdepth=24):
    if kind=='A':rect=(F(0),F(1,4),F(0),F(1));fun=chart_A
    elif kind=='B':rect=(F(0),F(1,4),F(0),F(1));fun=chart_B
    else:rect=(F(1,4),F(1),F(0),F(1,4));fun=chart_E
    todo=[(*rect,0,0)];good=[];count=0;start=time.time()
    while todo:
        a,b,c,d,da,db=todo.pop();count+=1
        try:res=fun((a,b),(c,d))
        except (ZeroDivisionError,PolynomialWidthError):res=None
        if res is not None and all(v.lo>0 for v in res):
            good.append([str(a),str(b),str(c),str(d)])
        else:
            if da+db>=maxdepth:
                print('FAILED',kind,a,b,c,d,None if res is None else [v.approx() for v in res],flush=True)
                raise AssertionError('cell failed')
            if da<=db:
                m=(a+b)/2;todo.extend([(a,m,c,d,da+1,db),(m,b,c,d,da+1,db)])
            else:
                m=(c+d)/2;todo.extend([(a,b,c,m,da,db+1),(a,b,m,d,da,db+1)])
        if count%100==0:print('PROGRESS',kind,count,len(good),len(todo),'seconds',time.time()-start,flush=True)
    return dict(chart=kind,cells=len(good),tested=count,seconds=time.time()-start,partition=good)

if __name__=='__main__':
    audit_endpoints();print('ERROR_BOUND',float(ERR),flush=True)
    kinds=sys.argv[1:] or ['A','B','E'];out=[]
    for kind in kinds:
        row=certify(kind);out.append(row);print('PASSED',{k:v for k,v in row.items() if k!='partition'},flush=True)
        row.update(polynomial_degree=M,common_error_bound=str(ERR))
        Path(__file__).with_name('all_contrast_endpoint_'+kind+'_M'+str(M)+'.json').write_text(json.dumps(row,indent=2)+'\n')
