"""Resolved corner B after the regular eta>=1/16 chart is certified.

Only the subset eta=s*delta<=1/16 needs a fresh sign check. Boxes wholly
above that threshold are inherited from E/Esmall; in crossing boxes the
Robin root intervals are intersected with eta<=1/16. Every exact root
parameter used by the subset remains enclosed.
"""
from pathlib import Path
from fractions import Fraction as F
import time,json
import all_contrast_endpoint_fast as fast
ep=fast.ep
CUT=F(1,16)
original_data=fast.data
def capped_data(lo,hi):
    assert lo<CUT
    return original_data(lo,min(hi,CUT))
ep.data=capped_data

if __name__=='__main__':
    folder=Path(__file__).resolve().parent
    inherited_files=['all_contrast_endpoint_fast_E.json','all_contrast_endpoint_fast_Esmall.json']
    for filename in inherited_files:
        row=json.loads((folder/filename).read_text())
        assert row['covered_fraction']=='1'
    ep.audit_endpoints();fast.audit_W_tail()
    todo=[(F(0),F(1,4),F(0),F(1),0,0)]
    good=[];inherited=[];count=0;covered=F(0);total=F(1,4);start=time.time()
    while todo:
        a,b,c,d,da,db=todo.pop();count+=1
        if a*c>=CUT:
            inherited.append([str(a),str(b),str(c),str(d)]);covered+=(b-a)*(d-c)
            continue
        try:res=ep.chart_B((a,b),(c,d))
        except (ZeroDivisionError,ep.PolynomialWidthError):res=None
        if res is not None and all(v.lo>0 for v in res):
            good.append([str(a),str(b),str(c),str(d)]);covered+=(b-a)*(d-c)
        else:
            assert da+db<30,(a,b,c,d,None if res is None else [v.approx() for v in res])
            if da<=db:
                m=(a+b)/2;todo.extend([(a,m,c,d,da+1,db),(m,b,c,d,da+1,db)])
            else:
                m=(c+d)/2;todo.extend([(a,b,c,m,da,db+1),(a,b,m,d,da,db+1)])
        if count%200==0:
            print('PROGRESS_B_COVERED',count,len(good),len(inherited),len(todo),'covered',float(covered/total),'seconds',time.time()-start,flush=True)
    assert covered==total
    row=dict(status='all_exact_resolved_B_assertions_passed',
             chart='B_with_regular_overlap',eta_cutoff=str(CUT),
             s_interval=['0','1/4'],delta_interval=['0','1'],covered_fraction='1',
             inherited_certificates=inherited_files,
             cells=len(good),inherited_cells=len(inherited),tested=count,seconds=time.time()-start,
             polynomial_degree=ep.M,common_error_bound=str(ep.ERR),W_tail_bound=str(fast.audit_W_tail()),
             partition=good,inherited_partition=inherited)
    Path(__file__).with_suffix('.json').write_text(json.dumps(row,indent=2)+'\n')
    print('PASSED_B_COVERED',{k:v for k,v in row.items() if k not in('partition','inherited_partition')},flush=True)
