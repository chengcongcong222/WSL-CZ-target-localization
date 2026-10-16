from pathlib import Path
p=Path('results/ENV_CAP_FULL_CAMPAIGN/scripts/e1_cold_audit.py');s=p.read_text(encoding='utf-8-sig')
old="  rays=[a for a in rayrows if a['case_id']==row['case_id']];raygood=len(rays)==2 and all(int(a['n_compared'])>=10 and float(a['max_return_range_difference_km'])<=.1 for a in rays)"
new="""  rays=[]
  for lev in ['R0','R1']:
   lines=(D/(lev+'.ray')).read_text().splitlines();ii=7;errors=[];aa=[float(v['angle_deg']) for v in ar];rs=[float(v['first_return_km']) for v in ar]
   while ii<len(lines):
    if not lines[ii].strip():ii+=1;continue
    angle=float(lines[ii]);n,ns,nb=map(int,lines[ii+1].split());ii+=2;pts=np.array([[float(v) for v in l.split()] for l in lines[ii:ii+n]]);ii+=n
    zz=pts[:,-1];hit=np.where((zz[:-1]>zs)&(zz[1:]<=zs))[0]
    if len(hit):
     jj=hit[0]
     if zz[:jj+2].min()>=1 and zz[:jj+2].max()<H-1:
      frac=(zs-zz[jj])/(zz[jj+1]-zz[jj]);r=(pts[jj,0]+frac*(pts[jj+1,0]-pts[jj,0]))/1000
      if aa[0]<=angle<=aa[-1] and 10<=r<=120:errors.append(abs(r-np.interp(angle,aa,rs)))
   saved=next(a for a in rayrows if a['case_id']==row['case_id'] and a['level']==lev)
   test(row['case_id']+'_'+lev+'_raw_ray_replay',len(errors)==int(saved['n_compared']) and abs(max(errors)-float(saved['max_return_range_difference_km']))<1e-8)
   rays.append({'n_compared':len(errors),'max_return_range_difference_km':max(errors) if errors else math.inf})
  raygood=all(a['n_compared']>=10 and a['max_return_range_difference_km']<=.1 for a in rays)"""
assert old in s;p.write_text(s.replace(old,new),encoding='utf-8',newline='\n')
