"""Saved finite-grid branch topology and error diagnostics, no search or observations."""
import json
import numpy as np,pandas as pd
import hla_h1_core as c
def components(ids):
    ids=np.asarray(ids,dtype=int);positions={int(v):i for i,v in enumerate(ids)};parent=np.arange(len(ids))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    strides=np.array([21*11*31,11*31,31,1])
    multi=np.array(np.unravel_index(ids,c.SHAPE)).T
    for i,node in enumerate(ids):
        for ax,stride in enumerate(strides):
            if multi[i,ax]+1<c.SHAPE[ax] and int(node+stride) in positions:
                a=find(i);b=find(positions[int(node+stride)])
                if a!=b:parent[b]=a
    groups={}
    for i,node in enumerate(ids):groups.setdefault(find(i),[]).append(node)
    return [np.array(v) for _,v in sorted(groups.items(),key=lambda v:min(v[1]))]
def run():
    configs=pd.read_csv(c.OUT/'CONFIGURATION_RESULTS.csv');grid=np.load(c.OUT/'GRID_STATES.npy');rows=[];summaries=[]
    for g in range(6):
        archive=np.load(c.OUT/f'accepted/G{g+1:02}.npz');memo={}
        for _,row in configs[configs.geometry==f'G{g+1:02}'].iterrows():
            ids=archive[f'{row.method}_{int(row.config)}_ids'];key=ids.tobytes()
            if key not in memo:memo[key]=components(ids)
            groups=memo[key];base=dict(geometry=row.geometry,method=row.method,config=int(row.config),sigma=row.sigma,source=row.source,replicate=int(row.replicate))
            summaries.append(dict(base,components=len(groups),accepted_states=len(ids),range_projection_components=int(row.range_branches),definition='axis-adjacent original finite grid; disconnected range projection reported separately'))
            for bi,group in enumerate(groups):
                states=grid[group];e=c.errors(states,c.GEOMS[g]);norm=np.linalg.norm(e/c.STEP,axis=1);nearest=int(norm.argmin())
                record=dict(base,branch=bi,count=len(group),truth_in_branch=bool((e==0).all(axis=1).any()),closest_joint_distance_grid_steps=float(norm.min()),furthest_joint_distance_grid_steps=float(norm.max()),nearest_state_id=int(group[nearest]))
                for j,ax in enumerate(c.AXES):
                    record['min_'+ax]=float(states[:,j].min());record['max_'+ax]=float(states[:,j].max());record['span_'+ax]=float(np.ptp(states[:,j]))
                    record['nearest_joint_point_error_'+ax]=float(e[nearest,j]);record['worst_branch_error_'+ax]=float(e[:,j].max())
                rows.append(record)
    pd.DataFrame(rows).to_csv(c.OUT/'BRANCH_COMPONENT_DIAGNOSTICS.csv',index=False,float_format='%.17g')
    pd.DataFrame(summaries).to_csv(c.OUT/'JOINT_BRANCH_COUNTS.csv',index=False)
    c.write_json(c.OUT/'BRANCH_DIAGNOSTIC_SCOPE.json',dict(rows=len(rows),configuration_rows=len(summaries),distance_definition='Euclidean parameter error divided by original grid step, evaluation-only. Closest joint state per branch is one state; componentwise minima are not combined into an oracle state.',topology='4D axis-adjacent finite-grid components, not continuous physical modes; original restricted angle domain has no +/-180 seam. Range intervals/gaps separately recorded. Empty configurations retained with zero branches.'))
if __name__=='__main__':run()
