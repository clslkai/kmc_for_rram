import argparse,csv,json,os
from pathlib import Path
from dataclasses import asdict
import numpy as np
from .model import AtomConfig,AtomSimulator,SOURCE

def export_xyz(path,cr,vac):
    """Only real occupied atoms exported; missing oxygen sites in separate NPZ."""
    xyz=np.concatenate([cr.hf,cr.oxygen[~vac]])
    symbols=['Hf']*len(cr.hf)+['O']*int((~vac).sum())
    cell=' '.join(f'{x:.8f}' for x in cr.cell.ravel())
    with path.open('w') as f:
        f.write(f'{len(xyz)}\nLattice="{cell}" Properties=species:S:1:pos:R:3 pbc="T T F"\n')
        for s,p in zip(symbols,xyz): f.write(s+' '+' '.join(f'{x:.8f}' for x in p)+'\n')

def plot(out,sim,records,states):
    os.environ.setdefault('MPLCONFIGDIR','/tmp/atomistic-mpl')
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Line3DCollection,Poly3DCollection
    cr=sim.crystal
    # Draw physical, within-cell Hf-O bonds, not wrap-spanning misleading long lines.
    d=np.linalg.norm(cr.hf[:,None,:]-cr.oxygen[None,:,:],axis=2)
    bi,bj=np.where(d<sim.cfg.bond_cutoff_A)
    fig=plt.figure(figsize=(13,8))
    for k,index,title in [(1,0,'Initial single crystal'),(2,-1,'After KMC events')]:
        ax=fig.add_subplot(2,2,k,projection='3d'); v=states[index]
        segments=[np.array([cr.hf[i],cr.oxygen[j]])/10 for i,j in zip(bi,bj) if not v[j]]
        ax.add_collection3d(Line3DCollection(segments,colors='#94a3b8',linewidths=.35,alpha=.22))
        ax.scatter(*(cr.hf/10).T,s=20,c='#4b76af',alpha=.6,label='Hf (fixed)')
        ax.scatter(*(cr.oxygen[~v]/10).T,s=7,c='#d14b4b',alpha=.35,label='O (occupied)')
        ax.scatter(*(cr.oxygen[v]/10).T,s=48,facecolors='none',edgecolors='#f59e0b',linewidths=1.7,label='O vacancy (missing atom)')
        a,b,c=cr.cell/10
        for off,color,label in [(np.zeros(3),'#64748b','blocking bottom'),(c,'#7c3aed','top oxygen reservoir')]:
            verts=np.array([off,off+a,off+a+b,off+b])
            ax.add_collection3d(Poly3DCollection([verts],facecolors=color,alpha=.15))
        ax.set(xlabel='x (nm)',ylabel='y (nm)',zlabel='z (nm)',title=title)
        ax.set_box_aspect((np.ptp(cr.oxygen[:,0]),np.ptp(cr.oxygen[:,1]),cr.height));ax.view_init(18,-60)
        ax.legend(fontsize=7,loc='upper left')
    time=np.array([r['time_s'] for r in records])*1e6
    ax=fig.add_subplot(2,2,3);ax.step(time,[r['vacancies'] for r in records],where='post',label='Oxide vacancies')
    ax.step(time,[r['reservoir_oxygen'] for r in records],'--',where='post',label='Reservoir O (same count)')
    ax.set(xlabel='Time (microseconds)',ylabel='Count');ax.legend()
    ax=fig.add_subplot(2,2,4);ax.semilogy(time,[r['conductance_S'] for r in records]);ax.set(xlabel='Time (microseconds)',ylabel='Effective graph conductance (S)')
    fig.suptitle('Crystallographic m-HfO2 oxygen-site KMC | real structure, illustrative kinetics',fontsize=13)
    fig.tight_layout();fig.savefig(out/'atomistic_overview.png',dpi=160);plt.close(fig)

def main():
    p=argparse.ArgumentParser(description='Monoclinic HfO2 oxygen-site lattice KMC, illustrative kinetics')
    p.add_argument('--repeats',type=int,nargs=3,default=None);p.add_argument('--seed',type=int,default=None)
    p.add_argument('--fraction',type=float,default=.03);p.add_argument('--temperature',type=float,default=None)
    p.add_argument('--pulse',type=float,nargs=2,action='append');p.add_argument('--max-events',type=int,default=1000)
    p.add_argument('--sample-every',type=int,default=10);p.add_argument('--config',type=Path)
    p.add_argument('--output',type=Path,default=Path('results/atomistic_demo'));p.add_argument('--no-plot',action='store_true')
    args=p.parse_args()
    if args.max_events<1 or args.sample_every<1:p.error('event limits must be positive')
    data=json.loads(args.config.read_text()) if args.config else {}
    if 'repeats' in data:data['repeats']=tuple(data['repeats'])
    if args.repeats is not None:data['repeats']=tuple(args.repeats)
    if args.seed is not None:data['seed']=args.seed
    if args.temperature is not None:data['temperature_K']=args.temperature
    cfg=AtomConfig(**data);sim=AtomSimulator(cfg);sim.initialize(args.fraction)
    pulses=args.pulse or [(1.,1e-6),(-1.,1e-6)]
    if any(not np.isfinite(v) or not np.isfinite(t) or t<=0 for v,t in pulses):p.error('invalid pulse')
    out=args.output;out.mkdir(parents=True,exist_ok=True)
    records=[];states=[];potentials=[];logs=[]
    def capture(v):records.append(sim.snapshot(v));states.append(sim.vacancy.copy());potentials.append(sim.phi.copy())
    export_xyz(out/'initial.extxyz',sim.crystal,sim.vacancy)
    end=0.;status='completed_waveform'
    for v,duration in pulses:
        capture(v);end+=duration
        while sim.time_s<end and sim.event_count<args.max_events:
            e=sim.step(v,end)
            if e:
                logs.append(dict(time_s=sim.time_s,kind=e.kind,source=e.source,target=e.target,delta_eV=e.delta_eV,rate_Hz=e.rate_Hz))
                if sim.event_count%args.sample_every==0:capture(v)
        capture(v)
        if sim.time_s<end:status='event_limit_before_waveform_end';break
    for name,rows in [('trace.csv',records),('events.csv',logs)]:
        with (out/name).open('w',newline='') as f:
            if rows:w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    export_xyz(out/'final.extxyz',sim.crystal,sim.vacancy)
    np.savez_compressed(out/'states.npz',time_s=[r['time_s'] for r in records],vacancy=states,potential_V=potentials,oxygen_A=sim.crystal.oxygen,hf_A=sim.crystal.hf,oxygen_kind=sim.crystal.kind,cell_A=sim.crystal.cell,edges=sim.crystal.edges)
    summary=dict(config=asdict(cfg),source=SOURCE,structure='P21/c monoclinic HfO2, Pathak combined refinement; idealized Hf occupancy 1',
                 status=status,event_counts=sim.counts,final=records[-1],oxygen_sites=sim.n,hf_atoms=len(sim.crystal.hf),bulk_O_coordination=sim.crystal.bulk_coordination,
                 pulses=pulses,warning='Real crystallographic sites; fixed lattice, isothermal illustrative rates, geometric candidate migration graph, uncalibrated effective transport. Not a validated SET/RESET simulation.')
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    from .viewer import write_viewer
    write_viewer(out,sim,records,states)
    if not args.no_plot:plot(out,sim,records,states)
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
