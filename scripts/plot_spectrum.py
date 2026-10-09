#!/usr/bin/env python3
"""Plot the paired degree-spectrum comparison."""
import argparse,csv,hashlib,json
from datetime import datetime,timezone
from pathlib import Path
from statistics import median
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
BLUE,RED,GRAY,GREEN='#1565B5','#B74D48','#646D78','#3D8061'
plt.rcParams.update({'font.family':'DejaVu Serif','mathtext.fontset':'stix','font.size':9,
    'axes.titlesize':10,'axes.labelsize':9.5,'xtick.labelsize':8,'ytick.labelsize':8,
    'legend.fontsize':7.3,'legend.frameon':False,'axes.spines.top':False,
    'axes.spines.right':False,'axes.linewidth':.8,'lines.linewidth':1.8,
    'lines.markersize':3.5,'grid.linewidth':.45,'pdf.fonttype':42,'ps.fonttype':42,
    'svg.fonttype':'none','svg.hashsalt':'fast-sbox-degree-spectrum','figure.facecolor':'white',
    'savefig.facecolor':'white','axes.axisbelow':True})

def read(path):
    with path.open() as f:return list(csv.DictReader(f))

def band(ax,rows,key,color,label,marker='o',style='-'):
    x=np.array([int(r['n']) for r in rows]);y=np.array([float(r[key]) for r in rows])
    ax.plot(x,y,color=color,label=label,marker=marker,linestyle=style,zorder=3)
    ax.fill_between(x,[float(r[key+'_q1']) for r in rows],[float(r[key+'_q3']) for r in rows],
                    color=color,alpha=.6,linewidth=0)

def spectrum(rows):
    fig,(a,b)=plt.subplots(1,2,figsize=(7.4,3.55))
    fig.subplots_adjust(left=.105,right=.985,bottom=.19,top=.90,wspace=.36)
    band(a,rows,'static_us',GREEN,'Static spectrum',marker='D')
    band(a,rows,'rank_us',BLUE,'Spectrum update')
    for name,label,marker in [('peigen','PEIGEN (3–8 bits)','s'),('bitwise','Bitwise (9–19 bits)','^')]:
        band(a,[r for r in rows if r['baseline']==name],'baseline_us',RED,label,marker,'--')
    a.set_ylabel('Time per evaluation (µs)');a.legend(loc='upper left')
    band(b,rows,'baseline_over_static',GREEN,'Static spectrum',marker='D')
    band(b,rows,'baseline_over_rank',BLUE,'Spectrum update')
    b.axhline(1,color=GRAY,linestyle=':',linewidth=.85)
    b.set_ylabel('Speedup over PEIGEN / bitwise (×)');b.legend(loc='upper left')
    for ax,title in [(a,'(a) Execution time ↓'),(b,'(b) Speedup ↑')]:
        ax.set_title(title,loc='left',pad=8);ax.set_yscale('log')
        ax.set_xlim(2.7,19.5);ax.set_xticks(range(3,20));ax.set_xlabel(r'S-box size $n$ ($n\times n$)')
        ax.axvline(8.5,color='#A6ADB5',linewidth=.75,linestyle=':');ax.grid(axis='y',alpha=.30)
    return fig

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--results',type=Path,default=ROOT/'results/spectrum')
    parser.add_argument('--output',type=Path,default=ROOT/'figures');args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True);sources={}
    plots=[('scaling.csv','03_degree_spectrum_update',spectrum)]
    for file,name,build in plots:
        path=args.results/file
        if not path.exists():continue
        sources[file]=hashlib.sha256(path.read_bytes()).hexdigest();fig=build(read(path));fig.canvas.draw()
        for ax in fig.axes:
            bbox=ax.get_tightbbox(fig.canvas.get_renderer())
            assert bbox.x0>=-1 and bbox.y0>=-1 and bbox.x1<=fig.bbox.x1+1 and bbox.y1<=fig.bbox.y1+1,(name,bbox)
        for ext in ['pdf','svg','png']:
            meta={'Creator':'Matplotlib; scripts/plot_spectrum.py'} if ext!='png' else {'Software':'Matplotlib; scripts/plot_spectrum.py'}
            if ext=='pdf':meta.update(CreationDate=datetime(2026,10,7,tzinfo=timezone.utc),ModDate=datetime(2026,10,7,tzinfo=timezone.utc))
            if ext=='svg':meta['Date']='2026-10-07'
            fig.savefig(args.output/(name+'.'+ext),dpi=300,metadata=meta)
        plt.close(fig)
    (args.output/'spectrum-manifest.json').write_text(json.dumps({'source_sha256':sources,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')

if __name__=='__main__':main()
