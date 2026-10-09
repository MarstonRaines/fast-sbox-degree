#!/usr/bin/env python3
"""Figure 4: time to termination, explicitly distinguishing failure to reach 7."""
from pathlib import Path
import hashlib,json
from datetime import datetime,timezone
import numpy as np
from matplotlib.collections import LineCollection
from plot_spectrum import plt,BLUE,RED,GREEN,GRAY
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/search/summary.json';OUT=ROOT/'figures';OUT.mkdir(parents=True,exist_ok=True)
rows=json.loads(SOURCE.read_text());fig,axes=plt.subplots(1,3,figsize=(7.4,3.55))
failed_bars=[]
fig.subplots_adjust(left=.08,right=.985,bottom=.21,top=.86,wspace=.48)
cases=[('freyre2020-s8','(a) Freyre S8'),('kuznetsov2023-hill-2','(b) Kuznetsov S-box2'),('kuznetsov2023-hill-3','(c) Kuznetsov S-box3')]
for ax,(name,title) in zip(axes,cases):
 case={r['method']:r for r in rows if r['input']==name};maximum=max(r['max_total_ms'] for r in case.values())
 for x,method,color in [(0,'peigen',RED),(1,'rank-static',GREEN),(2,'rank-shared',BLUE)]:
  r=case[method];vals=r['repetitions_total_ms'];y=r['median_total_ms'];ok=r['target_reached']
  bar=ax.bar(x,y,width=.58,color=color,alpha=.90,edgecolor=color,linewidth=.8)[0]
  if not ok:failed_bars.append((ax,bar))
  ax.errorbar(x,y,yerr=[[y-min(vals)],[max(vals)-y]],fmt='none',color=GRAY,capsize=3)
  ax.scatter(x+np.linspace(-.12,.12,5),vals,s=12,color='white',edgecolor=color,linewidth=.65,zorder=3)
  ax.text(x,max(vals)+maximum*.045,f'{y:.2f}\n'+('min. 7' if ok else 'min. 6*'),ha='center',fontsize=7,color=color,linespacing=1.4)
 ax.set_ylim(0,maximum*1.29);ax.set_xticks([0,1,2],['PEIGEN\nmin. only','Static\nspectrum','Spectrum\nupdate'],fontsize=7)
 ax.set_title(title,fontsize=9,pad=10);ax.set_ylabel('Time to termination (ms)');ax.grid(axis='y',alpha=.3)
fig.canvas.draw()
# Opaque vector strokes keep the failure hatch visible in PDF and PNG alike.
# Figure coordinates preserve a 45-degree angle and 11-point spacing at any DPI.
for ax,bar in failed_bars:
 box=bar.get_window_extent().transformed(fig.transFigure.inverted())
 rise=(box.x1-box.x0)*fig.get_figwidth()/fig.get_figheight()
 spacing=11/72/fig.get_figheight()
 segments=[[(box.x0,y),(box.x1,y+rise)] for y in np.arange(box.y0-rise,box.y1,spacing)]
 stripes=LineCollection(segments,colors=['#963A35'],linewidths=.9,transform=fig.transFigure,zorder=1.5)
 ax.add_collection(stripes,autolim=False);stripes.set_clip_path(bar)
fig.canvas.draw()
for ax in fig.axes:
 bbox=ax.get_tightbbox(fig.canvas.get_renderer());assert bbox.x0>=-1 and bbox.y0>=-1 and bbox.x1<=fig.bbox.x1+1 and bbox.y1<=fig.bbox.y1+1
for ext in ['pdf','svg','png']:
 meta={'Creator':'Matplotlib; scripts/plot_search.py'} if ext!='png' else {'Software':'Matplotlib; scripts/plot_search.py'}
 if ext=='pdf':meta.update(CreationDate=datetime(2026,10,8,tzinfo=timezone.utc),ModDate=datetime(2026,10,8,tzinfo=timezone.utc))
 if ext=='svg':meta['Date']='2026-10-08'
 fig.savefig(OUT/('04_search_time.'+ext),dpi=300,metadata=meta)
plt.close(fig)
(OUT/'search-manifest.json').write_text(json.dumps({'source':str(SOURCE.relative_to(ROOT)),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
