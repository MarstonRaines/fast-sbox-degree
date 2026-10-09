"""Validate every minimum/maximum record and derive paired statistics."""
import argparse,csv,hashlib,json
from collections import defaultdict
from pathlib import Path
from statistics import median
from statistics_helpers import scaling,quantile
ROOT=Path(__file__).resolve().parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def write(path,rows):
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
 p=argparse.ArgumentParser();p.add_argument('--measurements',type=Path,default=ROOT/'data/measurements/minmax');p.add_argument('--output',type=Path,default=ROOT/'results/minmax');a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
 protocol=json.loads((a.measurements/'protocol.json').read_text());assert json.loads((a.measurements/'COMPLETE.json').read_text())['passed']
 assert len(protocol['jobs'])==3520
 grouped=defaultdict(list);signatures={};hashes={}
 for job in protocol['jobs']:
  name=job['input'];path=a.measurements/f"{name}-{job['metric']}-{job['repetition']}-{job['method']}.json";record=json.loads(path.read_text());assert record['returncode']==0
  assert all(record[k]==v for k,v in job.items())
  assert record['input_sha256']==sha(ROOT/f'data/inputs/{name}.lut') and record['swaps_sha256']==sha(ROOT/f'data/inputs/{name}.swaps')
  r=record['result'];assert r['steps']==job['steps'] and r['kernel_ns']>0
  sig=(r['final_value'],r['final_box_hash'],r['digest']);key=(name,job['metric']);assert signatures.setdefault(key,sig)==sig
  grouped[name,job['metric'],job['method']].append(record);hashes[path.name]=sha(path)
 main=[]
 for (name,metric,method),records in sorted(grouped.items()):
  assert sorted(r['repetition'] for r in records)==list(range(5));r=records[0]
  vals=[x['result']['kernel_ns']/r['steps']/1000 for x in records]
  row=dict(case=name,n=r['n'],metric=metric,method=method,steps=r['steps'],group='random' if name.startswith('random_') else 'practical',repeats=5,final_value=r['result']['final_value'],final_box_hash=r['result']['final_box_hash'],kernel_us_per_step=median(vals),repeat_q1_us=quantile(vals,.25),repeat_q3_us=quantile(vals,.75))
  row.update({key+'_ms':median(x['result'][key+'_ns']/1e6 for x in records) for key in ['context','init','kernel']});row['total_ms']=median(sum(x['result'][k+'_ns'] for k in ['context','init','kernel'])/1e6 for x in records)
  main.append(row)
 lookup={(r['case'],r['metric'],r['method']):r for r in main};paired=[]
 for own in main:
  if own['method']!='proposed':continue
  baseline='peigen' if own['n']<=8 else 'bitwise';base=lookup[own['case'],own['metric'],baseline]
  row={k:own[k] for k in ['case','group','n','metric','steps','final_value']}
  row.update(baseline_us=base['kernel_us_per_step'],proposed_us=own['kernel_us_per_step'],paired_speedup=base['kernel_us_per_step']/own['kernel_us_per_step'],baseline_total_ms=base['total_ms'],proposed_total_ms=own['total_ms'],total_speedup=base['total_ms']/own['total_ms'],baseline=baseline);paired.append(row)
 write(a.output/'per-input-results.csv',main);write(a.output/'paired-results.csv',paired);write(a.output/'practical-results.csv',[r for r in paired if r['group']=='practical']);write(a.output/'scaling-results.csv',scaling(paired))
 (a.output/'validation.json').write_text(json.dumps(dict(passed=True,records=3520,paired_input_metrics=len(paired),all_five_repeats_retained=True,source_sha256=hashes),indent=2)+'\n')
 print('Validated 3520 min/max records; 352 paired input/metric groups.')
if __name__=='__main__':main()
