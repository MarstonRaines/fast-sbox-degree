#!/usr/bin/env python3
"""Validate degree-spectrum, finite-prefix and operation-count measurements."""
import argparse,csv,hashlib,json,math
from collections import defaultdict
from pathlib import Path
from statistics import median

ROOT=Path(__file__).resolve().parents[1]

def quantile(values,p):
    values=sorted(values);x=(len(values)-1)*p;lo=math.floor(x);hi=math.ceil(x)
    return values[lo]+(values[hi]-values[lo])*(x-lo)

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def write_csv(path,rows):
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--measurements',type=Path,default=ROOT/'data/measurements/spectrum')
    parser.add_argument('--output',type=Path,default=ROOT/'results/spectrum');args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    formal=args.measurements/'formal'
    assert json.loads((formal/'COMPLETE.json').read_text())['passed']
    protocol=json.loads((formal/'protocol.json').read_text())
    grouped=defaultdict(list);input_signatures={};raw=[];sources={}
    for job in protocol['jobs']:
        path=formal/f"{job['input']}--{job['repetition']}--{job['method']}.json"
        record=json.loads(path.read_text());assert record['returncode']==0
        r=record['result'];assert r['steps']==job['steps'] and r['n']==job['n']
        assert record['input_sha256']==digest(ROOT/f"data/inputs/{job['input']}.lut")
        assert record['swaps_sha256']==digest(ROOT/f"data/inputs/{job['input']}.swaps")
        assert r['total_ns']==r['context_ns']+r['init_ns']+r['kernel_ns']
        assert sum(r['histogram'])==(1<<r['m'])-1
        assert not any(r['histogram'][r['n']+1:])
        assert sum(d*v for d,v in enumerate(r['histogram']))==r['final_sum']
        signature=(r['digest'],r['final_box_hash'],tuple(r['histogram'][:r['n']+1]))
        if job['input'] in input_signatures:assert input_signatures[job['input']]==signature
        else:input_signatures[job['input']]=signature
        assert r['peak_rss_bytes']>0
        grouped[(job['input'],job['method'])].append(r)
        raw.append({**job,**{k:r[k] for k in ['context_ns','init_ns','kernel_ns','total_ns','peak_rss_bytes','context_allocated_bytes','state_allocated_bytes']}})
        sources[str(path.relative_to(args.measurements))]=digest(path)
    assert len(raw)==len(protocol['jobs'])
    three_routes=True
    assert len(raw)==2705
    medians={}
    for (name,method),records in grouped.items():
        assert len(records)==5
        medians[name,method]={'input':name,'method':method,'n':records[0]['n'],'steps':records[0]['steps'],
            'update_us':median(r['kernel_ns']/r['steps']/1000 for r in records),
            'context_ms':median(r['context_ns']/1e6 for r in records),
            'init_ms':median(r['init_ns']/1e6 for r in records),'total_ms':median(r['total_ns']/1e6 for r in records),
            'peak_rss_mib':median(r['peak_rss_bytes']/2**20 for r in records),
            'state_mib':median(r['state_allocated_bytes']/2**20 for r in records),
            'context_mib':median(r['context_allocated_bytes']/2**20 for r in records)}
    inputs=[]
    for name in sorted(input_signatures):
        rank=medians[name,'rank-shared'];n=rank['n']
        baseline=medians[name,'peigen' if n<=8 else 'bitwise']
        inputs.append({'input':name,'n':n,'steps':rank['steps'],'baseline':baseline['method'],
            'baseline_us':baseline['update_us'],'rank_us':rank['update_us'],
            'baseline_over_rank':baseline['update_us']/rank['update_us'],
            'baseline_total_ms':baseline['total_ms'],'rank_total_ms':rank['total_ms'],
            'baseline_total_over_rank':baseline['total_ms']/rank['total_ms']})
        if three_routes:
            static=medians[name,'rank-static']
            inputs[-1].update(static_us=static['update_us'],
                baseline_over_static=baseline['update_us']/static['update_us'],
                static_over_rank=static['update_us']/rank['update_us'],
                static_total_ms=static['total_ms'],
                baseline_total_over_static=baseline['total_ms']/static['total_ms'],
                static_total_over_rank=static['total_ms']/rank['total_ms'])
    scaling=[]
    numeric=[k for k in inputs[0] if k not in ['input','n','steps','baseline']]
    for n in range(3,20):
        rows=[r for r in inputs if r['n']==n and r['input'].startswith('random_')];assert len(rows)==10
        row={'n':n,'steps':rows[0]['steps'],'baseline':rows[0]['baseline']}
        for key in numeric:
            values=[r[key] for r in rows];row[key]=median(values);row[key+'_q1']=quantile(values,.25);row[key+'_q3']=quantile(values,.75)
        scaling.append(row)
    write_csv(args.output/'raw-times.csv',raw);write_csv(args.output/'per-method-medians.csv',list(medians.values()))
    write_csv(args.output/'per-input.csv',inputs);write_csv(args.output/'scaling.csv',scaling)
    write_csv(args.output/'practical.csv',[r for r in inputs if not r['input'].startswith('random_')])
    overlaps=[]
    for name in sorted(input_signatures):
        if (name,'peigen') in medians and (name,'bitwise') in medians:
            p=medians[name,'peigen'];b=medians[name,'bitwise']
            overlaps.append({'input':name,'peigen_us':p['update_us'],'bitwise_us':b['update_us'],'peigen_over_bitwise':p['update_us']/b['update_us']})
    write_csv(args.output/'overlap-8bit.csv',overlaps)
    cited=ROOT/'results/round2';cited.mkdir(parents=True,exist_ok=True)
    write_csv(cited/'overlap-8bit.csv',overlaps)
    follow=args.measurements/'followups'
    if (follow/'COMPLETE.json').exists():
        jobs=json.loads((follow/'protocol.json').read_text())['jobs'];prefix=defaultdict(list);counts=[];checks={}
        for job in jobs:
            path=follow/(job['id']+'.json');record=json.loads(path.read_text());assert record['returncode']==0;r=record['result']
            for relative,expected_hash in record['input_sha256'].items():assert digest(ROOT/relative)==expected_hash
            sources[str(path.relative_to(args.measurements))]=digest(path)
            if job['kind']=='prefix':
                for point in r['checkpoints']:
                    assert point['total_ns']==r['context_ns']+r['init_ns']+point['kernel_ns']
                    key=(job['input'],point['steps'])
                    signature=(point['digest'],tuple(point['histogram']))
                    if key in checks:assert checks[key]==signature
                    else:checks[key]=signature
                    prefix[job['input'],job['method'],point['steps']].append({**point,'context_ns':r['context_ns'],'init_ns':r['init_ns'],'peak_rss_bytes':r['peak_rss_bytes']})
            if job['kind']=='counts':
                expected=0;n=r['n'];pairs=[tuple(map(int,line.split()))for line in (ROOT/f"data/inputs/{job['input']}.swaps").read_text().splitlines() if line.strip()]
                for i in range(job['steps']):
                    a,b=pairs[i%len(pairs)];expected+=(1<<(n-a.bit_count()))+(1<<(n-b.bit_count()))-(2<<(n-(a|b).bit_count()))
                assert r['coefficient_updates']==expected
                counts.append({'input':job['input'],'steps':job['steps'],**{k:r[k] for k in ['scanned_rows','basis_xors','coefficient_updates','final_sum']}})
        prefix_rows=[]
        for (name,method,steps),rows in sorted(prefix.items()):
            assert len(rows)==5
            prefix_rows.append({'input':name,'method':method,'steps':steps,'context_ms':median(r['context_ns']/1e6 for r in rows),
                'init_ms':median(r['init_ns']/1e6 for r in rows),'kernel_ms':median(r['kernel_ns']/1e6 for r in rows),
                'total_ms':median(r['total_ns']/1e6 for r in rows),'min_total_ms':min(r['total_ns']/1e6 for r in rows),
                'max_total_ms':max(r['total_ns']/1e6 for r in rows),'peak_rss_mib':median(r['peak_rss_bytes']/2**20 for r in rows)})
        write_csv(args.output/'prefix-costs.csv',prefix_rows)
        write_csv(args.output/'rank-counts.csv',counts)
    (args.output/'validation.json').write_text(json.dumps({'passed':True,'formal_records':len(raw),'independent_inputs':len(input_signatures),
        'five_repeat_method_groups':len(medians),'paired_final_histograms_and_spectrum_sum_digests_checked':True,
        'followups_complete':(follow/'COMPLETE.json').exists(),'source_sha256':sources},indent=2)+'\n')
    print('Validated',len(raw),'formal records;',len(medians),'five-repeat groups.')

if __name__=='__main__':main()
