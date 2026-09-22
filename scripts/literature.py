#!/usr/bin/env python3
"""Check all three literature witnesses and timings; optionally replay C++ searches."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
from statistics import median
import subprocess

ROOT=Path(__file__).resolve().parents[1]


def read(p):
    return json.loads(p.read_text())


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def lut(p):
    values=list(map(int,p.read_text().split()))
    assert values[:2]==[8,8] and sorted(values[2:])==list(range(256))
    return values[2:]


def metrics(box):
    parity=[x.bit_count()&1 for x in range(256)]
    coefficients=box.copy()
    for bit in [1,2,4,8,16,32,64,128]:
        for u in range(256):
            if u&bit:
                coefficients[u]^=coefficients[u^bit]
    order=sorted(range(256),key=lambda u:u.bit_count(),reverse=True)
    degrees=[]
    maximum_walsh=0
    for mask in range(1,256):
        degrees.append(next((u.bit_count() for u in order if parity[mask&coefficients[u]]),0))
        w=[1-2*parity[mask&x] for x in box]
        for step in [1,2,4,8,16,32,64,128]:
            for base in range(0,256,2*step):
                for j in range(step):
                    a,b=w[base+j],w[base+step+j]
                    w[base+j],w[base+step+j]=a+b,a-b
        maximum_walsh=max(maximum_walsh,max(map(abs,w)))
    du=max(max(Counter(box[x]^box[x^a] for x in range(256)).values()) for a in range(1,256))
    basis={}
    for u in range(256):
        if u.bit_count()!=7:
            continue
        v=coefficients[u]
        while v:
            pivot=v.bit_length()-1
            if pivot in basis:
                v^=basis[pivot]
            else:
                basis[pivot]=v
                break
    return {'nl':128-maximum_walsh//2,'minimum':min(degrees),'maximum':max(degrees),
            'sum':sum(degrees),'histogram':dict(Counter(degrees)),'du':du,'rank7':len(basis)}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--check',action='store_true',help='Audit archived evidence (also the default)')
    ap.add_argument('--replay',action='store_true',help='Run both C++ evaluators from the original three LUTs')
    ap.add_argument('--output',type=Path,help='New directory for replay results')
    ap.add_argument('--cpu',type=int)
    ap.add_argument('--numa-node',type=int)
    args=ap.parse_args()
    if args.numa_node is not None and args.cpu is None:
        ap.error('--numa-node requires --cpu')
    origin=read(ROOT/'data/literature/sources-and-witnesses.json')
    source=ROOT/'data/measurements/literature.jsonl'
    rows=[json.loads(x) for x in source.read_text().splitlines()]
    audit=read(ROOT/'results/literature.json')
    assert len(rows)==30 and sha(source)==audit['records_sha256']
    assert {(r['case'],r['method'],r['repeat']) for r in rows}=={
        (c,m,k) for c in origin for m in ['peigen','proposed'] for k in range(5)}
    checked={}
    for case,info in origin.items():
        for key in ['original','final']:
            assert sha(ROOT/info[key+'_lut'])==info[key+'_sha256']
        rr=[r for r in rows if r['case']==case]
        fields=['attempts','nl_calls','degree_rejects','nl_rejects','candidate_digest','accepted','final_box','final_histogram']
        assert all(r['verified'] and r['total_ns']>0 and all(r[k]==rr[0][k] for k in fields) for r in rr)
        box=lut(ROOT/info['original_lut'])
        initial=metrics(box)
        assert initial['nl']==104 and initial['du']==8 and initial['minimum']==6
        assert initial['sum']==info['input_metrics']['degree_sum']
        states=[initial]
        for p,q,attempt,degree,total in rr[0]['accepted']:
            box[p],box[q]=box[q],box[p]
            current=metrics(box)
            assert current['nl']==104 and current['du']==8
            assert (current['minimum'],current['sum'])==(degree,total)
            assert (degree,total)>(states[-1]['minimum'],states[-1]['sum'])
            states.append(current)
        assert box==lut(ROOT/info['final_lut'])==rr[0]['final_box']
        assert states[-1]['histogram']=={7:255}
        summary=next(s for s in audit['summary'] if s['case']==case)
        for method in ['peigen','proposed']:
            values=[r['total_ns']/1e6 for r in rr if r['method']==method]
            assert sorted(values)==sorted(summary[method]['all_ms'])
            assert median(values)==summary[method]['median_ms']
        assert summary['speedup']==summary['peigen']['median_ms']/summary['proposed']['median_ms']
        checked[case]=states
    assert [x['rank7'] for x in checked['freyre2020-s8']]==[5,6,7,8]
    if args.replay:
        out=args.output or ROOT/'runs'/('literature-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
        out.mkdir(parents=True,exist_ok=False)
        prefix=[]
        if args.cpu is not None:
            prefix=['numactl','--physcpubind='+str(args.cpu)]
            if args.numa_node is not None:
                prefix.append('--membind='+str(args.numa_node))
        env={**os.environ,'OMP_NUM_THREADS':'1','OMP_DYNAMIC':'false'}
        for case,info in origin.items():
            expected=next(r for r in rows if r['case']==case)
            for method in ['peigen','proposed']:
                command=[*prefix,str(ROOT/'build/warm_cpu'),'50',str(ROOT/'build/postprocess'),
                         'bench',method,str(ROOT/info['original_lut'])]
                proc=subprocess.run(command,capture_output=True,text=True,env=env)
                proc.check_returncode()
                assert not proc.stderr
                row=json.loads(proc.stdout)
                assert row['verified'] and all(row[k]==expected[k] for k in fields)
                (out/(case+'_'+method+'.json')).write_text(proc.stdout)
        (out/'complete.json').write_text(json.dumps({'status':'pass','runs':6})+'\n')
    print(json.dumps({'status':'pass','records':30,'witnesses':checked,'replayed':args.replay},indent=2))


if __name__=='__main__':
    main()
