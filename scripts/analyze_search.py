#!/usr/bin/env python3
"""Verify raw records/accepted states and generate the search summaries."""
from pathlib import Path
import csv,hashlib,json,statistics
ROOT=Path(__file__).resolve().parents[1]
RAW=ROOT/'data/measurements/search'
OUT=ROOT/'results/search';OUT.mkdir(parents=True,exist_ok=True)
WEIGHTS = [bin(i).count("1") for i in range(256)]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_box(path):
    values = list(map(int, path.read_text().split()))
    assert values[:2] == [8, 8]
    box = values[2:]
    assert sorted(box) == list(range(256))
    return box


def properties(box):
    histogram = [0] * 9
    max_walsh = 0
    for mask in range(1, 256):
        bits = [WEIGHTS[mask & v] % 2 for v in box]
        coefficients = bits[:]
        walsh = [1 - 2 * b for b in bits]
        for bit in range(8):
            step = 1 << bit
            for u in range(256):
                if u & step:
                    coefficients[u] ^= coefficients[u ^ step]
            for base in range(0, 256, 2 * step):
                for offset in range(step):
                    i, j = base + offset, base + offset + step
                    walsh[i], walsh[j] = walsh[i] + walsh[j], walsh[i] - walsh[j]
        degree = max((WEIGHTS[u] for u, c in enumerate(coefficients) if c), default=0)
        histogram[degree] += 1
        max_walsh = max(max_walsh, max(map(abs, walsh)))
    du, witness = 0, None
    for a in range(1, 256):
        counts = [0] * 256
        # Count each unordered input pair once, with its two solutions.
        for x in range(256):
            if x < (x ^ a):
                counts[box[x] ^ box[x ^ a]] += 2
        assert sum(counts) == 256 and all(c % 2 == 0 for c in counts)
        peak = max(counts)
        if peak > du:
            du, witness = peak, [a, counts.index(peak)]
    return {
        "nonlinearity": 128 - max_walsh // 2,
        "differential_uniformity": du,
        "du_witness_a_b": witness,
        "minimum_degree": next(i for i, c in enumerate(histogram) if c),
        "spectrum_sum": sum(i * c for i, c in enumerate(histogram)),
        "degree_histogram": histogram,
    }

protocol=json.loads((RAW/'protocol.json').read_text())
complete=json.loads((RAW/'COMPLETE.json').read_text());assert complete['passed']
states={};groups={};verified_records=[]
for job in protocol['jobs']:
 path=RAW/(job['id']+'.json');record=json.loads(path.read_text());r=record['result']
 assert record['returncode']==0 and record['binary_sha256']==complete['binary_sha256']
 inp=ROOT/'data/literature/inputs'/(job['input']+'.lut');assert sha(inp)==record['input_sha256']
 box=read_box(inp);path_states=[]
 def check_state(box,minimum,spectrum_sum):
  key=hashlib.sha256(bytes(box)).hexdigest()
  if key not in states:states[key]={'box':box[:],**properties(box)}
  props=states[key];assert props['nonlinearity']==104 and props['differential_uniformity']==8
  assert props['minimum_degree']==minimum and props['spectrum_sum']==spectrum_sum
  return key
 path_states.append(check_state(box,r['initial_min'],r['initial_sum']))
 for a in r['accepted']:
  box[a['p']],box[a['q']]=box[a['q']],box[a['p']]
  path_states.append(check_state(box,a['minimum'],a['spectrum_sum']))
 assert box==r['final_box'] and states[path_states[-1]]['degree_histogram']==r['final_histogram'][:9]
 assert r['target_reached']==(r['final_min']==7)
 assert r['total_ns']==r['context_ns']+r['init_ns']+r['kernel_ns']
 assert r['init_ns']==r['degree_init_ns']+r['constraints_init_ns']
 assert r['attempts']==r['degree_rejects']+r['nl_rejects']+r['du_rejects']+len(r['accepted'])
 verified_records.append({'id':job['id'],'sha256':sha(path),'accepted_state_hashes':path_states})
 groups.setdefault((job['input'],job['method']),[]).append(record)
summary=[]
for (name,method),records in sorted(groups.items()):
 formal=[r for r in records if r['kind']=='search'];assert len(formal)==5
 assert {r['repetition'] for r in formal}==set(range(5));formal.sort(key=lambda r:r['repetition'])
 ref=formal[0]['result']
 stable=['attempts','nl_calls','du_calls','degree_rejects','nl_rejects','du_rejects','accepted','final_box','final_histogram','candidate_digest','termination']
 assert all(all(r['result'][k]==ref[k] for k in stable) for r in records)
 diag=next(r['result'] for r in records if r['kind']=='diagnostic')
 row={'input':name,'method':method,'objective':ref['objective'],
      **{k:ref[k] for k in ['initial_min','initial_sum','final_min','final_sum','final_nl','final_du','attempts','nl_calls','du_calls','degree_rejects','nl_rejects','du_rejects','target_reached','termination']},
      'accepted_count':len(ref['accepted']),'accepted':ref['accepted']}
 for field in ['total','kernel','context','degree_init','constraints_init']:
  values=[r['result'][field+'_ns']/1e6 for r in formal]
  row['median_'+field+'_ms']=statistics.median(values);row['min_'+field+'_ms']=min(values);row['max_'+field+'_ms']=max(values)
  row['repetitions_'+field+'_ms']=values
 row['diagnostic_stage_ms']={k:diag[k+'_ns']/1e6 for k in ['degree','nl','du','rollback','other']}
 row['diagnostic_nl_du_share']=(diag['nl_ns']+diag['du_ns'])/diag['kernel_ns']
 summary.append(row)
for name in sorted({r['input'] for r in summary}):
 rows={r['method']:r for r in summary if r['input']==name};a=rows['rank-static'];b=rows['rank-shared'];p=rows['peigen']
 assert a['accepted']==b['accepted'] and a['attempts']==b['attempts']
 b['static_over_update_ratio']=a['median_total_ms']/b['median_total_ms']
 b['peigen_over_update_time_to_target_ratio']=p['median_total_ms']/b['median_total_ms'] if p['target_reached'] else None
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
fields=[k for k in summary[0] if k not in ['accepted','diagnostic_stage_ms']]
with (OUT/'search.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader()
 for r in summary:w.writerow({k:json.dumps(v) if isinstance(v,list) else v for k,v in r.items()})
(OUT/'verification.json').write_text(json.dumps({'passed':True,'raw_records':len(verified_records),'formal_records':45,'unique_states':len(states),'states':states,'records':verified_records},indent=2)+'\n')
print('PASS:',len(verified_records),'records;',len(states),'independently recomputed states')
for r in summary:print(r['input'],r['method'],r['final_min'],round(r['median_total_ms'],3),r.get('static_over_update_ratio'),r.get('peigen_over_update_time_to_target_ratio'))
