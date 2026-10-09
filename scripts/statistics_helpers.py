"""Paired statistics for saved experimental records."""

from collections import defaultdict

from statistics import median

def key(row, names):
    return tuple(row[k] for k in names.split())

def group(rows, names):
    result = defaultdict(list)
    for row in rows:
        result[key(row, names)].append(row)
    return result

def quantile(values, fraction):
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = int(position)
    return ordered[lower] + (ordered[min(lower + 1, len(ordered) - 1)] - ordered[lower]) * (position - lower)

def scaling(rows):
    result=[]
    for _,values in sorted(group([r for r in rows if r['group']=='random'],'n metric').items()):
        assert len(values)==10
        first=values[0];r={k:first[k] for k in ['n','metric','steps']}
        r.update(independent_inputs=10,repeats_per_input=5,baseline=first['baseline'])
        for f in ['baseline_us','proposed_us','paired_speedup','baseline_total_ms','proposed_total_ms','total_speedup']:
            vs=[v[f] for v in values];r.update({f:median(vs),f+'_q1':quantile(vs,.25),f+'_q3':quantile(vs,.75)})
        result.append(r)
    return result
