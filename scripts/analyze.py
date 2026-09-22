#!/usr/bin/env python3
"""Audit published records and recompute every result table using the Python standard library."""
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
CHECK = False
TABLES = {}


@lru_cache(maxsize=None)
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def emit(name, rows, keys):
    path = ROOT / 'results' / name
    if CHECK:
        with path.open() as stream:
            reader = csv.DictReader(stream)
            old = list(reader)
        fields = reader.fieldnames
        index = {tuple(str(v) for v in key(r, keys)): r for r in rows}
        assert len(index) == len(rows) == len(old), name
        for saved in old:
            new = index[key(saved, keys)]
            for field in fields:
                a, b = saved[field], new[field]
                if str(a) == str(b):
                    continue
                try:
                    assert math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-12)
                except (ValueError, TypeError, AssertionError):
                    raise AssertionError((name, key(saved, keys), field, a, b)) from None
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
    TABLES[name] = len(rows)


def emit_json(name, value):
    path = ROOT / 'results' / name
    if CHECK:
        assert json.loads(path.read_text()) == value, name
    else:
        path.write_text(json.dumps(value, indent=2) + '\n')


def load_audit():
    rows = [json.loads(line) for line in (ROOT / 'data/measurements/serial.jsonl').read_text().splitlines()]
    tasks = json.loads((ROOT / 'protocol/serial-tasks.json').read_text())
    task_map = {(t['group'], t['id']): t for t in tasks}
    assert len(rows) == len(task_map) == len({(r['group'], r['id']) for r in rows})
    outputs, signatures = {}, {}
    for r in rows:
        task = task_map[r['group'], r['id']]
        assert all(r[k] == v for k, v in task.items() if k in r), task['id']
        assert sha(ROOT / task['input']) == r['input_sha256']
        if task['swaps'] != '-':
            assert sha(ROOT / task['swaps']) == r['swaps_sha256']
        assert r['output_sha256'] == task['expected_output_sha256']
        token = r['output_sha256']
        if token not in outputs:
            path = ROOT / 'data/outputs' / (token + '.json')
            assert sha(path) == token
            outputs[token] = json.loads(path.read_text())
            box = outputs[token]['final_box']
            assert sorted(box) == list(range(1 << r['n']))
            value = 14695981039346656037
            for x in box:
                value = ((value ^ x) * 1099511628211) & ((1 << 64) - 1)
            assert value == r['final_box_hash']
        result = outputs[token]
        assert result['final_histogram'] == r['final_histogram']
        assert len(result['accepted_swaps']) == r['accepted']
        assert sum(r['final_histogram']) == (1 << r['m']) - 1
        assert r['total_ns'] == r['context_ns'] + r['init_ns'] + r['kernel_ns']
        assert r['kernel_ns'] > 0 and r['steps'] > 0
        assert r['instrumented'] == (r['group'] == 'counts')
        workload = ('A' if r['group'] == 'counts' else r['group'], *key(r, 'case metric profile loops'))
        signature = key(r, 'steps initial_value final_value candidate_digest final_box_hash accepted output_sha256')
        assert signatures.setdefault(workload, signature) == signature, workload
    for values in group(rows, 'group case metric profile method binary threads').values():
        repeats = 1 if values[0]['group'] == 'counts' else 5
        assert sorted(v['repeat'] for v in values) == list(range(repeats))
    return rows


def serial_summary(rows, reference=False):
    result = []
    for _, values in sorted(group(rows, 'case metric profile method').items()):
        first = values[0]
        record = {k: first[k] for k in ['case', 'n', 'metric', 'method', 'steps']}
        record['group'] = first['input_group']
        if reference:
            record.update(m=first['m'], profile=first['profile'])
        else:
            record.update(repeats=len(values), **{k: first[k] for k in ['initial_value', 'final_value', 'accepted', 'final_box_hash']})
            record['final_histogram'] = json.dumps(first['final_histogram'])
        times = [r['kernel_ns'] / r['steps'] / 1000 for r in values]
        record.update(kernel_us_per_step=median(times), repeat_q1_us=quantile(times, .25), repeat_q3_us=quantile(times, .75))
        timers = ['context', 'init', 'kernel', 'total'] + (['input_conversion', 'state_init'] if reference else [])
        record.update({name + '_ms': median(r[name + '_ns'] / 1e6 for r in values) for name in timers})
        result.append(record)
    return result


def serial_pairs(rows, reference=False):
    lookup = {key(r, 'case metric') + (r.get('profile'), r['method']): r for r in rows}
    result = []
    for own in rows:
        if own['method'] != ('reference-proposed' if reference else 'proposed'):
            continue
        baseline = 'reference' if reference else 'peigen' if own['n'] <= 8 else 'bitwise'
        base = lookup[key(own, 'case metric') + (own.get('profile'), baseline)]
        fields = ['case', 'group', 'n', 'metric', 'steps'] + (['m', 'profile'] if reference else ['accepted', 'initial_value', 'final_value'])
        record = {k: own[k] for k in fields}
        record.update(baseline_us=base['kernel_us_per_step'], proposed_us=own['kernel_us_per_step'],
                      paired_speedup=base['kernel_us_per_step'] / own['kernel_us_per_step'],
                      baseline_total_ms=base['total_ms'], proposed_total_ms=own['total_ms'],
                      total_speedup=base['total_ms'] / own['total_ms'])
        if reference:
            factor = {'min': 2 ** own['m'] / own['m'], 'max': own['n'], 'spectrum': own['n'] + own['m']}[own['metric']]
            record.update(reference_growth_factor=factor, speedup_over_reference_factor=record['paired_speedup'] / factor)
        else:
            record['baseline'] = baseline
        result.append(record)
    return result


def serial_scaling(rows, reference=False):
    result = []
    keys = 'n metric profile' if reference else 'n metric'
    for _, values in sorted(group([r for r in rows if r['group'] == 'random'], keys).items()):
        assert len(values) == 10
        first = values[0]
        record = {k: first[k] for k in keys.split() + ['steps']}
        record.update(independent_inputs=10, repeats_per_input=5)
        record.update({('reference_growth_factor' if reference else 'baseline'): first['reference_growth_factor' if reference else 'baseline']})
        fields = ['baseline_us', 'proposed_us', 'paired_speedup', 'baseline_total_ms', 'proposed_total_ms', 'total_speedup']
        if reference:
            fields.append('speedup_over_reference_factor')
        for field in fields:
            points = [v[field] for v in values]
            record.update({field: median(points), field + '_q1': quantile(points, .25), field + '_q3': quantile(points, .75)})
        result.append(record)
    return result


def serial_tables(rows):
    ref = serial_summary([r for r in rows if r['group'] == 'A'], True)
    pairs = serial_pairs(ref, True)
    emit('serial/reference-per-input.csv', ref, 'case metric profile method')
    emit('serial/reference-paired.csv', pairs, 'case metric profile')
    emit('serial/reference-practical.csv', [r for r in pairs if r['group'] == 'practical'], 'case metric profile')
    emit('serial/reference-scaling.csv', serial_scaling(pairs, True), 'n metric profile')
    counted = [r for r in rows if r['group'] == 'counts']
    counts = [{**{k: r[k] for k in ['case', 'n', 'm', 'metric', 'profile', 'method', 'steps']},
               'group': r['input_group'], **r['counts']} for r in counted]
    emit('serial/operation-counts.csv', counts, 'case metric profile method')
    pair_map = {key(r, 'case metric profile'): r for r in pairs}
    mechanisms = []
    for k, values in sorted(group(counted, 'case metric profile').items()):
        by_method = {v['method']: v for v in values}
        own, base = by_method['reference-proposed'], by_method['reference']
        a, b = own['counts'], base['counts']
        fields = ['scalar_mobius_xor', 'scalar_component_xor', 'scalar_elimination_xor', 'scalar_coefficient_updates']
        ox, bx = sum(a[f] for f in fields), sum(b[f] for f in fields)
        guards = a['scalar_guard_kept'] + a['scalar_guard_fallback']
        timing = pair_map[k]
        mechanisms.append({**{f: timing[f] for f in ['case', 'group', 'n', 'm', 'metric', 'profile', 'steps']},
            'baseline_anf_xor_per_step': bx / own['steps'], 'proposed_anf_xor_per_step': ox / own['steps'],
            'logical_anf_xor_ratio': bx / ox, 'baseline_degree_tests_per_step': b['scalar_degree_checks'] / own['steps'],
            'proposed_degree_tests_per_step': a['scalar_degree_checks'] / own['steps'],
            'mean_delta_fraction': a['delta_coefficients'] / own['steps'] / (1 << own['n']) if own['metric'] != 'min' else '',
            'conditional_guard_keep_fraction': a['scalar_guard_kept'] / guards if guards else '',
            'actual_degree_drops': a['scalar_degree_drops'], 'measured_time_speedup': timing['paired_speedup'],
            'reference_growth_factor': timing['reference_growth_factor']})
    emit('serial/mechanism-per-input.csv', mechanisms, 'case metric profile')
    main = serial_summary([r for r in rows if r['group'] == 'B'])
    search = serial_summary([r for r in rows if r['group'] == 'search'])
    optimized_pairs = serial_pairs(main)
    emit('serial/optimized/per-input-results.csv', main, 'case metric method')
    emit('serial/optimized/search-per-method.csv', search, 'case metric method')
    emit('serial/optimized/search-results.csv', serial_pairs(search), 'case metric')
    emit('serial/optimized/practical-results.csv', [r for r in optimized_pairs if r['group'] == 'practical'], 'case metric')
    emit('serial/optimized/scaling-results.csv', serial_scaling(optimized_pairs), 'n metric')
    initial = []
    for r in main:
        if r['group'] == 'practical' and r['method'] == 'proposed' and r['metric'] == 'min':
            h = json.loads(r['final_histogram'])
            degrees = [i for i, count in enumerate(h) if count]
            initial.append(dict(case=r['case'], n=r['n'], minimum=min(degrees), maximum=max(degrees),
                                spectrum_sum=sum(d * count for d, count in enumerate(h)), histogram=h))
    emit_json('serial/optimized/practical-initial-metrics.json', initial)


def main():
    global CHECK
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Verify tables without overwriting them')
    CHECK = parser.parse_args().check
    manifest = json.loads((ROOT / 'data/input-manifest.json').read_text())
    assert len(manifest['inputs']) == 177
    for entry in manifest['inputs']:
        assert sha(ROOT / entry['input']) == entry['input_sha256']
        assert sha(ROOT / entry['swaps']) == entry['swaps_sha256']
    serial = load_audit()
    assert len(serial) == 14460
    serial_tables(serial)
    main_rows=[r for r in serial if r['group']=='B' and r['input_group']=='random']
    assert {(r['n'],r['metric']) for r in main_rows} == {(n,m) for n in range(3,20) for m in ['min','max','spectrum']}
    audit = dict(status='pass', timing_records=13000, separate_count_records=1460, input_cases=177,
                 archived_record_sha256={'serial': sha(ROOT / 'data/measurements/serial.jsonl')},
                 result_tables=TABLES, same_workload_outputs_equal=True, complete_repeat_sets=True)
    if not CHECK:
        (ROOT / 'results/audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
