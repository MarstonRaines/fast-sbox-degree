#!/usr/bin/env python3
"""Replay saved workloads with the release build; archived measurements stay immutable."""
import argparse
import hashlib
import json
import os
import platform
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite', choices=['serial', 'parallel'], required=True)
    parser.add_argument('--group', choices=['A', 'B', 'search', 'counts'])
    parser.add_argument('--case', help='Exact input ID, for example random_n08_00')
    parser.add_argument('--metric', choices=['min', 'max', 'spectrum'])
    parser.add_argument('--method', choices=['reference', 'reference-proposed', 'peigen', 'bitwise', 'proposed'])
    parser.add_argument('--threads', type=int)
    parser.add_argument('--limit', type=int, help='Run the first N matching tasks (smoke checks)')
    parser.add_argument('--cpu-start', type=int, help='First of contiguous physical CPUs; verify topology first')
    parser.add_argument('--numa-node', type=int, help='Bind memory using numactl; requires --cpu-start')
    parser.add_argument('--output', type=Path, help='New, empty output directory')
    parser.add_argument('--list', action='store_true', help='Print matching tasks without running')
    args = parser.parse_args()
    if args.numa_node is not None and args.cpu_start is None:
        parser.error('--numa-node requires --cpu-start')
    manifest = ROOT / f'protocol/{args.suite}-tasks.json'
    tasks = json.loads(manifest.read_text())
    for key in ['group', 'case', 'metric', 'method', 'threads']:
        value = getattr(args, key)
        if value is not None:
            tasks = [t for t in tasks if t[key] == value]
    if args.limit is not None:
        if args.limit < 1:
            parser.error('--limit must be positive')
        tasks = tasks[:args.limit]
    if not tasks:
        parser.error('No matching tasks')
    if args.list:
        print(json.dumps(tasks, indent=2))
        return
    output = (args.output or ROOT / 'runs' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')).resolve()
    if output.exists() and any(output.iterdir()):
        parser.error('Output directory must be empty; existing records are never overwritten')
    (output / 'traces').mkdir(parents=True, exist_ok=True)
    binaries = {}
    for task in tasks:
        name = ('counts' if task['group'] == 'counts' else
                'serial' if args.suite == 'serial' or task['binary'] == 'coam-serial-source' else 'parallel')
        binary = ROOT / 'build' / name
        if not binary.is_file():
            parser.error(f'Missing {binary.name}; run make all')
        binaries[name] = sha(binary)
    metadata = {'created_utc': datetime.now(timezone.utc).isoformat(),
                'platform': platform.platform(), 'suite': args.suite,
                'task_manifest_sha256': sha(manifest), 'binary_sha256': binaries,
                'cpu_start': args.cpu_start, 'numa_node': args.numa_node,
                'tasks': len(tasks), 'purpose': 'Rerun using publication sources; independent of archived timings.'}
    (output / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    (output / 'tasks.json').write_text(json.dumps(tasks, indent=2) + '\n')
    archived = {(r['group'], r['id']): r for line in
                (ROOT / f'data/measurements/{args.suite}.jsonl').read_text().splitlines()
                if (r := json.loads(line))}
    for index, task in enumerate(tasks, 1):
        threads = task['threads']
        name = ('counts' if task['group'] == 'counts' else
                'serial' if args.suite == 'serial' or task['binary'] == 'coam-serial-source' else 'parallel')
        env = {**os.environ, 'OMP_NUM_THREADS': str(threads), 'OMP_DYNAMIC': 'false',
               'OMP_PROC_BIND': 'close', 'OMP_WAIT_POLICY': 'PASSIVE',
               'COAM_OMP_MIN_WORK': '32768', 'OMP_PLACES': 'cores'}
        command = []
        if args.cpu_start is not None:
            cpus = list(range(args.cpu_start, args.cpu_start + threads))
            env['OMP_PLACES'] = ','.join('{%d}' % cpu for cpu in cpus)
            command = ['numactl', '--physcpubind=' + ','.join(map(str, cpus))]
            if args.numa_node is not None:
                command.append('--membind=' + str(args.numa_node))
        trace = output / 'traces' / (task['id'] + '.json')
        command += [str(ROOT / 'build' / name), 'bench', task['method'], task['metric'],
                    task['mode'], str(ROOT / task['input']),
                    '-' if task['swaps'] == '-' else str(ROOT / task['swaps']), str(task['loops']), str(trace)]
        start = time.monotonic()
        completed = subprocess.run(command, env=env, capture_output=True, text=True)
        if completed.returncode:
            (output / 'failure.txt').write_text(completed.stdout + completed.stderr)
            raise SystemExit(f'Task failed: {task["id"]}; see failure.txt')
        result = json.loads(completed.stdout)
        expected = archived[task['group'], task['id']]
        for field in ['steps', 'loops', 'initial_value', 'final_value', 'candidate_digest',
                      'final_box_hash', 'final_histogram', 'accepted']:
            if result[field] != expected[field]:
                raise RuntimeError(f'Saved {field} mismatch: {task["id"]}')
        full = json.loads(trace.read_text())
        content = json.dumps({k: full[k] for k in ['final_box', 'final_histogram', 'accepted_swaps']},
                             separators=(',', ':'), sort_keys=True) + '\n'
        if hashlib.sha256(content.encode()).hexdigest() != task['expected_output_sha256']:
            raise RuntimeError('Saved output mismatch: ' + task['id'])
        result.update(task)
        result.update(process_seconds=time.monotonic() - start, release_binary=name)
        with (output / 'records.jsonl').open('a') as file:
            file.write(json.dumps(result) + '\n')
        print(f'{index}/{len(tasks)} {task["id"]}: output verified', flush=True)
    (output / 'complete.json').write_text(json.dumps({'status': 'pass', 'tasks': len(tasks)}) + '\n')


if __name__ == '__main__':
    main()
