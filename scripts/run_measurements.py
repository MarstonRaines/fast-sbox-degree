#!/usr/bin/env python3
"""Run a saved experimental protocol into a separate output directory."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
BATCHES = {
    'minmax': 'data/measurements/minmax',
    'spectrum': 'data/measurements/spectrum/formal',
    'supplemental': 'data/measurements/spectrum/followups',
    'search': 'data/measurements/search',
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def job_command(batch, job):
    if batch == 'minmax':
        args = [f"build/minmax-{job['variant']}", job['method'], job['metric'],
                f"data/inputs/{job['input']}.lut", f"data/inputs/{job['input']}.swaps", str(job['steps'])]
        name = f"{job['input']}-{job['metric']}-{job['repetition']}-{job['method']}.json"
    elif batch == 'spectrum':
        args = [f"build/rank-{job['variant']}", 'bench', job['method'],
                f"data/inputs/{job['input']}.lut", f"data/inputs/{job['input']}.swaps", str(job['steps'])]
        name = f"{job['input']}--{job['repetition']}--{job['method']}.json"
    elif batch == 'supplemental':
        args = ['build/' + job['binary'], *job['args']]
        name = job['id'] + '.json'
    else:
        args = ['build/rank-search', job['mode'], job['method'], f"data/literature/inputs/{job['input']}.lut"]
        name = job['id'] + '.json'
    return args, name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('batch', choices=BATCHES)
    parser.add_argument('--cpu', type=int, default=0)
    parser.add_argument('--numa', type=int, default=0)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    protocol = json.loads((ROOT / BATCHES[args.batch] / 'protocol.json').read_text())
    jobs = protocol['jobs']
    queue = []
    for job in jobs:
        command, name = job_command(args.batch, job)
        assert all((ROOT / a).is_file() for a in command if a.endswith(('.lut', '.swaps')))
        queue.append((job, command, name))
    if args.dry_run:
        print(json.dumps({'batch': args.batch, 'jobs': len(queue), 'cpu': args.cpu, 'numa': args.numa,
                          'first_command': queue[0][1], 'last_command': queue[-1][1]}, indent=2))
        return
    for _, command, _ in queue:
        if not (ROOT / command[0]).is_file():
            parser.error('Build the programs with: make -j3 all diagnostics')
    output = ROOT / 'runs' / args.batch
    output.mkdir(parents=True, exist_ok=False)
    protocol['binding'] = {'cpu': args.cpu, 'numa': args.numa, 'threads': 1}
    def save(name, obj):
        (output / name).write_text(json.dumps(obj, indent=2) + '\n')
    save('protocol.json', protocol)
    sources = [ROOT / 'Makefile', Path(__file__), *(ROOT / 'src').rglob('*.cpp'), *(ROOT / 'src').rglob('*.hpp')]
    save('source-hashes.json', {str(p.relative_to(ROOT)): sha(p) for p in sources})
    env = {**os.environ, 'OMP_NUM_THREADS': '1', 'LC_ALL': 'C'}
    for index, (job, program, name) in enumerate(queue, 1):
        command = ['numactl', f'--physcpubind={args.cpu}', f'--membind={args.numa}', *program]
        started = datetime.datetime.now(datetime.timezone.utc).isoformat()
        start = time.monotonic()
        run = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True)
        record = {**job, 'command': command, 'started_utc': started,
                  'process_seconds': time.monotonic() - start, 'returncode': run.returncode,
                  'binary_sha256': sha(ROOT / program[0])}
        inputs = {p: sha(ROOT / p) for p in program if p.endswith(('.lut', '.swaps'))}
        if args.batch == 'supplemental':
            record['input_sha256'] = inputs
        else:
            record['input_sha256'] = next(v for p, v in inputs.items() if p.endswith('.lut'))
            if args.batch != 'search':
                record['swaps_sha256'] = next(v for p, v in inputs.items() if p.endswith('.swaps'))
        if run.returncode:
            record['stderr'] = run.stderr
            save(name, record)
            raise RuntimeError(f'{name} failed; see {output}')
        record['result'] = json.loads(run.stdout)
        save(name, record)
        if index % 100 == 0 or index == len(queue):
            print(f'{index}/{len(queue)}', flush=True)
    complete = {'passed': True, 'records': len(queue)}
    if args.batch == 'search':
        complete['binary_sha256'] = sha(ROOT / 'build/rank-search')
    save('COMPLETE.json', complete)


if __name__ == '__main__':
    main()
