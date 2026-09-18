#!/usr/bin/env python3
"""Render publication figures directly from the audited result tables."""
import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from statistics import median

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter, NullFormatter, MaxNLocator

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures'
BLUE, RED, TEAL, GOLD, GRAY = '#0F4D92', '#B64342', '#42949E', '#A36A16', '#767676'
METRICS = ['min', 'max', 'spectrum']
TITLES = ['Minimum degree', 'Maximum degree update', 'Spectrum update']
DATE = datetime(2026, 9, 18, tzinfo=timezone.utc)
SOURCES = {}

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.titlesize': 9,
    'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'legend.fontsize': 7.5, 'legend.frameon': False,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .7,
    'lines.linewidth': 1.4, 'lines.markersize': 3.4, 'grid.linewidth': .4,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'svg.hashsalt': 'fast-sbox-degree-v1', 'figure.facecolor': 'white',
    'savefig.facecolor': 'white', 'axes.axisbelow': True,
})


def read(name):
    path = ROOT / 'results' / name
    SOURCES['results/' + name] = hashlib.sha256(path.read_bytes()).hexdigest()
    with path.open() as stream:
        return list(csv.DictReader(stream))


def vector(rows, field):
    return np.array([float(r[field]) for r in rows])


def band(ax, rows, x, field, color, label=None, marker='o', style='-', suffix=('_q1', '_q3')):
    y = vector(rows, field)
    ax.plot(x, y, color=color, marker=marker, linestyle=style, label=label)
    if field + suffix[0] in rows[0]:
        ax.fill_between(x, vector(rows, field + suffix[0]), vector(rows, field + suffix[1]),
                        color=color, alpha=.12, linewidth=0)


def panel(ax, index, title):
    ax.set_title(f'({chr(97 + index)}) {title}', loc='left', pad=8)
    ax.grid(axis='y', alpha=.25)


def width_axis(ax, boundary=False):
    ax.set_xlim(2.7, 16.3)
    ax.set_xticks([3, 6, 8, 10, 13, 16])
    ax.set_xlabel('S-box width n (bits)')
    if boundary:
        ax.axvline(8.5, color='#BEBEBE', linewidth=.7, linestyle=':')


def thread_axis(ax):
    ax.set_xscale('log', base=2)
    ax.set_xlim(.9, 71)
    ax.set_xticks([1, 2, 4, 8, 16, 32, 64])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_xlabel('Physical-thread cap')


def save(fig, name):
    OUT.mkdir(exist_ok=True)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in fig.axes:
        box = ax.get_tightbbox(renderer)
        assert box.x0 >= -1 and box.y0 >= -1 and box.x1 <= fig.bbox.x1 + 1 and box.y1 <= fig.bbox.y1 + 1, (name, box)
    for extension in ['pdf', 'svg', 'png']:
        metadata = {'Creator': 'Matplotlib; scripts/plot_results.py'} if extension != 'png' else {'Software': 'Matplotlib; scripts/plot_results.py'}
        if extension == 'pdf':
            metadata.update(CreationDate=DATE, ModDate=DATE)
        elif extension == 'svg':
            metadata['Date'] = '2026-09-18'
        fig.savefig(OUT / f'{name}.{extension}', dpi=300, metadata=metadata)
    plt.close(fig)


def optimized_scaling():
    data = read('serial/optimized/scaling-results.csv')
    fig, axes = plt.subplots(2, 3, figsize=(7.4, 4.8))
    fig.subplots_adjust(left=.085, right=.985, bottom=.11, top=.88, hspace=.48, wspace=.36)
    for j, metric in enumerate(METRICS):
        rows = sorted([r for r in data if r['metric'] == metric], key=lambda r: int(r['n']))
        x = vector(rows, 'n')
        ax = axes[0, j]
        band(ax, rows, x, 'proposed_us', BLUE, 'Proposed', 'o')
        for baseline, color, marker, label in [('peigen', GOLD, 's', 'PEIGEN'), ('bitwise', RED, '^', 'Packed traditional')]:
            subset = [r for r in rows if r['baseline'] == baseline]
            band(ax, subset, vector(subset, 'n'), 'baseline_us', color, label, marker)
        ax.set_yscale('log')
        panel(ax, j, TITLES[j])
        width_axis(ax, True)
        ax = axes[1, j]
        band(ax, rows, x, 'paired_speedup', BLUE, 'Core', 'o')
        band(ax, rows, x, 'total_speedup', TEAL, 'Batch total', 's', '--')
        ax.axhline(1, color=GRAY, linewidth=.8, linestyle=':')
        ax.set_yscale('log')
        if metric == 'spectrum':
            ax.set_yticks([.5, 1, 2, 4])
            ax.yaxis.set_major_formatter(ScalarFormatter())
            ax.yaxis.set_minor_formatter(NullFormatter())
        panel(ax, j + 3, 'Baseline / proposed')
        width_axis(ax, True)
    axes[0, 0].set_ylabel('Core time per step (µs)')
    axes[1, 0].set_ylabel('Time ratio (×)')
    handles, labels = axes[0, 0].get_legend_handles_labels()
    h, l = axes[1, 0].get_legend_handles_labels()
    fig.legend(handles + h, labels + l, loc='upper center', bbox_to_anchor=(.52, .985), ncol=5, columnspacing=1.4)
    save(fig, '01_optimized_scaling')


def reference_scaling():
    data = read('serial/reference-scaling.csv')
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.8))
    fig.subplots_adjust(left=.085, right=.985, bottom=.19, top=.78, wspace=.38)
    for j, metric in enumerate(METRICS):
        ax = axes[j]
        for profile, color, label in [('static' if metric == 'min' else 'natural', BLUE, 'Static / natural swaps'),
                                      ('high', TEAL, 'High-update swaps')]:
            rows = sorted([r for r in data if r['metric'] == metric and r['profile'] == profile], key=lambda r: int(r['n']))
            if not rows:
                continue
            band(ax, rows, vector(rows, 'n'), 'paired_speedup', color, label,
                 marker='o' if color == BLUE else 's')
        x = np.arange(3, 17)
        factor = 2. ** x / x if metric == 'min' else x if metric == 'max' else 2 * x
        ax.plot(x, factor, color=GRAY, linestyle='--', linewidth=1.1, label='Reference growth factor')
        ax.set_yscale('log')
        panel(ax, j, TITLES[j])
        width_axis(ax)
        label = r'$2^n/n$' if metric == 'min' else r'$n$' if metric == 'max' else r'$2n$'
        ax.text(.04, .94, 'Growth factor: ' + label, transform=ax.transAxes, fontsize=7, va='top')
    axes[0].set_ylabel('Scalar reference / proposed (×)')
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.52, .995), ncol=3)
    save(fig, '02_reference_scaling')


def practical():
    data = read('serial/optimized/practical-results.csv')
    cases = ['PRESENT', 'AES', 'CLEFIA_S0', 'CLEFIA_S1', 'MISTY1_S9', 'MISTY1_FI_key0000']
    labels = ['PRESENT (4)', 'AES (8)', 'CLEFIA S0 (8)', 'CLEFIA S1 (8)', 'MISTY1 S9 (9)', 'MISTY1 FI (16)']
    fig, axes = plt.subplots(1, 3, sharey=True, figsize=(7.4, 3.0))
    fig.subplots_adjust(left=.17, right=.985, bottom=.20, top=.80, wspace=.32)
    for j, metric in enumerate(METRICS):
        rows = [next(r for r in data if r['case'] == case and r['metric'] == metric) for case in cases]
        ax = axes[j]
        for field, color, marker, shift, label in [('paired_speedup', BLUE, 'o', -.11, 'Core'),
                                                  ('total_speedup', TEAL, 's', .11, 'Batch total')]:
            ax.scatter(vector(rows, field), np.arange(6) + shift, color=color, marker=marker, s=23, label=label)
        ax.axvline(1, color=GRAY, linewidth=.8, linestyle=':')
        ax.set_xscale('log')
        if metric == 'spectrum':
            ax.set_xticks([1, 2, 4])
            ax.xaxis.set_major_formatter(ScalarFormatter())
            ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_ylim(5.5, -.5)
        ax.set_yticks(np.arange(6), labels)
        ax.grid(axis='x', alpha=.25)
        panel(ax, j, TITLES[j])
        ax.set_xlabel('Baseline / proposed (×)')
    fig.legend(*axes[0].get_legend_handles_labels(), loc='upper center', bbox_to_anchor=(.56, .98), ncol=2)
    save(fig, '03_practical_instances')


def search():
    data = read('serial/optimized/search-results.csv')
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.9))
    fig.subplots_adjust(left=.085, right=.985, bottom=.18, top=.84, wspace=.34)
    for j, metric in enumerate(METRICS):
        rows = [r for r in data if r['metric'] == metric and r['group'] == 'random']
        assert len(rows) == 10
        ax = axes[j]
        base, own = vector(rows, 'baseline_total_ms'), vector(rows, 'proposed_total_ms')
        for a, b in zip(base, own):
            ax.plot([0, 1], [a, b], color='#ABB4C0', alpha=.65, linewidth=.8)
        ax.scatter(np.zeros(10), base, color=RED, s=15, zorder=3)
        ax.scatter(np.ones(10), own, color=BLUE, s=15, zorder=3)
        ax.set_ylim(0, max(base) * 1.10)
        ax.set_xlim(-.23, 1.23)
        ax.set_xticks([0, 1], ['PEIGEN', 'Proposed'])
        panel(ax, j, ['Minimum degree', 'Maximum degree', 'Spectrum sum'][j])
        ratio = median(vector(rows, 'total_speedup'))
        ax.text(.5, 1.015, f'{ratio:.2f}× median paired speedup', ha='center', transform=ax.transAxes, fontsize=7)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
    axes[0].set_ylabel('Complete search total time (ms)')
    save(fig, '04_complete_search')


def parallel_scaling():
    data = read('parallel/hardware-scaling.csv')
    fig, axes = plt.subplots(2, 3, figsize=(7.4, 4.8))
    fig.subplots_adjust(left=.085, right=.985, bottom=.11, top=.88, hspace=.50, wspace=.36)
    for j, metric in enumerate(METRICS):
        for i, field in enumerate(['core', 'total']):
            ax = axes[i, j]
            for sample, style in [('00', '-'), ('01', '--')]:
                for method, color, marker in [('proposed', BLUE, 'o'), ('bitwise', RED, '^')]:
                    rows = sorted([r for r in data if r['group'] == 'B' and r['binary'] == 'coam'
                                   and r['case'] == 'random_n16_' + sample and r['method'] == method
                                   and r['metric'] == metric], key=lambda r: int(r['threads']))
                    assert len(rows) == 7
                    x = vector(rows, 'threads')
                    ax.plot(x, vector(rows, field + '_speedup'), color=color, linestyle=style, marker=marker)
                    ax.fill_between(x, vector(rows, field + '_q1'), vector(rows, field + '_q3'), color=color, alpha=.07, linewidth=0)
            ax.axhline(1, color=GRAY, linewidth=.7, linestyle=':')
            ax.set_ylim(bottom=0)
            thread_axis(ax)
            panel(ax, i * 3 + j, TITLES[j])
    axes[0, 0].set_ylabel('Core hardware speedup (×)')
    axes[1, 0].set_ylabel('Batch-total hardware speedup (×)')
    handles = [Line2D([], [], color=BLUE, marker='o', label='Proposed'),
               Line2D([], [], color=RED, marker='^', label='Packed traditional'),
               Line2D([], [], color=GRAY, label='Input 00'),
               Line2D([], [], color=GRAY, linestyle='--', label='Input 01')]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.52, .985), ncol=4)
    save(fig, '05_parallel_scaling')


def parallel_comparison():
    data = read('parallel/algorithm-comparison.csv')
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.9))
    fig.subplots_adjust(left=.085, right=.985, bottom=.19, top=.79, wspace=.34)
    for j, metric in enumerate(METRICS):
        ax = axes[j]
        for sample, style in [('00', '-'), ('01', '--')]:
            rows = sorted([r for r in data if r['group'] == 'B' and r['case'] == 'random_n16_' + sample
                           and r['baseline'] == 'bitwise' and r['metric'] == metric], key=lambda r: int(r['threads']))
            assert len(rows) == 7
            for field, color, marker in [('paired_algorithm_speedup', BLUE, 'o'), ('paired_batch_total_speedup', TEAL, 's')]:
                ax.plot(vector(rows, 'threads'), vector(rows, field), color=color, marker=marker, linestyle=style)
        ax.axhline(1, color=GRAY, linewidth=.8, linestyle=':')
        ax.set_yscale('log')
        panel(ax, j, TITLES[j])
        thread_axis(ax)
    axes[0].set_ylabel('Packed traditional / proposed (×)')
    handles = [Line2D([], [], color=BLUE, marker='o', label='Core'),
               Line2D([], [], color=TEAL, marker='s', label='Batch total'),
               Line2D([], [], color=GRAY, label='Input 00'),
               Line2D([], [], color=GRAY, linestyle='--', label='Input 01')]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.52, .995), ncol=4)
    save(fig, '06_equal_thread_comparison')


def main():
    functions = {f.__name__: f for f in [optimized_scaling, reference_scaling, practical, search, parallel_scaling, parallel_comparison]}
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--only', choices=list(functions), help='Regenerate one figure after a local layout change')
    args = parser.parse_args()
    for name, function in functions.items():
        if args.only is None or name == args.only:
            function()
            print('Rendered', name)
    path = OUT / 'manifest.json'
    old_sources = json.loads(path.read_text()).get('data_sha256', {}) if path.exists() else {}
    manifest = {'data_sha256': {**old_sources, **SOURCES}, 'matplotlib': matplotlib.__version__,
                'numpy': np.__version__, 'font': 'DejaVu Sans', 'width_inches': 7.4,
                'png_dpi': 300, 'vector_text': 'PDF embedded TrueType; SVG live text',
                'exports_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                   for p in sorted(OUT.iterdir()) if p.suffix in {'.png', '.svg', '.pdf'}}}
    path.write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()
