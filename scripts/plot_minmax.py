#!/usr/bin/env python3
"""Plot measured time and directly labelled speedups from the audited CSV tables."""
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
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import Bbox

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures'
BLUE, RED, GRAY = '#1565B5', '#B74D48', '#646D78'
METRICS = ['min', 'max']
TITLES = {'min': 'Minimum algebraic degree', 'max': 'Maximum algebraic degree update',
          'spectrum': 'Algebraic degree spectrum update'}
DATE = datetime(2026, 10, 7, tzinfo=timezone.utc)
SOURCES = {}

plt.rcParams.update({
    'font.family': 'DejaVu Serif', 'mathtext.fontset': 'stix', 'font.size': 9,
    'axes.titlesize': 10, 'axes.labelsize': 9.5, 'xtick.labelsize': 8,
    'ytick.labelsize': 8, 'legend.fontsize': 8, 'legend.frameon': False,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .8,
    'lines.linewidth': 1.9, 'lines.markersize': 4, 'grid.linewidth': .45,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'svg.hashsalt': 'fast-sbox-degree', 'figure.facecolor': 'white',
    'savefig.facecolor': 'white', 'axes.axisbelow': True,
})


def read(name):
    path = ROOT / 'results' / name
    SOURCES['results/' + name] = hashlib.sha256(path.read_bytes()).hexdigest()
    with path.open() as stream:
        return list(csv.DictReader(stream))


def vector(rows, field):
    return np.array([float(r[field]) for r in rows])


def band(ax, rows, x, field, color, label, marker='o', style='-', suffix=('_q1', '_q3')):
    y = vector(rows, field)
    ax.plot(x, y, color=color, marker=marker, linestyle=style, label=label,
            linewidth=2.1 if color == BLUE else 1.65, zorder=3 if color == BLUE else 2)
    if field + suffix[0] in rows[0]:
        ax.fill_between(x, vector(rows, field + suffix[0]), vector(rows, field + suffix[1]),
                        color=color, alpha=.60, linewidth=0)


def panel(ax, index, title):
    ax.set_title(f'({chr(97 + index)}) {title}', loc='left', pad=9)
    ax.grid(axis='y', alpha=.30)


def width_axis(ax, boundary=False):
    ax.set_xlim(2.7, 19.5)
    ax.set_xticks(range(3, 20))
    ax.set_xlabel(r'S-box size $n$ ($n\times n$)')
    if boundary:
        ax.axvline(8.5, color='#A6ADB5', linewidth=.75, linestyle=':')


def ratio_axis(ax, label='Speedup of our method (×)'):
    ax.axhline(1, color=GRAY, linewidth=.85, linestyle=':')
    ax.set_ylabel(label)
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6, min_n_ticks=4))


def pair():
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.7))
    fig.subplots_adjust(left=.092, right=.976, bottom=.205, top=.795, wspace=.32)
    return fig, axes


def label_point(ax, x, y, color=BLUE, offset=(0, 9), ha='center', bold=False):
    ax.annotate(f'{y:.2f}×', (x, y), xytext=offset, textcoords='offset points',
                ha=ha, va='bottom' if offset[1] >= 0 else 'top', color=color,
                fontsize=9 if bold else 7.5, fontweight='bold' if bold else 'normal',
                bbox=dict(facecolor='white', edgecolor='none', alpha=.88, pad=.7), zorder=5)


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in fig.axes:
        box = ax.get_tightbbox(renderer)
        assert box.x0 >= -1 and box.y0 >= -1 and box.x1 <= fig.bbox.x1 + 1 and box.y1 <= fig.bbox.y1 + 1, (name, box)
    for extension in ['pdf', 'svg', 'png']:
        metadata = {'Creator': 'Matplotlib; scripts/plot_minmax.py'} if extension != 'png' else {'Software': 'Matplotlib; scripts/plot_minmax.py'}
        if extension == 'pdf':
            metadata.update(CreationDate=DATE, ModDate=DATE)
        elif extension == 'svg':
            metadata['Date'] = '2026-10-07'
        crop = Bbox.from_extents(0, (28 if name=='04_search_time' else 22)/72, 7.4, 3.55 if name=='04_search_time' else 3.7) if name[:2] in ['01','02','03','04'] else None
        fig.savefig(OUT / f'{name}.{extension}', dpi=300, metadata=metadata, bbox_inches=crop)
    plt.close(fig)


def ordinary(metric):
    rows = sorted([r for r in read('minmax/scaling-results.csv') if r['metric'] == metric],
                  key=lambda r: int(r['n']))
    assert [int(r['n']) for r in rows] == list(range(3, 20))
    x = vector(rows, 'n')
    fig, (time, speed) = pair()
    band(time, rows, x, 'proposed_us', BLUE, 'Our method')
    for base, label, marker in [('peigen', 'PEIGEN (3–8 bits)', 's'),
                                 ('bitwise', 'Bitwise method (9–19 bits)', '^')]:
        subset = [r for r in rows if r['baseline'] == base]
        band(time, subset, vector(subset, 'n'), 'baseline_us', RED, label, marker, '--')
    time.set_yscale('log')
    time.set_ylabel('Execution time per ' + ('calculation' if metric == 'min' else 'update') + ' (µs)')
    time.legend(loc='upper left', fontsize=7.5)
    panel(time, 0, 'Execution time ↓')
    width_axis(time, True)
    band(speed, rows, x, 'paired_speedup', BLUE, 'Speedup')
    ratio_axis(speed)
    speed.set_ylim(0, max(vector(rows, 'paired_speedup_q3')) * 1.29)
    panel(speed, 1, 'Speedup ↑')
    width_axis(speed, True)
    for r in rows:
        n = int(r['n'])
        if n in ([3, 8, 12, 16, 19] if metric != 'spectrum' else [4, 8, 12, 14, 16, 19]):
            offset = (-3, 10) if n == 19 else (0, 8)
            if metric == 'spectrum' and n in [3, 6, 14]:
                offset = (0, -17)
            label_point(speed, n, float(r['paired_speedup']), offset=offset,
                        ha='right' if n == 19 else 'center', bold=n == 19)
    return fig


BUILDERS = {'01_minimum_degree': lambda: ordinary('min'), '02_maximum_degree_update': lambda: ordinary('max')}


def main():
    global OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--only', choices=list(BUILDERS), help='Regenerate one plot after a layout edit')
    parser.add_argument('--output', type=Path, default=OUT, help='Output directory; default: figures/')
    args = parser.parse_args()
    OUT = args.output
    for name, build in BUILDERS.items():
        if args.only is None or name == args.only:
            save(build(), name)
            print('Rendered', name)
    path = OUT / 'minmax-manifest.json'
    manifest = {'data_sha256': SOURCES,
                'plot_script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'matplotlib': matplotlib.__version__,
                'numpy': np.__version__, 'font': 'DejaVu Serif', 'width_inches': 7.4,
                'png_dpi': 300, 'vector_text': 'PDF embedded TrueType; SVG live text',
                'exports_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                   for p in sorted(OUT.iterdir()) if p.stem in BUILDERS and p.suffix in {'.png', '.svg', '.pdf'}}}
    path.write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()
