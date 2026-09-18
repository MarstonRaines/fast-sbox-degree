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

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures'
BLUE, RED, GRAY = '#1565B5', '#B74D48', '#646D78'
METRICS = ['min', 'max', 'spectrum']
TITLES = {'min': 'Minimum algebraic degree', 'max': 'Maximum algebraic degree update',
          'spectrum': 'Complete algebraic degree spectrum update'}
DATE = datetime(2026, 9, 19, tzinfo=timezone.utc)
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
                        color=color, alpha=.10, linewidth=0)


def panel(ax, index, title):
    ax.set_title(f'({chr(97 + index)}) {title}', loc='left', pad=9)
    ax.grid(axis='y', alpha=.30)


def width_axis(ax, boundary=False):
    ax.set_xlim(2.7, 16.5)
    ax.set_xticks(range(3, 17))
    ax.set_xlabel(r'S-box size $n$ ($n\times n$)')
    if boundary:
        ax.axvline(8.5, color='#A6ADB5', linewidth=.75, linestyle=':')


def ratio_axis(ax, label='Speedup of our method (×)'):
    ax.axhline(1, color=GRAY, linewidth=.85, linestyle=':')
    ax.set_ylabel(label)
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6, min_n_ticks=4))


def pair(title, note):
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.7))
    fig.subplots_adjust(left=.092, right=.976, bottom=.205, top=.795, wspace=.32)
    fig.suptitle(title, y=.97, fontsize=12, fontweight='bold')
    fig.text(.5, .045, note, ha='center', fontsize=7.5, color='#404750')
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
        metadata = {'Creator': 'Matplotlib; scripts/plot_results.py'} if extension != 'png' else {'Software': 'Matplotlib; scripts/plot_results.py'}
        if extension == 'pdf':
            metadata.update(CreationDate=DATE, ModDate=DATE)
        elif extension == 'svg':
            metadata['Date'] = '2026-09-19'
        fig.savefig(OUT / f'{name}.{extension}', dpi=300, metadata=metadata)
    plt.close(fig)


def ordinary(metric):
    rows = sorted([r for r in read('serial/optimized/scaling-results.csv') if r['metric'] == metric],
                  key=lambda r: int(r['n']))
    assert [int(r['n']) for r in rows] == list(range(3, 17))
    x = vector(rows, 'n')
    fig, (time, speed) = pair(TITLES[metric],
        '10 random inputs per size; 5 repeats per input. Bands: middle 50% across inputs.')
    band(time, rows, x, 'proposed_us', BLUE, 'Our method')
    for base, label, marker in [('peigen', 'PEIGEN (3–8 bits)', 's'),
                                 ('bitwise', 'Bitwise baseline (9–16 bits)', '^')]:
        subset = [r for r in rows if r['baseline'] == base]
        band(time, subset, vector(subset, 'n'), 'baseline_us', RED, label, marker, '--')
    time.set_yscale('log')
    time.set_ylabel('Core time per ' + ('evaluation' if metric == 'min' else 'update') + ' (µs)')
    time.legend(loc='upper left', fontsize=7.5)
    panel(time, 0, 'Computation time ↓')
    width_axis(time, True)
    band(speed, rows, x, 'paired_speedup', BLUE, 'Core computation')
    band(speed, rows, x, 'total_speedup', GRAY, 'Including initialization', 's', '--')
    ratio_axis(speed)
    speed.set_ylim(0, max(vector(rows, 'paired_speedup_q3')) * 1.29)
    speed.legend(loc='upper left', fontsize=7.5)
    panel(speed, 1, 'Speedup over the baseline ↑')
    width_axis(speed, True)
    for r in rows:
        n = int(r['n'])
        if n in [3, 6, 8, 12, 16] or (metric == 'spectrum' and n in [4, 14]):
            offset = (-3, 10) if n == 16 else (0, 8)
            if metric == 'spectrum' and n in [3, 6, 14]:
                offset = (0, -17)
            label_point(speed, n, float(r['paired_speedup']), offset=offset,
                        ha='right' if n == 16 else 'center', bold=n == 16)
    label_point(speed, 16, float(rows[-1]['total_speedup']), GRAY, (-3, -15), 'right')
    return fig


def search():
    data = read('serial/optimized/search-results.csv')
    fig, ax = plt.subplots(figsize=(7.4, 3.75))
    fig.subplots_adjust(left=.105, right=.975, bottom=.19, top=.77)
    fig.suptitle('Faster completion of the same 8-bit search', y=.965, fontsize=12, fontweight='bold')
    for j, metric in enumerate(METRICS):
        rows = sorted([r for r in data if r['metric'] == metric and r['group'] == 'random'], key=lambda r: r['case'])
        assert len(rows) == 10
        for field, shift, color, label in [('baseline_total_ms', -.19, RED, 'PEIGEN'),
                                           ('proposed_total_ms', .19, BLUE, 'Our method')]:
            values = vector(rows, field)
            y = median(values)
            q1, q3 = np.quantile(values, [.25, .75])
            ax.bar(j + shift, y, width=.32, color=color, alpha=.90, label=label if j == 0 else None)
            ax.errorbar(j + shift, y, yerr=[[y - q1], [q3 - y]], color='#313840', fmt='none', capsize=3, linewidth=.9)
            ax.scatter(j + shift + np.linspace(-.075, .075, 10), values, s=9,
                       color='white', edgecolor=color, linewidth=.6, zorder=3)
            ax.text(j + shift, max(values) + 3.2, f'{y:.2f} ms', ha='center', fontsize=8.5,
                    color=color, fontweight='bold' if color == BLUE else 'normal')
        ratio = median(vector(rows, 'total_speedup'))
        ax.text(j, 137, f'{ratio:.2f}× faster', ha='center', color=BLUE, fontsize=12, fontweight='bold')
    ax.set_ylim(0, 153)
    ax.set_xticks(np.arange(3), ['Minimum degree', 'Maximum degree', 'Spectrum sum'])
    ax.set_ylabel('Total search time (ms)')
    ax.grid(axis='y', alpha=.3)
    ax.legend(loc='upper center', bbox_to_anchor=(.5, 1.23), ncol=2)
    fig.text(.5, .078, 'Time includes initialization. Same candidates, accepted swaps and final result.',
             ha='center', fontsize=7.5, color='#404750')
    fig.text(.5, .025, '10 starts × 5 repeats; bars: median; whiskers: middle 50%; dots: per-start medians.',
             ha='center', fontsize=7.5, color='#404750')
    return fig


def practical():
    data = read('serial/optimized/practical-results.csv')
    cases = [('PRESENT', 'PRESENT', '4', 'PEIGEN'), ('AES', 'AES', '8', 'PEIGEN'),
             ('CLEFIA_S0', 'CLEFIA S0', '8', 'PEIGEN'), ('CLEFIA_S1', 'CLEFIA S1', '8', 'PEIGEN'),
             ('MISTY1_S9', 'MISTY1 S9', '9', 'Bitwise'), ('MISTY1_FI_key0000', 'MISTY1 FI', '16', 'Bitwise')]
    cells = []
    for case, label, n, baseline in cases:
        row = [label, n, baseline]
        for metric in METRICS:
            r = next(r for r in data if r['case'] == case and r['metric'] == metric)
            row.append(f"{float(r['paired_speedup']):.2f}× ({float(r['total_speedup']):.2f}×)")
        cells.append(row)
    fig, ax = plt.subplots(figsize=(7.4, 2.8))
    fig.subplots_adjust(left=.025, right=.975, bottom=.20, top=.79)
    ax.axis('off')
    fig.suptitle('Practical instances: speedup over the optimized baseline', y=.965, fontsize=11, fontweight='bold')
    fig.text(.5, .825, 'Core speedup (speedup including initialization)', ha='center', fontsize=8.5, color=GRAY)
    table = ax.table(cellText=cells, colLabels=['Instance', 'Bits', 'Baseline', 'Minimum degree',
                     'Maximum update', 'Spectrum update'], colWidths=[.17, .055, .105, .24, .215, .215],
                     cellLoc='center', loc='center', bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    for (i, j), cell in table.get_celld().items():
        cell.set_edgecolor('#D5DCE4')
        cell.set_linewidth(.45)
        if i == 0:
            cell.set_facecolor('#EAF1F8')
            cell.get_text().set_fontweight('bold')
        else:
            cell.set_facecolor('#F6F8FB' if i % 2 == 0 else 'white')
            if j >= 3:
                cell.get_text().set_color(BLUE)
                cell.get_text().set_fontweight('bold')
    fig.text(.5, .075, 'Five repeats per instance. Ratios below 1 favor the baseline. FI uses fixed subkey 0.',
             ha='center', fontsize=7.5, color='#404750')
    fig.text(.5, .025, 'Update measurements follow fixed transpositions starting from each named instance.',
             ha='center', fontsize=7.5, color='#404750')
    return fig


BUILDERS = {
    '01_minimum_degree': lambda: ordinary('min'),
    '02_maximum_degree_update': lambda: ordinary('max'),
    '03_degree_spectrum_update': lambda: ordinary('spectrum'),
    '04_search_time': search,
    'table01_practical_instances': practical,
}


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
    path = OUT / 'manifest.json'
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
