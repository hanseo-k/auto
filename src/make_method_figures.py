"""분석 방법 PPT용 설명 그림 생성 (영문 라벨로 폰트 문제 회피)."""
import os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT_DIR = '/Users/gimhanseo/Desktop/공프/자동분석폴더/figures'
os.makedirs(OUT_DIR, exist_ok=True)


def fig_hampel_mechanism():
    fig, ax = plt.subplots(figsize=(10, 6.5), dpi=160)
    np.random.seed(0)
    base_vals = 33 + np.random.normal(0, 1.2, (5, 5))
    base_vals[2, 2] = 8

    for i in range(5):
        for j in range(5):
            x, y = j, 4 - i
            v = base_vals[i, j]
            color, ec, lw = 'lightblue', 'navy', 0.8
            if (i, j) == (2, 2):
                color, ec, lw = 'lightcoral', 'red', 2.5
            elif abs(i - 2) <= 1 and abs(j - 2) <= 1:
                color, ec, lw = '#FFE4B5', '#D2691E', 1.8
            ax.add_patch(Rectangle((x - 0.45, y - 0.45), 0.9, 0.9,
                                    facecolor=color, edgecolor=ec, lw=lw))
            ax.text(x, y, f'{v:.1f}', ha='center', va='center',
                    fontsize=11, fontweight='bold')

    ax.annotate('', xy=(3.4, 2), xytext=(2.5, 2),
                arrowprops=dict(arrowstyle='->', color='red', lw=2))
    ax.annotate('', xy=(1.6, 2), xytext=(2.5, 2),
                arrowprops=dict(arrowstyle='->', color='red', lw=2))
    ax.annotate('', xy=(2, 1.6), xytext=(2, 2.5),
                arrowprops=dict(arrowstyle='->', color='red', lw=2))

    ax.text(5.5, 3.7,
            'Suspect die (red): 8\n'
            'Median of 8 neighbors: 33.2\n'
            'MAD of neighbors: ~1.0\n\n'
            '|8 - 33.2| = 25.2\n'
            'Threshold (3 x MAD): 3.0\n\n'
            '25.2 > 3.0  -->  outlier',
            fontsize=11, va='top', ha='left', family='monospace',
            bbox=dict(boxstyle='round', facecolor='#fff8dc',
                      edgecolor='#daa520', lw=1.5))

    ax.set_xlim(-0.7, 9); ax.set_ylim(-1.0, 5)
    ax.set_aspect('equal'); ax.axis('off')
    ax.set_title('Hampel Filter Mechanism: 8-Neighbor Comparison',
                 fontsize=13, fontweight='bold')

    ax.add_patch(Rectangle((0, -0.85), 0.3, 0.3,
                            facecolor='lightcoral', edgecolor='red', lw=2))
    ax.text(0.4, -0.7, 'suspect die', fontsize=10, va='center')
    ax.add_patch(Rectangle((2.2, -0.85), 0.3, 0.3,
                            facecolor='#FFE4B5', edgecolor='#D2691E', lw=1.8))
    ax.text(2.6, -0.7, '8 neighbors', fontsize=10, va='center')
    ax.add_patch(Rectangle((4.6, -0.85), 0.3, 0.3,
                            facecolor='lightblue', edgecolor='navy', lw=0.8))
    ax.text(5.0, -0.7, 'other dies (not used)', fontsize=10, va='center')

    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, 'fig_hampel_mechanism.png'),
                bbox_inches='tight')
    plt.close()


def fig_robust_vs_3sigma():
    np.random.seed(1)
    data_clean = np.random.normal(50, 5, 13)
    data = np.append(data_clean, [120])

    mean = np.mean(data); std = np.std(data)
    median = np.median(data); mad = np.median(np.abs(data - median))

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), dpi=160)

    ax = axes[0]
    ax.scatter(range(len(data)), data, s=80, c='steelblue',
               edgecolor='navy', zorder=3)
    ax.scatter([13], [120], s=140, c='red', edgecolor='darkred', zorder=4)
    ax.axhline(mean, color='blue', ls='-', label=f'mean = {mean:.1f}')
    ax.axhspan(mean - 3*std, mean + 3*std, color='blue', alpha=0.15,
               label=f'3-sigma = +/-{3*std:.1f}')
    ax.set_title('3-sigma method\n(outlier inflates mean & sigma)',
                 fontsize=12, fontweight='bold', color='darkred')
    ax.set_ylabel('value'); ax.set_xlabel('data index')
    ax.legend(loc='lower right', fontsize=9); ax.grid(alpha=0.3)
    ax.set_ylim(20, 145)
    inside = (data >= mean - 3*std) & (data <= mean + 3*std)
    ax.text(0.5, 130, f'Inside bound: {inside.sum()}/{len(data)}\n'
                       f'Outlier missed (X)',
            fontsize=10, color='darkred', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='#ffe4e1',
                      edgecolor='red'))

    ax = axes[1]
    ax.scatter(range(len(data)), data, s=80, c='steelblue',
               edgecolor='navy', zorder=3)
    ax.scatter([13], [120], s=140, c='red', edgecolor='darkred', zorder=4)
    ax.axhline(median, color='green', ls='-', label=f'median = {median:.1f}')
    ax.axhspan(median - 3*mad, median + 3*mad, color='green', alpha=0.15,
               label=f'3 x MAD = +/-{3*mad:.1f}')
    ax.set_title('MAD method (Hampel)\n(median & MAD are robust)',
                 fontsize=12, fontweight='bold', color='darkgreen')
    ax.set_ylabel('value'); ax.set_xlabel('data index')
    ax.legend(loc='lower right', fontsize=9); ax.grid(alpha=0.3)
    ax.set_ylim(20, 145)
    inside = (data >= median - 3*mad) & (data <= median + 3*mad)
    ax.text(0.5, 130, f'Inside bound: {inside.sum()}/{len(data)}\n'
                       f'Outlier detected (V)',
            fontsize=10, color='darkgreen', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='#e8f5e8',
                      edgecolor='green'))

    plt.suptitle('Same data (13 normal + 1 outlier): detection comparison',
                 fontsize=12, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, 'fig_robust_vs_3sigma.png'),
                bbox_inches='tight')
    plt.close()


def fig_physical_bounds():
    metrics = [
        ('ER (dB)',   0,   45,  -5,  60),
        ('IL (dB)',  -20,   0, -30,   5),
        ('V_pi (V)',  5,   80,   0, 100),
    ]
    fig, axes = plt.subplots(3, 1, figsize=(11, 5.5), dpi=160)
    for i, (ax, (name, lo, hi, vmin, vmax)) in enumerate(zip(axes, metrics)):
        ax.axhspan(0.3, 0.7, xmin=(lo - vmin)/(vmax - vmin),
                   xmax=(hi - vmin)/(vmax - vmin),
                   color='lightgreen', alpha=0.7,
                   label='valid (normal)' if i == 0 else None)
        ax.axhspan(0.3, 0.7, xmin=0,
                   xmax=(lo - vmin)/(vmax - vmin),
                   color='lightcoral', alpha=0.5,
                   label='invalid (outlier)' if i == 0 else None)
        ax.axhspan(0.3, 0.7, xmin=(hi - vmin)/(vmax - vmin),
                   xmax=1, color='lightcoral', alpha=0.5)
        ax.axvline(lo, color='green', lw=2.5)
        ax.axvline(hi, color='green', lw=2.5)
        ax.text(lo, 0.85, f'{lo}', ha='center', fontsize=11,
                color='darkgreen', fontweight='bold')
        ax.text(hi, 0.85, f'{hi}', ha='center', fontsize=11,
                color='darkgreen', fontweight='bold')
        ax.text((lo + hi) / 2, 0.5, f'valid: [{lo}, {hi}]',
                ha='center', va='center', fontsize=12, fontweight='bold')
        ax.set_xlim(vmin, vmax); ax.set_ylim(0, 1.05)
        ax.set_yticks([])
        ax.set_xlabel(name, fontsize=11, fontweight='bold')
        for sp in ['top', 'right', 'left']:
            ax.spines[sp].set_visible(False)
    axes[0].legend(loc='upper right', fontsize=9, ncol=2)
    plt.suptitle('Physical Bounds: Valid vs Outlier Regions',
                 fontsize=13, fontweight='bold', y=1.0)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, 'fig_physical_bounds.png'),
                bbox_inches='tight')
    plt.close()


fig_hampel_mechanism()
fig_robust_vs_3sigma()
fig_physical_bounds()
print('Figures saved to', OUT_DIR)
