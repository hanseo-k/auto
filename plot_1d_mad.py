"""1D 분포 그래프 (MAD 기반 박스) — ER · IL · Vpi 각각.

기존 plot_1d.py의 IQR 박스 대신 MAD 기반 박스로 교체.
  - 박스:   median ± 1.4826×MAD  (±1σ 등가)
  - 수염:   median ± 3×1.4826×MAD (±3σ 등가 = outlier 판정 경계)
  → 수염 밖 = outlier와 일치하므로 산점도 속 빈 동그라미가 항상 수염 밖에 위치
"""
import os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.lines import Line2D

from outlier_detect import PHYSICAL_BOUNDS


WAFER_BAND_COLOR = {
    ('D07', 'C'): '#4C72B0', ('D08', 'C'): '#7AA0CB',
    ('D08', 'O'): '#DD8452', ('D23', 'O'): '#E8A87C', ('D24', 'O'): '#F4B999',
}


def _ordered_groups(df):
    pairs = sorted({(w, b) for w, b in zip(df['Wafer'], df['Band'])},
                   key=lambda x: (x[1] != 'C', x[0]))
    return pairs


def _mad_stats(values):
    """median, sigma_robust(=1.4826×MAD) 반환. 값 3개 미만이면 None."""
    v = values[~np.isnan(values)]
    if len(v) < 3:
        return None
    med = float(np.median(v))
    mad = float(np.median(np.abs(v - med)))
    sigma = 1.4826 * mad
    return med, sigma


def _draw_mad_box(ax, x_center, med, sigma, color, width=0.4):
    """MAD 기반 박스+수염을 직접 그림."""
    box_lo = med - sigma        # ±1σ
    box_hi = med + sigma
    whi_lo = med - 3 * sigma    # ±3σ (outlier 경계)
    whi_hi = med + 3 * sigma
    hw = width / 2

    # 박스
    rect = Rectangle((x_center - hw, box_lo), width, box_hi - box_lo,
                      facecolor=color, alpha=0.55, edgecolor='black', lw=1.2, zorder=2)
    ax.add_patch(rect)

    # median 선
    ax.plot([x_center - hw, x_center + hw], [med, med],
            color='black', lw=2.0, zorder=3)

    # 수염 (세로선)
    ax.plot([x_center, x_center], [box_hi, whi_hi], color='black', lw=1.2, zorder=2)
    ax.plot([x_center, x_center], [whi_lo, box_lo], color='black', lw=1.2, zorder=2)

    # 수염 끝 가로선 (cap)
    cap_w = hw * 0.6
    ax.plot([x_center - cap_w, x_center + cap_w], [whi_hi, whi_hi],
            color='black', lw=1.2, zorder=2)
    ax.plot([x_center - cap_w, x_center + cap_w], [whi_lo, whi_lo],
            color='black', lw=1.2, zorder=2)


def plot_1d_mad(df, value_col, label, save_path):
    out_col = f'is_outlier_{value_col}'
    lo, hi = PHYSICAL_BOUNDS[value_col]
    groups = _ordered_groups(df)

    fig, ax = plt.subplots(figsize=(11, 6), dpi=140)

    # 물리 신뢰 영역
    ax.axhspan(lo, hi, color='lightgreen', alpha=0.18, zorder=0,
               label=f'physical bound [{lo}, {hi}]')

    rng = np.random.default_rng(42)
    for i, (w, b) in enumerate(groups):
        x_pos = i + 1
        sub = df[(df['Wafer'] == w) & (df['Band'] == b)]
        color = WAFER_BAND_COLOR.get((w, b), '#888')

        trusted_vals = sub[~sub[out_col]][value_col].to_numpy(dtype=float)
        stats = _mad_stats(trusted_vals)
        if stats is not None:
            med, sigma = stats
            _draw_mad_box(ax, x_pos, med, sigma, color)

        # 산점도
        z = sub[value_col].to_numpy(dtype=float)
        out = sub[out_col].to_numpy(dtype=bool)
        x = rng.normal(x_pos, 0.06, size=len(z))
        m_trust = ~out & ~np.isnan(z)
        ax.scatter(x[m_trust], z[m_trust], facecolor=color, edgecolor='black',
                   s=42, lw=0.5, alpha=0.95, zorder=4)
        m_out = out & ~np.isnan(z)
        ax.scatter(x[m_out], z[m_out], facecolor='none', edgecolor=color,
                   s=70, lw=1.6, alpha=0.95, zorder=5)

    # 밴드 경계선
    band_seq = [b for _, b in groups]
    if 'C' in band_seq and 'O' in band_seq:
        n_C = band_seq.count('C')
        ax.axvline(n_C + 0.5, color='gray', ls='--', alpha=0.5, lw=1)

    ax.set_xticks(range(1, len(groups) + 1))
    ax.set_xticklabels([f'{w}\n[{b}]' for w, b in groups])
    ax.set_ylabel(label, fontsize=12)
    ax.set_title(f'{label}  +MAD  (box=±1σ, whisker=±3σ, ○=outlier)',
                 fontsize=13, pad=10)
    ax.grid(alpha=0.3, axis='y')

    legend_handles = [
        Line2D([0], [0], color='lightgreen', linewidth=8, alpha=0.4,
               label=f'physical bound [{lo}, {hi}]'),
        Line2D([0], [0], marker='o', color='gray', markerfacecolor='gray',
               markersize=7, label='trusted', linestyle='None'),
        Line2D([0], [0], marker='o', color='gray', markerfacecolor='none',
               markersize=9, markeredgewidth=1.6, label='outlier', linestyle='None'),
    ]
    ax.legend(handles=legend_handles, loc='best', fontsize=9, framealpha=0.85)

    plt.tight_layout()
    plt.savefig(save_path, bbox_inches='tight')
    plt.close(fig)
    print(f'1D+MAD plot saved: {save_path}')


def plot_all(df, run_dir):
    metrics = [
        ('ER_dB', 'Extinction Ratio (dB)'),
        ('IL_dB', 'Insertion Loss (dB)'),
        ('Vpi_V', 'V_pi (V)'),
    ]
    for col, label in metrics:
        out = os.path.join(run_dir, f'1d_mad_{col}.png')
        plot_1d_mad(df, col, label, out)
