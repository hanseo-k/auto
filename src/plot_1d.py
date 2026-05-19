"""1D 분포 그래프 — ER · IL · Vpi 각각.

특징:
  - 4 웨이퍼를 같은 패널에 나란히 (y축 동일)
  - 신뢰 다이: 채운 마커
  - Outlier 다이: 속이 빈 동그라미 (데이터에서 빼지 않음)
  - 물리적 신뢰 구간: 반투명 네모 박스로 배경에 표시
  - 박스플롯은 신뢰 데이터만으로 그려서 통계 왜곡 방지
"""
import os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from outlier_detect import PHYSICAL_BOUNDS


WAFER_BAND_COLOR = {
    ('D07', 'C'): '#4C72B0', ('D08', 'C'): '#7AA0CB',
    ('D08', 'O'): '#DD8452', ('D23', 'O'): '#E8A87C', ('D24', 'O'): '#F4B999',
}


def _ordered_groups(df):
    """C-band 먼저, 그 다음 O-band. 각 밴드 안에서 wafer 이름 정렬."""
    pairs = sorted({(w, b) for w, b in zip(df['Wafer'], df['Band'])},
                   key=lambda x: (x[1] != 'C', x[0]))
    return pairs


def plot_1d(df, value_col, label, save_path):
    out_col = f'is_outlier_{value_col}'
    lo, hi = PHYSICAL_BOUNDS[value_col]
    groups = _ordered_groups(df)

    fig, ax = plt.subplots(figsize=(11, 6), dpi=140)

    # 신뢰 영역 (반투명 박스)
    ax.axhspan(lo, hi, color='lightgreen', alpha=0.18, zorder=0,
               label=f'physical bound [{lo}, {hi}]')

    box_data, labels, colors = [], [], []
    for w, b in groups:
        sub = df[(df['Wafer'] == w) & (df['Band'] == b)]
        trusted = sub[~sub[out_col]][value_col].dropna().to_numpy(dtype=float)
        box_data.append(trusted)
        labels.append(f'{w}\n[{b}]')
        colors.append(WAFER_BAND_COLOR.get((w, b), '#888'))

    bp = ax.boxplot(box_data, tick_labels=labels, patch_artist=True,
                    widths=0.5, showmeans=True,
                    meanprops=dict(marker='D', mfc='white', mec='black', ms=6),
                    zorder=2)
    for patch, c in zip(bp['boxes'], colors):
        patch.set_facecolor(c); patch.set_alpha(0.55)

    # 산점도 — 모든 다이 (신뢰는 채움, outlier는 속 빈)
    z_col = f'robust_z_{value_col}'
    rng = np.random.default_rng(42)
    for i, (w, b) in enumerate(groups):
        sub = df[(df['Wafer'] == w) & (df['Band'] == b)]
        z = sub[value_col].to_numpy(dtype=float)
        out = sub[out_col].to_numpy(dtype=bool)
        rz = (sub[z_col].to_numpy(dtype=float)
              if z_col in sub else np.full_like(z, np.nan))
        x = rng.normal(i + 1, 0.06, size=len(z))
        c = colors[i]
        m = ~out & ~np.isnan(z)
        ax.scatter(x[m], z[m], facecolor=c, edgecolor='black', s=42,
                   lw=0.5, alpha=0.95, zorder=4)
        m = out & ~np.isnan(z)
        ax.scatter(x[m], z[m], facecolor='none', edgecolor=c, s=70,
                   lw=1.6, alpha=0.95, zorder=5)

    # 밴드 경계선만 (라벨은 x축 틱에 [C]/[O]로 이미 표시됨)
    band_seq = [b for _, b in groups]
    if 'C' in band_seq and 'O' in band_seq:
        n_C = band_seq.count('C')
        ax.axvline(n_C + 0.5, color='gray', ls='--', alpha=0.5, lw=1)

    ax.set_ylabel(label, fontsize=12)
    ax.set_title(f'{label}  (○ = outlier, green box = physical bound)',
                 fontsize=13, pad=10)
    ax.grid(alpha=0.3, axis='y')
    ax.legend(loc='best', fontsize=9, framealpha=0.85)

    plt.tight_layout()
    plt.savefig(save_path, bbox_inches='tight')
    plt.close(fig)
    print(f'1D plot saved: {save_path}')


def plot_all(df, run_dir):
    metrics = [
        ('ER_dB', 'Extinction Ratio (dB)'),
        ('IL_dB', 'Insertion Loss (dB)'),
        ('Vpi_V', 'V_pi (V)'),
    ]
    for col, label in metrics:
        out = os.path.join(run_dir, f'1d_{col}.png')
        plot_1d(df, col, label, out)
