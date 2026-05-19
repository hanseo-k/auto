"""신뢰도 맵 — 각 다이의 Robust Z-score를 색으로.

Robust Z = (x − median) / (1.4826 × MAD)   (per Wafer-Band 그룹)
|Z| > 3 인 다이 = outlier (3-sigma의 robust 등가)

색: 발산형 cmap (median=0 중심, 음/양 양쪽으로)
숫자: |z| 값 표기
"""
import os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_trust_map(df, value_col, label, save_path):
    z_col = f'robust_z_{value_col}'
    out_col = f'is_outlier_{value_col}'
    if z_col not in df.columns:
        print(f'[skip] {z_col} not found')
        return

    pairs = sorted({(w, b) for w, b in zip(df['Wafer'], df['Band'])},
                   key=lambda x: (x[1] != 'C', x[0]))
    n = len(pairs)
    fig, axes = plt.subplots(1, n, figsize=(3.6 * n, 4.2), dpi=140)
    if n == 1:
        axes = [axes]
    fig.suptitle(f'Robust Z-score Map  —  {label}\n'
                 f'(|z| > 3  →  outlier;   z = (x − median) / (1.4826 × MAD))',
                 fontsize=12, fontweight='bold', y=1.04)

    last = None
    for ax, (w, b) in zip(axes, pairs):
        sub = df[(df['Wafer'] == w) & (df['Band'] == b)]
        x = sub['Col'].to_numpy(dtype=float)
        y = sub['Row'].to_numpy(dtype=float)
        rz = sub[z_col].to_numpy(dtype=float)
        is_out = sub[out_col].to_numpy(dtype=bool)

        # 발산형 컬러: median=0 기준, ±5 범위로 클립
        sc = ax.scatter(x, y, c=rz, cmap='RdBu_r', vmin=-5, vmax=5,
                        s=600, marker='s', ec='black', lw=0.8)
        last = sc

        # 숫자 라벨 (|z| 표시)
        for xi, yi, zv, ov in zip(x, y, rz, is_out):
            if np.isnan(zv):
                continue
            text_c = 'white' if abs(zv) > 2.5 else 'black'
            extra = ' !' if ov else ''
            ax.text(xi, yi, f'{zv:+.1f}{extra}', ha='center', va='center',
                    fontsize=9, fontweight='bold', color=text_c)

        ax.set_xlim(-5.5, 5.5); ax.set_ylim(-5.5, 5.5)
        ax.set_aspect('equal')
        ax.set_xticks(range(-5, 6)); ax.set_yticks(range(-5, 6))
        ax.set_xlabel('Column'); ax.set_ylabel('Row')
        ax.set_title(f'{w}  [{b}-band]', fontsize=10, fontweight='bold')
        ax.grid(alpha=0.25)

    cbar = fig.colorbar(last, ax=axes, shrink=0.85, pad=0.02,
                        ticks=range(-5, 6))
    cbar.set_label('Robust Z-score (median = 0)\n'
                   '|z| > 3 = outlier  ("!" mark)')
    plt.savefig(save_path, bbox_inches='tight')
    plt.close(fig)
    print(f'Trust map saved: {save_path}')


def plot_all(df, run_dir):
    metrics = [
        ('ER_dB', 'ER'),
        ('IL_dB', 'IL'),
        ('Vpi_V', 'V_pi'),
    ]
    for col, label in metrics:
        out = os.path.join(run_dir, f'trust_map_{col}.png')
        plot_trust_map(df, col, label, out)
