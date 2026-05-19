"""웨이퍼맵 — 연속 surface (Delaunay 보간 컨투어).

다이 좌표(Row,Col)에서 Delaunay 삼각분할로 부드러운 색 채움.
4 웨이퍼 1×4 패널, 같은 metric은 컬러 스케일 통일.
검은 점 = 실측 위치. Outlier 다이는 속이 빈 마커로 구분.
"""
import os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.ticker import MaxNLocator


def plot_wafer_map(df, value_col, label, save_path):
    """연속 surface 웨이퍼맵. (Wafer, Band)별 한 패널씩.
    같은 다이가 두 밴드로 측정됐으면 (D08-C, D08-O) 별도 패널."""
    pairs = sorted({(w, b) for w, b in zip(df['Wafer'], df['Band'])},
                   key=lambda x: (x[1] != 'C', x[0]))
    vals = df[value_col].dropna()
    if vals.empty:
        print(f'[skip] no values for {value_col}')
        return
    # 컬러 범위: 실제 데이터 min/max에 딱 맞춤 (그라데이션 보존)
    lo, hi = float(vals.min()), float(vals.max())
    vlim = (lo, hi)

    n = len(pairs)
    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 5), dpi=140)
    if n == 1:
        axes = [axes]
    fig.suptitle(f'Wafer Map (continuous surface) — {label}',
                 fontsize=14, fontweight='bold', y=1.02)
    out_col = f'is_outlier_{value_col}'
    last = None
    for ax, (w, b) in zip(axes, pairs):
        sub = df[(df['Wafer'] == w) & (df['Band'] == b)]
        last = _panel(ax, sub, value_col, out_col, vlim,
                      f'{w}  [{b}-band]')
    if last is not None:
        cbar = fig.colorbar(last, ax=axes, shrink=0.85, pad=0.02)
        cbar.set_label(label)
        cbar.locator = MaxNLocator(integer=True)
        cbar.update_ticks()
    plt.savefig(save_path, bbox_inches='tight')
    plt.close(fig)
    print(f'Wafer map saved: {save_path}')


def _panel(ax, sub, col, out_col, vlim, title):
    x = sub['Col'].to_numpy(dtype=float)
    y = sub['Row'].to_numpy(dtype=float)
    z = sub[col].to_numpy(dtype=float)
    is_out = sub[out_col].to_numpy(dtype=bool) if out_col in sub else \
             np.zeros_like(z, dtype=bool)
    valid = ~np.isnan(z)
    x, y, z, is_out = x[valid], y[valid], z[valid], is_out[valid]
    if len(x) < 3:
        ax.set_title(title); return None

    tri = mtri.Triangulation(x, y)
    levels = np.linspace(vlim[0], vlim[1], 20)
    tcf = ax.tricontourf(tri, z, levels=levels, cmap='turbo',
                         vmin=vlim[0], vmax=vlim[1], extend='both')
    ax.tricontour(tri, z, levels=levels, colors='white',
                  linewidths=0.3, alpha=0.4)

    # 신뢰 다이: 채운 검은 점 / outlier: 속 빈 동그라미
    trust = ~is_out
    if trust.any():
        ax.scatter(x[trust], y[trust], c='black', s=24, zorder=5,
                   ec='white', lw=0.6, label='trusted')
    if is_out.any():
        ax.scatter(x[is_out], y[is_out], facecolors='none',
                   edgecolors='black', s=80, lw=1.8, zorder=6,
                   label='outlier')

    ax.set_xlim(-6, 6); ax.set_ylim(-6, 6)
    ax.set_aspect('equal')
    ax.set_xlabel('Column'); ax.set_ylabel('Row')
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.grid(alpha=0.2)
    if is_out.any():
        ax.legend(fontsize=8, loc='upper right', framealpha=0.85)
    return tcf


def plot_all(df, run_dir):
    """ER, IL, Vpi 모두 그리기."""
    metrics = [
        ('ER_dB', 'Extinction Ratio (dB)'),
        ('IL_dB', 'Insertion Loss (dB)'),
        ('Vpi_V', 'V_pi (V)'),
    ]
    for col, label in metrics:
        out = os.path.join(run_dir, f'wafer_map_{col}.png')
        plot_wafer_map(df, col, label, out)
