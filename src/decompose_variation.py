"""Variation Decomposition — Xing et al. ACS Photonics 2023 영감.

각 (Wafer × Band) 그룹에서 다이 metric 값을 두 컴포넌트로 분해:

    observed(R, C)  =  systematic(R, C)  +  random(R, C)
                      (2D 다항식 fit)      (잔차)

- systematic : 웨이퍼 위치(Row, Col) 의 매끄러운 함수.  공정의 spatial
               gradient (lithography focus, etch uniformity, deposition
               thickness 분포 등) 로 해석.
- random     : systematic 으로 설명 안 되는 잔여.  측정 노이즈 + 진짜
               다이별 random 공정 변동의 합.

기존 `zscore_map.py` 는 둘을 통째로 합쳐서 보여주므로
"왜 이 다이가 outlier 인가?" 의 두 원인을 구분 못 함.

Reference:
    Y. Xing, J. Dong, U. Khan, W. Bogaerts,
    "Capturing the Effects of Spatial Process Variations in Silicon
     Photonic Circuits", ACS Photonics 10, 928 (2023).

────────────────────────────────────────────────────────────────────────
출력 그림 (per metric):
    decompose_<col>.png  :  3 열 × n_wafer_band 행
                            [observed | systematic | random]
────────────────────────────────────────────────────────────────────────
"""
import os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

from plot_common import ordered_groups


# 2D 다항식 fit 차수
#   1 = 평면 (a + b·R + c·C)
#   2 = 곡면 (1 + R + C + R² + R·C + C²)
# 14 다이로는 degree 2 가 적당 (parameter 6개, 자유도 8 남음)
POLY_DEGREE = 2


def _poly_features(rows, cols, degree=POLY_DEGREE):
    """(R, C) → 다항식 feature 행렬 X.  X·β = z 형태로 fit 가능."""
    feats = [np.ones_like(rows, dtype=float)]
    for d in range(1, degree + 1):
        for i in range(d + 1):
            j = d - i
            feats.append((rows ** i) * (cols ** j))
    return np.stack(feats, axis=1)


def fit_systematic(rows, cols, values, degree=POLY_DEGREE):
    """(R, C, z) → systematic trend + 잔차.

    반환:
        coeffs    : 다항식 계수 (β)
        trend     : 각 다이 위치에서의 systematic 값
        residual  : observed − trend
        r2        : 결정계수 (systematic 이 분산을 얼마나 설명하는가)
    """
    rows = np.asarray(rows, dtype=float)
    cols = np.asarray(cols, dtype=float)
    values = np.asarray(values, dtype=float)
    valid = ~np.isnan(values)
    if valid.sum() < (degree + 1) * (degree + 2) // 2 + 1:
        return None  # 다이 수 부족
    X = _poly_features(rows[valid], cols[valid], degree)
    beta, *_ = np.linalg.lstsq(X, values[valid], rcond=None)
    # 모든 다이에 대해 예측
    X_all = _poly_features(rows, cols, degree)
    trend = X_all @ beta
    residual = values - trend
    # R²
    ss_res = np.nansum(residual[valid] ** 2)
    ss_tot = np.nansum((values[valid] - values[valid].mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    return {
        'coeffs':   beta,
        'trend':    trend,
        'residual': residual,
        'r2':       float(r2),
    }


def plot_decomposition(df, value_col, label, save_path, degree=POLY_DEGREE):
    """3 열 × n행 패널: [observed | systematic | random]."""
    pairs = ordered_groups(df)
    n = len(pairs)
    if n == 0:
        return

    fig, axes = plt.subplots(n, 3, figsize=(11, 3.0 * n), dpi=130,
                             squeeze=False)
    fig.suptitle(
        f'Variation Decomposition — {label}\n'
        f'(observed = systematic spatial trend + random residual; '
        f'2D poly deg={degree})',
        fontsize=12, fontweight='bold', y=0.995,
    )

    # 색 스케일: observed 와 systematic 은 같은 범위, residual 은 ±max(|res|)
    all_vals = df[value_col].dropna().to_numpy()
    if len(all_vals) == 0:
        plt.close(fig); return
    obs_lo, obs_hi = float(all_vals.min()), float(all_vals.max())

    for row_idx, (w, b) in enumerate(pairs):
        sub = df[(df['Wafer'] == w) & (df['Band'] == b)]
        rows = sub['Row'].to_numpy(dtype=float)
        cols = sub['Col'].to_numpy(dtype=float)
        vals = sub[value_col].to_numpy(dtype=float)

        fit = fit_systematic(rows, cols, vals, degree)
        if fit is None:
            for k in range(3):
                axes[row_idx][k].set_title(f'{w}[{b}]: 다이 부족')
                axes[row_idx][k].axis('off')
            continue
        trend = fit['trend']
        residual = fit['residual']
        r2 = fit['r2']
        res_max = float(np.nanmax(np.abs(residual))) or 1.0

        # 각 패널: observed, systematic, random
        for col_idx, (z, title, cmap, vlim) in enumerate([
            (vals,     'observed',                   'turbo',  (obs_lo, obs_hi)),
            (trend,    f'systematic (R²={r2:.2f})',  'turbo',  (obs_lo, obs_hi)),
            (residual, 'random (residual)',          'RdBu_r', (-res_max, res_max)),
        ]):
            ax = axes[row_idx][col_idx]
            sc = ax.scatter(cols, rows, c=z, cmap=cmap,
                            vmin=vlim[0], vmax=vlim[1],
                            s=260, marker='s', ec='black', lw=0.5)
            # 숫자 라벨
            for xi, yi, zv in zip(cols, rows, z):
                if np.isnan(zv): continue
                txt_color = 'white' if (abs(zv - (vlim[0] + vlim[1]) / 2)
                                        > 0.5 * (vlim[1] - vlim[0])) else 'black'
                ax.text(xi, yi, f'{zv:+.1f}' if col_idx == 2 else f'{zv:.1f}',
                        ha='center', va='center', fontsize=6.5, color=txt_color)
            ax.set_xlim(-5.5, 5.5); ax.set_ylim(-5.5, 5.5)
            ax.set_aspect('equal'); ax.grid(alpha=0.25)
            if row_idx == 0:
                ax.set_title(title, fontsize=10, fontweight='bold')
            if col_idx == 0:
                ax.set_ylabel(f'{w}[{b}]', fontsize=10, fontweight='bold')
            ax.set_xticks([]); ax.set_yticks([])
            plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(save_path, bbox_inches='tight')
    plt.close(fig)
    print(f'Decomposition saved: {save_path}')


def plot_all(df, run_dir):
    """ER, IL, Vpi 모두 분해."""
    metrics = [
        ('ER_dB', 'ER'),
        ('IL_dB', 'IL'),
        ('Vpi_V', 'V_pi'),
    ]
    for col, label in metrics:
        out = os.path.join(run_dir, f'decompose_{col}.png')
        plot_decomposition(df, col, label, out)


def variation_summary(df):
    """그룹별 systematic / random 분산 비율을 표로 출력 (CLI 진단용)."""
    metrics = ['ER_dB', 'IL_dB', 'Vpi_V']
    rows = []
    for (w, b), sub in df.groupby(['Wafer', 'Band']):
        rec = {'Wafer': w, 'Band': b, 'n': len(sub)}
        for col in metrics:
            fit = fit_systematic(sub['Row'], sub['Col'], sub[col])
            if fit is None:
                rec[f'{col}_R2'] = float('nan')
                rec[f'{col}_random_std'] = float('nan')
                continue
            rec[f'{col}_R2'] = round(fit['r2'], 2)
            valid_res = fit['residual'][~np.isnan(sub[col].to_numpy())]
            rec[f'{col}_random_std'] = round(float(np.std(valid_res)), 3)
        rows.append(rec)
    import pandas as pd
    return pd.DataFrame(rows)
