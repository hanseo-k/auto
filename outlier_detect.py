"""Outlier 검출 — 물리 한계 + Hampel filter (이웃 비교).

3-sigma는 정규성·표본 크기 가정 때문에 안 씀.
대신:
  Layer 1: 물리적 한계 바운드 (신뢰성 가장 높음, n 무관)
  Layer 2: Hampel filter (공간 이웃과의 robust 비교)

각 다이마다 boolean 플래그 ('is_outlier')를 붙여 반환.
실제 데이터는 보존하고 시각화에서만 다르게 표시 (속이 빈 마커 등).
"""
import numpy as np
import pandas as pd  # noqa


# ──────────────────────────────────────────────
# Si MZM 물리 한계 — 우리 디바이스 타입 기준
#   디바이스 사양: standard TE-mode Si depletion MZM,
#                  single-arm drive, MMI splitter, single-stage,
#                  no polarization filter / no push-pull / no cascaded
#
# Reference:
#   Patel et al., Opt. Express 23, 14263 (2015) — 41GHz Si MZM
#   Witzens, Proc. IEEE 106, 2158 (2018) — 종합 review
#   Soref & Bennett, IEEE JQE 23, 123 (1987) — plasma dispersion
# ──────────────────────────────────────────────
PHYSICAL_BOUNDS = {
    'ER_dB':  (0.0,   45.0),   # standard TE Si depletion MZM with MMI splitter
    'IL_dB':  (-20.0,  0.0),   # 표준 IL (peak transmission) 운용 한계
    'Vpi_V':  (5.0,   80.0),   # V_π·L 1~3 V·cm × 짧은 phase shifter
}


def physical_outlier(values, col):
    """값이 PHYSICAL_BOUNDS 밖이면 True."""
    lo, hi = PHYSICAL_BOUNDS[col]
    v = pd.Series(values)
    return ((v < lo) | (v > hi)).fillna(True)


def robust_z_outlier(df, col, k=3.0):
    """Per-wafer-band Robust Z-score (Modified Z-score) outlier detection.

    공식:
        σ_robust = 1.4826 × MAD            (정규분포에서 σ 등가)
        z'       = (x − median) / σ_robust  (3-sigma의 robust 버전)
        outlier  if  |z'| > k               (k=3 이면 3-sigma 동등 임계)

    웨이퍼당 14다이로는 지역적 비교(Hampel) 의미가 적으므로,
    같은 (Wafer, Band) 그룹 전체의 median/MAD를 reference로 사용.

    Reference: Iglewicz & Hoaglin, "How to Detect and Handle Outliers",
               ASQC Quality Press, 1993.

    반환: (is_outlier_series, robust_z_series)
    """
    out = pd.Series(False, index=df.index)
    z_score = pd.Series(np.nan, index=df.index, dtype=float)
    for (wafer, band), grp in df.groupby(['Wafer', 'Band']):
        vals = grp[col].dropna()
        if len(vals) < 3:
            continue
        m = vals.median()
        mad = (vals - m).abs().median()
        sigma_robust = 1.4826 * mad
        if sigma_robust == 0:
            continue
        for idx in grp.index:
            x = grp.loc[idx, col]
            if pd.isna(x):
                out.at[idx] = True
                continue
            z = (x - m) / sigma_robust
            z_score.at[idx] = z
            if abs(z) > k:
                out.at[idx] = True
    return out, z_score


def mark_outliers(df, k_threshold=3.0):
    """ER, IL, Vpi 각각에 대해 outlier 플래그 컬럼 생성.

    k_threshold: Robust Z-score 임계값 (기본 3.0 = 3-sigma 동등).
                 더 관대하게 하려면 3.5, 4.0 등으로.

    각 컬럼: 물리 바운드 OR Robust Z 둘 중 하나라도 outlier면 True.
    """
    df = df.copy()
    for col in ['ER_dB', 'IL_dB', 'Vpi_V']:
        phys = physical_outlier(df[col], col)
        rz, z_val = robust_z_outlier(df, col, k=k_threshold)
        df[f'is_outlier_{col}'] = phys | rz
        df[f'robust_z_{col}'] = z_val.round(2)
    df['is_trusted'] = ~(
        df['is_outlier_ER_dB'] |
        df['is_outlier_IL_dB'] |
        df['is_outlier_Vpi_V']
    )
    return df
