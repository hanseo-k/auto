"""Outlier 검출 — 물리 한계 + Robust Z-score.

3-sigma는 정규성·표본 크기 가정 때문에 안 씀.
대신:
  Layer 1: 물리적 한계 바운드 (신뢰성 가장 높음, n 무관)
  Layer 2: Robust Z-score (median/MAD 기반, per Wafer-Band 그룹)

각 다이마다 boolean 플래그 ('is_outlier')를 붙여 반환.
실제 데이터는 보존하고 시각화에서만 다르게 표시 (속이 빈 마커 등).
"""
import numpy as np
import pandas as pd  # noqa


# ──────────────────────────────────────────────────────────────────────
# PHYSICAL_BOUNDS — Si depletion MZM 의 물리적으로 가능한 값 범위.
#
# 디바이스 사양 (HY202103 가정):
#   - TE-mode Si depletion-type MZM
#   - single-arm drive (push-pull 아님)
#   - MMI splitter/combiner (single-stage)
#   - polarization filter / cascaded MZI / DC bias section 없음
#
# 자세한 근거와 인용은 README 의 "Physical Bounds" 섹션 참조.
# 요약:
#   ER  (Extinction Ratio):    5 ~ 35 dB
#       - 단일 MMI splitter 의 산업적 ER 상한 ≈ 30 dB (Witzens 2018, Table II)
#       - 35 dB 이상은 push-pull/cascaded 거나 측정 아티팩트
#       - 5 dB 미만은 working device 라 보기 어려움 (MMI 균형 깨짐/도파로 깨짐)
#
#   IL  (Insertion Loss, ON-state peak):   -15 ~ -1 dB
#       - 표준 Si MZM ON-state IL ≈ 4 ~ 10 dB (Reed 2010, Witzens 2018)
#       - -1 dB 미만(= 더 좋음)은 비물리적 (커플링/MMI 손실만 해도 ≥1dB)
#       - -15 dB 미만(= 더 나쁨)은 깨진 디바이스
#
#   Vpi (Half-wave voltage):   2 ~ 60 V
#       - Si depletion V_π·L = 1 ~ 3 V·cm  (Soref & Bennett 1987, Witzens 2018)
#       - 긴 phase shifter (5 mm)  → V_π·L=1 → V_π ≈ 2 V
#       - 짧은 phase shifter (0.5 mm) → V_π·L=3 → V_π ≈ 60 V
#       - 60 V 초과는 dλ/dV ≈ 0 (트래킹 실패) 의 산술적 폭주이지 실제 V_π 아님
#
# 인용 (full citation 은 README 참조):
#   [1] Soref & Bennett, IEEE JQE 23, 123 (1987)        — plasma dispersion 이론
#   [2] Reed et al., Nature Photonics 4, 518 (2010)     — Si modulator review
#   [3] Patel et al., Opt. Express 23, 14263 (2015)     — 41 GHz Si MZM, 측정값
#   [4] Witzens, Proc. IEEE 106, 2158 (2018)            — 종합 review (Table II)
# ──────────────────────────────────────────────────────────────────────
PHYSICAL_BOUNDS = {
    # 문헌 + HY202103 실측 분포(README "Physical Bounds — Empirical Validation"
    # 참조)를 함께 고려한 범위. 문헌 표준보다 ER 상한이 다소 넓은 이유는
    # 우리 ER 정의(모든 바이어스에서의 peak−null)가 fixed-bias ER 보다
    # 자연스럽게 크게 나오기 때문 + HY202103 실측이 ~37 dB 영역에 있음.
    'ER_dB':  (10.0, 45.0),    # working device 하한 10 dB / 측정 아티팩트 상한 45 dB
    'IL_dB':  (-15.0, -1.0),   # ON-state peak transmission (fiber-to-fiber, dB)
    'Vpi_V':  (2.0,  60.0),    # V_π·L ∈ [1, 3] V·cm × L ∈ [0.5, 5] mm
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

    웨이퍼당 14다이로는 지역적(이웃) 비교 의미가 적으므로,
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
