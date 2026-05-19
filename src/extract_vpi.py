"""V_π 추출 — 파장 시프트 방법.

원리:
    V_π = FSR / (2 · |dλ_null/dV|)

단계:
    1) V=-2V 스펙트럼에서 prominent null 위치 식별 → FSR 계산
    2) 각 null을 V≤0 구간에서 parabolic fit으로 정밀 추적 → dλ/dV
    3) Slope sanity check: |dλ/dV| < MIN_SLOPE_PM_PER_V 이면 NaN
    4) V_π 계산

순바이어스(+0.5V)는 forward 전류로 동작 모드 변하므로 제외.
Null 트래킹이 이웃 null로 점프할 위험은 윈도우 크기 < FSR/2로 방지.

──────────────────────────────────────────────────────────────────────
MIN_SLOPE_PM_PER_V (slope filter) — 망가진 측정 검출
──────────────────────────────────────────────────────────────────────
바이어스를 인가했는데도 null 파장이 거의 움직이지 않으면 (dλ/dV ≈ 0)
V_π = FSR / (2·|0|) → 무한대로 폭주.  이건 디바이스가 그런 게 아니라
측정 자체가 망가졌다는 신호 (probe contact 불량, 케이블 단선, SW 버그).

HY202103 의 정상 디바이스 dλ/dV 범위 (실측 검증):
  - C-band: -210 pm/V 부근
  - O-band:  -120 ~ -180 pm/V
  - 정상 측정의 |min| ≈ 107 pm/V

2019-05-31 의 망가진 28개 측정:
  - dλ/dV 의 |min| ≈ 0.06 pm/V (= 정상의 ~1/1700)
  - V_π 값이 1062 ~ 78633 V 로 폭주

→ 정상(107)과 망가짐(4.7)의 안전한 분기점으로 10 pm/V 채택.
  sensitivity test (`src/sensitivity_test.py`):
    - min_slope=10 → 망가진 28개 100% NaN, 정상 70개 0개 false positive
"""

MIN_SLOPE_PM_PER_V = 10.0   # |dλ/dV| 이 이 값 미만이면 측정 망가진 것으로 간주
import numpy as np
from scipy.signal import find_peaks


def _parabolic_null(L, IL, lam_guess, half):
    m = (L >= lam_guess - half) & (L <= lam_guess + half)
    if m.sum() < 5:
        return float('nan')
    Lw, Iw = L[m], IL[m]
    i = int(np.argmin(Iw))
    lo, hi = max(0, i - 5), min(len(Lw), i + 6)
    if hi - lo < 3:
        return float(Lw[i])
    a, b, _ = np.polyfit(Lw[lo:hi], Iw[lo:hi], 2)
    if a <= 0:
        return float(Lw[i])
    return float(-b / (2 * a))


def extract_vpi(die):
    """반환: dict { 'fsr_nm', 'dlam_dV_pm_per_V', 'vpi_V' }"""
    sweeps = die['sweeps']
    biases = sorted(sweeps.keys())
    if len(biases) < 3:
        return _nan_result()

    L0, I0 = sweeps[biases[0]]
    distance = max(1, int(1.0 / (L0[1] - L0[0])))
    peaks, _ = find_peaks(-I0, prominence=10, distance=distance)
    deep = sorted([float(L0[p]) for p in peaks if I0[p] < -25])
    if len(deep) < 2:
        return _nan_result()
    fsr = float(np.median(np.diff(deep)))

    rev_biases = [v for v in biases if v <= 0.0]
    if len(rev_biases) < 3:
        return {'fsr_nm': round(fsr, 4),
                'dlam_dV_pm_per_V': float('nan'), 'vpi_V': float('nan')}

    half = min(0.4, fsr * 0.35)
    slopes = []
    for lam0 in deep:
        positions = []
        for v in rev_biases:
            L, I = sweeps[v]
            positions.append(_parabolic_null(L, I, lam0, half))
        positions = np.array(positions)
        valid = ~np.isnan(positions)
        if valid.sum() < 3:
            continue
        # 점프 감지: 총 시프트가 윈도우보다 크면 reject
        total_shift = abs(positions[valid][-1] - positions[valid][0])
        if total_shift > half * 1.5:
            continue
        s, _ = np.polyfit(np.array(rev_biases)[valid], positions[valid], 1)
        slopes.append(s)

    if not slopes:
        return {'fsr_nm': round(fsr, 4),
                'dlam_dV_pm_per_V': float('nan'), 'vpi_V': float('nan')}

    slopes = np.array(slopes)
    if len(slopes) >= 3:
        med = np.median(slopes)
        mad = np.median(np.abs(slopes - med))
        slopes = slopes[np.abs(slopes - med) <= 3 * mad + 1e-6]

    dlam_dV = float(np.mean(slopes))            # nm/V
    dlam_dV_pm = dlam_dV * 1000                  # pm/V (보고용)

    # Slope filter: 너무 작으면 측정 망가진 것 → V_π 폭주 방지
    if abs(dlam_dV_pm) < MIN_SLOPE_PM_PER_V:
        return {'fsr_nm': round(fsr, 4),
                'dlam_dV_pm_per_V': round(dlam_dV_pm, 2),
                'vpi_V': float('nan')}

    vpi = fsr / (2 * abs(dlam_dV))
    return {
        'fsr_nm': round(fsr, 4),
        'dlam_dV_pm_per_V': round(dlam_dV_pm, 2),
        'vpi_V': round(vpi, 2),
    }


def _nan_result():
    return {'fsr_nm': float('nan'),
            'dlam_dV_pm_per_V': float('nan'),
            'vpi_V': float('nan')}
