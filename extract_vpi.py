"""V_π 추출 — 파장 시프트 방법.

원리:
    V_π = FSR / (2 · |dλ_null/dV|)

단계:
    1) V=-2V 스펙트럼에서 prominent null 위치 식별 → FSR 계산
    2) 각 null을 V≤0 구간에서 parabolic fit으로 정밀 추적 → dλ/dV
    3) V_π 계산

순바이어스(+0.5V)는 forward 전류로 동작 모드 변하므로 제외.
Null 트래킹이 이웃 null로 점프할 위험은 윈도우 크기 < FSR/2로 방지.
"""
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

    dlam_dV = float(np.mean(slopes))   # nm/V
    if dlam_dV == 0:
        return {'fsr_nm': round(fsr, 4),
                'dlam_dV_pm_per_V': float('nan'), 'vpi_V': float('nan')}
    vpi = fsr / (2 * abs(dlam_dV))
    return {
        'fsr_nm': round(fsr, 4),
        'dlam_dV_pm_per_V': round(dlam_dV * 1000, 2),
        'vpi_V': round(vpi, 2),
    }


def _nan_result():
    return {'fsr_nm': float('nan'),
            'dlam_dV_pm_per_V': float('nan'),
            'vpi_V': float('nan')}
