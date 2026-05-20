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
반환 dict 의 'vpi_status' — 명시적 상태 플래그 (꼬리표)
──────────────────────────────────────────────────────────────────────
    'ok'           : V_π 정상 추출
    'slope_filter' : |dλ/dV| < MIN_SLOPE_PM_PER_V — 측정 망가짐 (V_π=NaN)
    'no_nulls'     : deep null < 2 개 — FSR 추출 실패 (모두 NaN)
    'no_slopes'    : 모든 null tracking 이 점프로 reject (V_π=NaN)
    'few_biases'   : reverse-bias 데이터 부족 (3개 미만)
    'no_sweeps'    : sweep 데이터 자체가 없음
이 플래그는 CSV/XLSX 생성 모듈에서 사람 친화적 reason 으로 변환됨.

──────────────────────────────────────────────────────────────────────
MIN_SLOPE_PM_PER_V (slope filter) — 망가진 측정 검출 근거
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

import numpy as np
from scipy.signal import find_peaks


MIN_SLOPE_PM_PER_V = 10.0   # |dλ/dV| 이 이 값 미만이면 측정 망가진 것으로 간주

# 사람 친화적 reason 메시지 매핑 (CSV/XLSX 모듈에서 활용)
STATUS_MESSAGE = {
    'ok':           '',
    'slope_filter': 'broken: |dλ/dV| < {thr:.0f} pm/V (slope filter, |dλ/dV|={dl:.2f})',
    'no_nulls':     'broken: deep null < 2 (FSR 추출 실패)',
    'no_slopes':    'broken: 모든 null tracking jump-reject',
    'few_biases':   'broken: reverse-bias 데이터 < 3개',
    'no_sweeps':    'broken: sweep 데이터 없음',
}


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
    """반환: dict { 'fsr_nm', 'dlam_dV_pm_per_V', 'vpi_V', 'vpi_status' }

    vpi_status 가 'ok' 가 아닌 경우 vpi_V 는 NaN.
    dλ/dV 와 FSR 은 가능한 한 보존 (진단용).
    """
    sweeps = die['sweeps']
    biases = sorted(sweeps.keys())
    if len(biases) < 3:
        return _result(status='no_sweeps')

    L0, I0 = sweeps[biases[0]]
    distance = max(1, int(1.0 / (L0[1] - L0[0])))
    peaks, _ = find_peaks(-I0, prominence=10, distance=distance)
    deep = sorted([float(L0[p]) for p in peaks if I0[p] < -25])
    if len(deep) < 2:
        return _result(status='no_nulls')
    fsr = float(np.median(np.diff(deep)))

    rev_biases = [v for v in biases if v <= 0.0]
    if len(rev_biases) < 3:
        return _result(fsr=fsr, status='few_biases')

    half = min(0.4, fsr * 0.35)
    slopes = []
    r2_values = []   # 각 null fit 의 R² (선형성 지표)
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
        vb = np.array(rev_biases)[valid]
        yy = positions[valid]
        s, b = np.polyfit(vb, yy, 1)
        slopes.append(s)
        # R² (선형성): 1 = 완벽 직선, 0 = 무관계
        pred = s * vb + b
        ss_res = np.sum((yy - pred) ** 2)
        ss_tot = np.sum((yy - yy.mean()) ** 2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')
        r2_values.append(r2)

    if not slopes:
        return _result(fsr=fsr, status='no_slopes')

    slopes = np.array(slopes)
    r2_arr = np.array(r2_values)
    if len(slopes) >= 3:
        med = np.median(slopes)
        mad = np.median(np.abs(slopes - med))
        keep = np.abs(slopes - med) <= 3 * mad + 1e-6
        slopes = slopes[keep]
        r2_arr = r2_arr[keep]

    dlam_dV = float(np.mean(slopes))            # nm/V
    dlam_dV_pm = dlam_dV * 1000                  # pm/V (보고용)

    # 선형성 R² (살아남은 null 들의 median — 한 null 이 우연히 좋게/나쁘게 나오는 효과 완화)
    linearity_r2 = float(np.nanmedian(r2_arr)) if len(r2_arr) > 0 else float('nan')

    # Slope filter: 너무 작으면 측정 망가진 것 → V_π 폭주 방지
    if abs(dlam_dV_pm) < MIN_SLOPE_PM_PER_V:
        return _result(fsr=fsr, dlam_dV_pm=dlam_dV_pm,
                       r2=linearity_r2, status='slope_filter')

    vpi = fsr / (2 * abs(dlam_dV))
    return _result(fsr=fsr, dlam_dV_pm=dlam_dV_pm, vpi=vpi,
                   r2=linearity_r2, status='ok')


def _result(fsr=None, dlam_dV_pm=None, vpi=None, r2=None, status='ok'):
    """반환 dict 만드는 helper."""
    return {
        'fsr_nm':           round(fsr, 4) if fsr is not None else float('nan'),
        'dlam_dV_pm_per_V': round(dlam_dV_pm, 2) if dlam_dV_pm is not None else float('nan'),
        'vpi_V':            round(vpi, 2) if vpi is not None else float('nan'),
        'linearity_R2':     round(r2, 4) if r2 is not None and not np.isnan(r2) else float('nan'),
        'vpi_status':       status,
    }


def status_to_reason(status, dlam_dV_pm=None):
    """vpi_status → 사람 친화적 reason 문자열."""
    if status == 'ok':
        return ''
    msg = STATUS_MESSAGE.get(status, f'unknown status: {status}')
    if status == 'slope_filter':
        return msg.format(thr=MIN_SLOPE_PM_PER_V,
                          dl=abs(dlam_dV_pm) if dlam_dV_pm is not None else float('nan'))
    return msg
