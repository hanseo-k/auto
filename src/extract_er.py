"""ER (Extinction Ratio) 추출 — 고정 윈도우 기반.

정의:
    ER = peak T_dev − null T_dev   (밴드별 고정 윈도우 안, 모든 바이어스에서)

윈도우 (밴드별 고정, FSR ~ 1개 이상 포함):
    C-band:  1546 ~ 1560 nm
    O-band:  1306 ~ 1320 nm  (C-band과 동일 비례로 설정)

이전엔 λ_c ± half_window 방식이라 다이마다 윈도우가 살짝 달랐음.
이제 절대 좌표로 고정 → 다이 간 비교 일관성 ↑.
"""
import numpy as np
from xml_loader import t_dev


# 밴드별 고정 윈도우 (사용자 지정)
ER_WINDOW_NM = {
    'C': (1546.0, 1560.0),
    'O': (1306.0, 1320.0),
}


def extract_er(die):
    """die: xml_loader.load_die 출력 dict.
    밴드별 고정 윈도우 안에서 모든 바이어스 × 모든 파장의 T_dev 중
    최댓값(peak)과 최솟값(null)의 차이.
    """
    band = die['band']
    if band not in ER_WINDOW_NM:
        return float('nan')
    lo, hi = ER_WINDOW_NM[band]
    ref_L, ref_IL = die['ref_L'], die['ref_IL']

    peak_t = -np.inf
    null_t = np.inf
    for V, (L, IL_mzm) in die['sweeps'].items():
        T = t_dev(L, IL_mzm, ref_L, ref_IL)
        mask = (L >= lo) & (L <= hi)
        if mask.sum() == 0:
            continue
        peak_t = max(peak_t, float(T[mask].max()))
        null_t = min(null_t, float(T[mask].min()))
    if peak_t == -np.inf or null_t == np.inf:
        return float('nan')
    return round(peak_t - null_t, 3)
