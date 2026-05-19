"""ER (Extinction Ratio) 추출 — 고정 윈도우 기반.

정의:
    ER = peak T_dev − null T_dev   (밴드별 고정 윈도우 안, 모든 바이어스에서)

윈도우 (밴드별 고정, FSR 보다 약간 넓게 → 다이별로 최소 1개 null 안정 포함):
    C-band:  1545 ~ 1561 nm  (16 nm 폭, FSR ≈ 14.3 nm)
    O-band:  1306 ~ 1320 nm  (14 nm 폭, FSR ≈  9.8 nm → 1.4개 null)

이전엔 λ_c ± half_window 방식이라 다이마다 윈도우가 살짝 달랐음.
이제 절대 좌표로 고정 → 다이 간 비교 일관성 ↑.

──────────────────────────────────────────────────────────────────────
ER 윈도우 폭 결정 근거 (HY202103 sensitivity test 기반)
──────────────────────────────────────────────────────────────────────
C-band FSR ≈ 14.3 nm 인데 윈도우 14 nm 는 너무 좁아 다이별로 null 이
0.97 ~ 0.99 개 (즉 가끔 0 개만 포함) → ER 안정성 저하.

`src/sensitivity_test.py` 결과:
  | 윈도우 폭 | median  | std    |
  | 14 nm     | 36.87   | 1.47   |
  | 16 nm     | 37.11   | 1.24   ← std 16% 개선 ✅
  | 18 nm     | 37.31   | 1.40
  | 22 nm     | 38.03   | 1.36
  | 36 nm     | 39.53   | 2.48   ← 노이즈/엣지 효과로 폭증, outlier 발생

→ 16 nm 가 sweet spot.  너무 넓히면 ALIGN reference 가 신뢰성 떨어지는
   밴드 엣지를 포함해서 노이즈가 늘어남.

O-band 는 FSR 9.8 nm 이라 14 nm 윈도우에 이미 1.4 개 null 안정 포함
→ 변경 불필요 (16 nm 로 늘려도 std 1.34 → 1.30 으로 미미한 개선).
"""
import numpy as np
from xml_loader import t_dev


# 밴드별 고정 윈도우 (sensitivity test 검증)
ER_WINDOW_NM = {
    'C': (1545.0, 1561.0),    # 16 nm — C-band FSR ≈ 14.3 nm 보다 약간 넓게
    'O': (1306.0, 1320.0),    # 14 nm — O-band FSR ≈  9.8 nm, 이미 충분
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
