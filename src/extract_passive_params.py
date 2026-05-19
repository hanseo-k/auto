"""Passive 파라미터 추출 — coupler split ratio, MZM-section propagation loss.

APL Photonics 2024 (Optical and geometric parameter extraction across 300mm
PIC wafers) 의 방법론에 영감 받음: 단일 MZI 측정에서 여러 파라미터 동시 추출.

────────────────────────────────────────────────────────────────────────
물리 모델
────────────────────────────────────────────────────────────────────────
비균형 MZI 의 transfer function:
    T(λ) = |a·exp(jφ₁) + b·exp(jφ₂)|²
         = a² + b² + 2ab·cos(Δφ)

여기서 a, b 는 두 팔의 amplitude (splitter + arm loss 포함).

극값:
    T_max = (a + b)²    (constructive, Δφ = 0)
    T_min = (a − b)²    (destructive, Δφ = π)

ER (linear) = T_max / T_min = ((a+b)/(a-b))²

amplitude 비 k = b/a (≤ 1) 로 두면:
    √ER_lin = (1+k)/(1-k)
    k       = (√ER_lin − 1) / (√ER_lin + 1)

Power split ratio (50:50 이상적): k²
Imbalance (dB):                   −10·log₁₀(k²)  = −20·log₁₀(k)

────────────────────────────────────────────────────────────────────────
ALIGN 차감 후 의미
────────────────────────────────────────────────────────────────────────
T_dev = IL_mzm − IL_ref  →  ALIGN 의 커플링/도파로 손실 제거
                            남는 건 MZM section 의 loss + ER 패턴

T_dev_peak (가장 큰 값) ≈ −(MZM section 손실, dB)
    이상적 lossless MZM 이면 T_dev_peak = 0
    실제는 MMI splitter + combiner + phase shifter arm 손실의 합

→ propagation_loss_dB = −T_dev_peak  (양수, dB)
"""
import numpy as np


def coupler_split_ratio_from_er(er_dB):
    """ER (dB) → amplitude 비 k, power 비 k², imbalance (dB).

    가정:
      - MZI 양 쪽 splitter/combiner 균형이 같다 (대칭 구조)
      - arm loss 가 두 팔에서 같다 (or 둘이 합쳐서 a, b 안에 흡수됨)
      - ER 이 splitter imbalance 에만 의해 한정된다 (=  measurement noise floor 가 아니다)
    """
    if er_dB is None or np.isnan(er_dB) or er_dB <= 0:
        return {'amplitude_ratio_k': float('nan'),
                'power_split_ratio': float('nan'),
                'imbalance_dB': float('nan')}
    er_lin = 10 ** (er_dB / 10)
    sqrt_er = np.sqrt(er_lin)
    k = (sqrt_er - 1) / (sqrt_er + 1)       # b/a
    power_ratio = k ** 2                     # b²/a²
    imbalance_dB = -20 * np.log10(k) if k > 0 else float('inf')
    return {
        'amplitude_ratio_k':  round(float(k), 5),
        'power_split_ratio':  round(float(power_ratio), 5),
        'imbalance_dB':       round(float(imbalance_dB), 2),
    }


def extract_propagation_loss(die):
    """T_dev_peak 의 최댓값 → MZM section 손실 (dB).

    ALIGN 레퍼런스 차감 후의 transfer function 에서 모든 바이어스 중
    최대 transmission 을 찾고, 그게 0 보다 낮은 정도를 손실로 본다.

    음수면 손실, 양수면 (드물지만) 캘리브레이션 오차.
    """
    from xml_loader import t_dev
    ref_L, ref_IL = die['ref_L'], die['ref_IL']
    band = die['band']
    # ER 윈도우와 동일한 범위에서 평가 (일관성)
    from extract_er import ER_WINDOW_NM
    if band not in ER_WINDOW_NM:
        return float('nan')
    lo, hi = ER_WINDOW_NM[band]

    peak_overall = -np.inf
    for V, (L, IL_mzm) in die['sweeps'].items():
        T = t_dev(L, IL_mzm, ref_L, ref_IL)
        mask = (L >= lo) & (L <= hi)
        if mask.sum() == 0:
            continue
        peak_overall = max(peak_overall, float(T[mask].max()))
    if peak_overall == -np.inf:
        return float('nan')
    # 손실은 양수로 보고 (T_dev_peak = -0.5 dB → 손실 0.5 dB)
    return round(-peak_overall, 3)


def coupler_imbalance_per_bias(die):
    """각 바이어스에서의 splitter imbalance (dB) 를 dict 로 반환.

    각 바이어스의 스펙트럼 안에서 peak-null 차이 = bias-별 ER 을 계산하고,
    그로부터 imbalance_dB 를 동일 공식으로 유도한다.

    반환: { V_bias (float) : imbalance_dB (float) }
    """
    from xml_loader import t_dev
    from extract_er import ER_WINDOW_NM

    band = die['band']
    if band not in ER_WINDOW_NM:
        return {}
    lo, hi = ER_WINDOW_NM[band]
    ref_L, ref_IL = die['ref_L'], die['ref_IL']

    result = {}
    for V, (L, IL_mzm) in die['sweeps'].items():
        T = t_dev(L, IL_mzm, ref_L, ref_IL)
        mask = (L >= lo) & (L <= hi)
        if mask.sum() == 0:
            result[float(V)] = float('nan')
            continue
        er_dB = float(T[mask].max() - T[mask].min())
        if er_dB <= 0 or np.isnan(er_dB):
            result[float(V)] = float('nan')
            continue
        er_lin = 10 ** (er_dB / 10)
        sqrt_er = np.sqrt(er_lin)
        k = (sqrt_er - 1) / (sqrt_er + 1)
        if k <= 0:
            result[float(V)] = float('nan')
            continue
        result[float(V)] = round(float(-20 * np.log10(k)), 3)
    return result


def extract_passive_params(die, er_dB):
    """ER 값 + 다이 spectrum 으로부터 passive 파라미터 dict.

    반환:
        amplitude_ratio_k       : b/a, ≤ 1
        power_split_ratio       : k², 1 이면 완벽 50:50
        imbalance_dB            : 전체 ER 기반 imbalance, 0 dB 이상적
        mzm_loss_dB             : MZM section 손실 (양수, dB)
        imbalance_per_bias_dB   : { V_bias : imbalance_dB }  (바이어스별 상세)
    """
    out = coupler_split_ratio_from_er(er_dB)
    out['mzm_loss_dB'] = extract_propagation_loss(die)
    out['imbalance_per_bias_dB'] = coupler_imbalance_per_bias(die)
    return out
