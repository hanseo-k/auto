"""IL (Insertion Loss) — V=-1V 고정, 윈도우 내 peak transmission.
   ALIGN reference 차감 안 함 (XML 값이 이미 device IL로 가정).

정의:
    IL = max over (λ in λ_c ± 5nm) of IL_mzm(λ, V=-1V)
       = MZM이 ON 상태일 때의 raw transmission

값은 dB, 보통 음수.
표준 Si MZM IL 범위: 약 -1 ~ -10 dB (커플러 포함 시 더 큼).
"""
import numpy as np

BIAS_FIXED = -1.0


def extract_il(die, half_window_nm=5.0, bias=BIAS_FIXED):
    lam_c = die['lam_c']

    if bias not in die['sweeps']:
        biases = list(die['sweeps'].keys())
        bias = min(biases, key=lambda v: abs(v - BIAS_FIXED))

    L, IL_mzm = die['sweeps'][bias]
    mask = (L >= lam_c - half_window_nm) & (L <= lam_c + half_window_nm)
    if mask.sum() == 0:
        return float('nan')
    return round(float(IL_mzm[mask].max()), 3)
