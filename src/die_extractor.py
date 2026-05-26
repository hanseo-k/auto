"""다이별 파라미터 추출 — 공통 함수 `process_die`.

run.py (dedup main 분석) 와 analyze_by_date.py (날짜별 보존) 가
같은 metric 셋을 추출하도록 통일.  분석 도메인이 달라져도 컬럼
구성은 동일하게 유지된다.

반환 dict 키 (Date 컬럼은 caller 가 추가):
    Wafer, Band, Row, Col, Width_nm, Length_um,
    ER_dB, IL_dB, Vpi_V, Vpi_L_V_cm,
    FSR_nm, dlam_dV_pm_per_V, R2_dlam_vs_V, quality_grade, vpi_status,
    amplitude_ratio_k, power_split_ratio, imbalance_dB, mzm_loss_dB,
    imbalance_V<bias>_dB ...
"""
import pandas as pd

from xml_loader import load_die
from extract_er import extract_er
from extract_il import extract_il
from extract_vpi import extract_vpi, linearity_grade
from extract_passive_params import extract_passive_params


def process_die(xml_path):
    """단일 XML 경로 → 한 행 dict.  파싱 실패 시 None."""
    die = load_die(xml_path)
    if die is None:
        return None
    er       = extract_er(die)
    il       = extract_il(die)
    vpi_info = extract_vpi(die)
    passive  = extract_passive_params(die, er)

    vpi_V = vpi_info['vpi_V']
    L_um  = die.get('length_um')
    if (vpi_V is not None and not pd.isna(vpi_V)
            and L_um is not None):
        vpi_L_V_cm = round(vpi_V * L_um * 1e-4, 4)
    else:
        vpi_L_V_cm = float('nan')

    row = {
        'Wafer':    die['wafer'],
        'Band':     die['band'],
        'Row':      die['row'],
        'Col':      die['col'],
        'Width_nm':  die['width_nm'],
        'Length_um': die['length_um'],
        'ER_dB':    er,
        'IL_dB':    il,
        'Vpi_V':    vpi_V,
        'Vpi_L_V_cm':       vpi_L_V_cm,
        'FSR_nm':           vpi_info['fsr_nm'],
        'dlam_dV_pm_per_V': vpi_info['dlam_dV_pm_per_V'],
        'R2_dlam_vs_V':     vpi_info['linearity_R2'],
        'quality_grade':    linearity_grade(vpi_info['linearity_R2'],
                                            vpi_info['vpi_status']),
        'vpi_status':       vpi_info['vpi_status'],
        'amplitude_ratio_k': passive['amplitude_ratio_k'],
        'power_split_ratio': passive['power_split_ratio'],
        'imbalance_dB':      passive['imbalance_dB'],
        'mzm_loss_dB':       passive['mzm_loss_dB'],
    }
    for V in sorted(passive['imbalance_per_bias_dB'].keys()):
        row[f'imbalance_V{V:+.2f}_dB'] = passive['imbalance_per_bias_dB'][V]
    return row
