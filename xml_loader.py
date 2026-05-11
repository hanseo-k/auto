"""XML 파싱 공통 모듈.

각 다이의 XML 파일에서:
  - 식별자 (Wafer, Row, Col, Width)
  - ALIGN 레퍼런스 스펙트럼 (커플러+도파로 손실)
  - MZM 스펙트럼 6개 (바이어스별)
  - IV 측정값
을 한 dict로 묶어 반환.

Band 자동 판정: LMZC = C-band(1550), LMZO = O-band(1310).
"""
import xml.etree.ElementTree as ET
import numpy as np
import glob, os, re


BAND_OF_FILE = {'LMZC': 'C', 'LMZO': 'O'}
BAND_CENTER = {'C': 1550.0, 'O': 1310.0}


def _arr(node):
    return np.array([float(x) for x in node.text.split(',')])


def _parse_width(modulator_name):
    """'MZMCTE_LULAB_450_500' → 450"""
    m = re.search(r'_LULAB_(\d+)_', modulator_name)
    return int(m.group(1)) if m else None


def load_die(xml_path):
    """단일 XML → 다이 데이터 dict.

    반환 dict 구조:
        wafer, band, row, col, width_nm, lam_c (밴드 중심파장)
        ref_L, ref_IL  (ALIGN 스펙트럼)
        sweeps : {V_bias: (L, IL_mzm)}  # 6개 바이어스
        iv_V, iv_I  (IV 측정)
    None을 반환하면 파싱 실패 (skip).
    """
    fname = os.path.basename(xml_path)
    band_tag = next((tag for tag in BAND_OF_FILE if f'_DCM_{tag}' in fname), None)
    if band_tag is None:
        return None
    band = BAND_OF_FILE[band_tag]
    lam_c = BAND_CENTER[band]

    try:
        root = ET.parse(xml_path).getroot()
    except ET.ParseError:
        return None
    tsi = root.find('TestSiteInfo')
    if tsi is None:
        return None

    out = {
        'wafer': tsi.attrib.get('Wafer', ''),
        'band': band,
        'row': int(tsi.attrib['DieRow']),
        'col': int(tsi.attrib['DieColumn']),
        'lam_c': lam_c,
        'width_nm': None,
        'ref_L': None, 'ref_IL': None,
        'sweeps': {}, 'iv_V': None, 'iv_I': None,
    }

    # ALIGN reference
    for m in root.iter('Modulator'):
        if 'ALIGN' not in m.attrib.get('Name', '').upper():
            continue
        ws = m.find('PortCombo/WavelengthSweep')
        if ws is None:
            continue
        out['ref_L'] = _arr(ws.find('L'))
        out['ref_IL'] = _arr(ws.find('IL'))
        break

    # MZM modulator (non-ALIGN)
    for mod in root.iter('Modulator'):
        name = mod.attrib.get('Name', '')
        if 'ALIGN' in name.upper():
            continue
        out['width_nm'] = _parse_width(name)
        for pc in mod.findall('PortCombo'):
            iv = pc.find('IVMeasurement')
            if iv is not None:
                out['iv_V'] = _arr(iv.find('Voltage'))
                out['iv_I'] = _arr(iv.find('Current'))
            for ws in pc.findall('WavelengthSweep'):
                if 'DCBias' not in ws.attrib:
                    continue
                v = round(float(ws.attrib['DCBias']), 2)
                out['sweeps'][v] = (_arr(ws.find('L')), _arr(ws.find('IL')))
            if out['sweeps']:
                break
        if out['sweeps']:
            break
    return out if out['sweeps'] else None


def find_all_xmls(root_dir):
    """HY202103 폴더 아래 모든 LMZC, LMZO XML 찾고 (Wafer, Row, Col, Band)로 dedupe.
    같은 다이가 같은 밴드로 여러 번 측정됐으면 마지막(최신 날짜) 채택.
    한 다이가 C와 O 둘 다 측정됐으면 두 항목 모두 유지."""
    out = {}
    for tag, band in BAND_OF_FILE.items():
        files = sorted(glob.glob(os.path.join(root_dir, '**', f'*_DCM_{tag}.xml'),
                                  recursive=True))
        for fp in files:
            m = re.search(r'(D\d+)_\((-?\d+),(-?\d+)\)', os.path.basename(fp))
            if not m:
                continue
            key = (m.group(1), int(m.group(2)), int(m.group(3)), band)
            out[key] = fp  # 같은 (다이, 밴드) 내에서만 최신으로 덮어씀
    return list(out.values())


def t_dev(L_mzm, IL_mzm, ref_L, ref_IL):
    """T_dev(λ) = IL_mzm(λ) − IL_ref(λ).
    ALIGN을 MZM 파장 그리드에 보간하여 차감."""
    if ref_L is None or ref_IL is None:
        return IL_mzm.copy()  # 차감 없음 (fallback)
    ref_interp = np.interp(L_mzm, ref_L, ref_IL)
    return IL_mzm - ref_interp
