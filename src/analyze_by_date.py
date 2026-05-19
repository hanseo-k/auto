"""날짜별 분석 — HY202103 폴더의 모든 측정 날짜를 분리해서 CSV 생성.

`xml_loader.find_all_xmls` 는 같은 (다이, 밴드) 에 대해 최신 측정만 남기지만,
이 스크립트는 모든 날짜의 측정을 보존해서 시간 변화를 추적할 수 있게 함.

출력:
    res/csv/data_by_date.csv
    컬럼: Date, Wafer, Band, Row, Col, Width_nm, ER_dB, IL_dB, Vpi_V

실행:
    python3 src/analyze_by_date.py
"""
import sys, os, re, glob
PROGRAM_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROGRAM_ROOT, 'src'))

import pandas as pd
import multiprocessing

from xml_loader import load_die, BAND_OF_FILE
from extract_er import extract_er
from extract_il import extract_il
from extract_vpi import extract_vpi


DATA_ROOT = '/Users/gimhanseo/Desktop/공프/HY202103'


def find_all_xmls_with_dates(root_dir):
    """모든 LMZC/LMZO XML 반환 (dedup 없음).

    경로에 포함된 YYYYMMDD_HHMMSS 폴더에서 날짜를 추출.
    반환: [(xml_path, 'YYYY-MM-DD'), ...]
    """
    items = []
    for tag in BAND_OF_FILE:
        pattern = os.path.join(root_dir, '**', f'*_DCM_{tag}.xml')
        for fp in sorted(glob.glob(pattern, recursive=True)):
            m = re.search(r'(\d{4})(\d{2})(\d{2})_\d{6}', fp)
            if not m:
                continue
            date_str = f'{m.group(1)}-{m.group(2)}-{m.group(3)}'
            items.append((fp, date_str))
    return items


def _process(args):
    xml_path, date_str = args
    die = load_die(xml_path)
    if die is None:
        return None
    return {
        'Date':     date_str,
        'Wafer':    die['wafer'],
        'Band':     die['band'],
        'Row':      die['row'],
        'Col':      die['col'],
        'Width_nm': die['width_nm'],
        'ER_dB':    extract_er(die),
        'IL_dB':    extract_il(die),
        'Vpi_V':    extract_vpi(die)['vpi_V'],
    }


def analyze(data_root=DATA_ROOT):
    """모든 날짜 측정을 추출해서 DataFrame 반환."""
    items = find_all_xmls_with_dates(data_root)
    with multiprocessing.Pool() as pool:
        results = pool.map(_process, items)
    rows = [r for r in results if r is not None]
    df = pd.DataFrame(rows).sort_values(['Date', 'Band', 'Wafer', 'Row', 'Col'])
    return df.reset_index(drop=True)


def main():
    print('=' * 60)
    print(' 날짜별 분석 시작')
    print('=' * 60)

    items = find_all_xmls_with_dates(DATA_ROOT)
    print(f'\n발견된 다이-측정: {len(items)}개')
    print('병렬 추출 중...')

    df = analyze(DATA_ROOT)
    print(f'        → {len(df)}개 측정 처리 완료')

    out_dir = os.path.join(PROGRAM_ROOT, 'res', 'csv')
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, 'data_by_date.csv')
    df.to_csv(out_path, index=False, encoding='utf-8-sig')
    print(f'\nCSV 저장: {out_path}')

    # 날짜별 요약
    print('\n[날짜별 측정 수 (Wafer/Band)]')
    summary = df.groupby(['Date', 'Wafer', 'Band']).size().unstack(
        level=['Wafer', 'Band'], fill_value=0)
    print(summary)

    print('\n[날짜별 중앙값]')
    medians = df.groupby('Date')[['ER_dB', 'IL_dB', 'Vpi_V']].median().round(2)
    print(medians)


if __name__ == '__main__':
    main()
