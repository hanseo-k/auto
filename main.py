"""
========================================================================
 메인 분석 프로그램
========================================================================

데이터 흐름:
    1) HY202103 폴더에서 모든 LMZC/LMZO XML 자동 탐색
    2) 다이별로 ER, IL, V_π 추출 (각각 별도 모듈)
    3) Outlier 검출 (물리 한계 + Hampel filter)
    4) 실행 시각 폴더 만들고 CSV 저장
    5) 웨이퍼맵 (연속 surface) 생성 — ER/IL/Vpi
    6) 1D 분포 그래프 생성 — outlier는 속 빈 마커

VSCode ▶ (F5) 로 그대로 실행 가능. 절대경로 사용.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd

from xml_loader import find_all_xmls, load_die
from extract_er import extract_er
from extract_il import extract_il
from extract_vpi import extract_vpi
from outlier_detect import mark_outliers
from csv_export import make_run_dir, export_csv
import wafer_map
import plot_1d
import plot_1d_mad
import trust_map


DATA_ROOT = '/Users/gimhanseo/Desktop/공프/HY202103'


def process_die(xml_path):
    """단일 XML → 한 행 dict (ER, IL, Vpi 포함)."""
    die = load_die(xml_path)
    if die is None:
        return None
    er = extract_er(die)
    il = extract_il(die)
    vpi_info = extract_vpi(die)
    return {
        'Wafer':    die['wafer'],
        'Band':     die['band'],
        'Row':      die['row'],
        'Col':      die['col'],
        'Width_nm': die['width_nm'],
        'ER_dB':    er,
        'IL_dB':    il,
        'Vpi_V':    vpi_info['vpi_V'],
        'FSR_nm':           vpi_info['fsr_nm'],
        'dlam_dV_pm_per_V': vpi_info['dlam_dV_pm_per_V'],
    }


def main():
    print('=' * 60)
    print(' MZM 4-wafer 분석 시작')
    print('=' * 60)

    # 1) XML 수집
    xmls = find_all_xmls(DATA_ROOT)
    print(f'\n[1/6] 발견된 다이 XML: {len(xmls)}개')

    # 2) 다이별 추출
    print('[2/6] ER, IL, V_π 추출 중...')
    rows = []
    for fp in xmls:
        row = process_die(fp)
        if row is not None:
            rows.append(row)
    df = pd.DataFrame(rows).sort_values(['Band', 'Wafer', 'Row', 'Col'])
    df = df.reset_index(drop=True)
    print(f'        → {len(df)}개 다이 처리 완료')

    # 3) Outlier 검출
    print('[3/6] Outlier 검출 (물리 바운드 + Hampel)...')
    df = mark_outliers(df)
    n_trusted = int(df['is_trusted'].sum())
    print(f'        → 신뢰 다이 {n_trusted}/{len(df)}개')
    for col in ['ER_dB', 'IL_dB', 'Vpi_V']:
        n_out = int(df[f'is_outlier_{col}'].sum())
        print(f'         {col}: outlier {n_out}개')

    # 4) 실행 폴더 + CSV
    print('[4/6] 결과 폴더 생성 + CSV 저장...')
    run_dir = make_run_dir()
    print(f'        → {run_dir}')
    export_csv(df, run_dir)

    # 5) 웨이퍼맵
    print('[5/7] 웨이퍼맵 생성 (연속 surface)...')
    wafer_map.plot_all(df, run_dir)

    # 6) 1D 분포 그래프
    print('[6/7] 1D 분포 그래프 생성...')
    plot_1d.plot_all(df, run_dir)
    plot_1d_mad.plot_all(df, run_dir)

    # 7) 신뢰도 맵
    print('[7/7] 신뢰도 맵 생성 (Hampel 이웃 수)...')
    trust_map.plot_all(df, run_dir)

    print('\n' + '=' * 60)
    print(' 완료!  결과 위치:')
    print(f'   {run_dir}')
    print('=' * 60)

    # 요약 통계 (신뢰 데이터만)
    trusted = df[df['is_trusted']]
    print('\n[중앙값 요약 (신뢰 데이터만)]')
    summary = trusted.groupby(['Band', 'Wafer']).agg(
        n=('Row', 'count'),
        ER=('ER_dB', 'median'),
        IL=('IL_dB', 'median'),
        Vpi=('Vpi_V', 'median'),
    ).round(2)
    print(summary)


if __name__ == '__main__':
    main()
