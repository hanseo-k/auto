"""
========================================================================
 메인 분석 프로그램 — MZM 4-wafer 분석
========================================================================

데이터 흐름:
    1) HY202103 폴더에서 모든 LMZC/LMZO XML 자동 탐색
    2) 다이별로 ER, IL, V_π 추출 — 멀티코어 병렬 처리
    3) Outlier 검출 (물리 한계 + Robust Z-score)
    4) 실행 시각 폴더 만들고 CSV 저장
    5) 웨이퍼맵 / 1D분포 / 1D+MAD / Z-score맵 — 12개 동시 생성
    6) 최신 결과를 res/csv, res/figures 로 복사 후 GitHub 자동 push

폴더 구조:
    src/         - 분석 모듈
    data/        - 입력 데이터 (HY202103 별도 위치 참조)
    res/csv/     - 최신 CSV (git 추적)
    res/figures/ - 최신 그림 (git 추적)
    res/<ts>/    - 실행 시각별 보관 (gitignored)
    doc/         - 방법론 그림 등
"""
import sys, os
PROGRAM_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(PROGRAM_ROOT, 'src'))

import pandas as pd
import multiprocessing
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor

from xml_loader import find_all_xmls, load_die
from extract_er import extract_er
from extract_il import extract_il
from extract_vpi import extract_vpi, linearity_grade
from extract_passive_params import extract_passive_params
from outlier_detect import mark_outliers
from csv_export import make_run_dir, export_csv
import wafer_map
import plot_1d
import plot_1d_mad
import zscore_map
import decompose_variation
import analyze_by_date
import evaluation


DATA_ROOT = os.path.join(PROGRAM_ROOT, 'data', 'HY202103')


def process_die(xml_path):
    """단일 XML → 한 행 dict (ER, IL, Vpi 포함)."""
    die = load_die(xml_path)
    if die is None:
        return None
    er = extract_er(die)
    il = extract_il(die)
    vpi_info = extract_vpi(die)
    passive  = extract_passive_params(die, er)
    # V_pi·L : 디바이스 표준 spec (V·cm)
    vpi_V = vpi_info['vpi_V']
    L_um = die.get('length_um')
    vpi_L_V_cm = (vpi_V * L_um * 1e-4) if (vpi_V is not None
                                            and not pd.isna(vpi_V)
                                            and L_um is not None) else float('nan')

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
        'Vpi_L_V_cm':       round(vpi_L_V_cm, 4) if not pd.isna(vpi_L_V_cm) else float('nan'),
        'FSR_nm':           vpi_info['fsr_nm'],
        'dlam_dV_pm_per_V': vpi_info['dlam_dV_pm_per_V'],
        'R2_dlam_vs_V':     vpi_info['linearity_R2'],
        'quality_grade':    linearity_grade(vpi_info['linearity_R2'],
                                            vpi_info['vpi_status']),
        'vpi_status':       vpi_info['vpi_status'],
        # ── Splitter 영역 ─────────────────────────────────────────
        'amplitude_ratio_k': passive['amplitude_ratio_k'],
        'power_split_ratio': passive['power_split_ratio'],
        'imbalance_dB':      passive['imbalance_dB'],
        'mzm_loss_dB':       passive['mzm_loss_dB'],
    }
    # 바이어스별 imbalance — 'splitter 영역' 의 상세 컬럼 (V 순서대로)
    for V in sorted(passive['imbalance_per_bias_dB'].keys()):
        row[f'imbalance_V{V:+.2f}_dB'] = passive['imbalance_per_bias_dB'][V]
    return row


METRICS = [
    ('ER_dB', 'Extinction Ratio (dB)'),
    ('IL_dB', 'Insertion Loss (dB)'),
    ('Vpi_V', 'V_pi (V)'),
]


# ── 플롯 작업 래퍼 — metric별로 쪼개서 12개 작업으로 (코어 풀 활용) ──
def _run_plot(args):
    """(plot_type, col, label, df, run_dir) → 개별 그래프 1장 생성."""
    import os
    plot_type, col, label, df, run_dir = args
    if plot_type == 'wafer':
        import wafer_map as _m
        _m.plot_wafer_map(df, col, label, os.path.join(run_dir, f'wafer_map_{col}.png'))
    elif plot_type == '1d':
        import plot_1d as _m
        _m.plot_1d(df, col, label, os.path.join(run_dir, f'1d_{col}.png'))
    elif plot_type == '1d_mad':
        import plot_1d_mad as _m
        _m.plot_1d_mad(df, col, label, os.path.join(run_dir, f'1d_mad_{col}.png'))
    elif plot_type == 'zscore':
        import zscore_map as _m
        _m.plot_zscore_map(df, col, label, os.path.join(run_dir, f'zscore_map_{col}.png'))
    elif plot_type == 'decompose':
        import decompose_variation as _m
        _m.plot_decomposition(df, col, label, os.path.join(run_dir, f'decompose_{col}.png'))


def _sync_to_res(run_dir):
    """run_dir의 최신 결과를 res/csv, res/figures 로 분리 복사."""
    csv_dir = os.path.join(PROGRAM_ROOT, 'res', 'csv')
    fig_dir = os.path.join(PROGRAM_ROOT, 'res', 'figures')
    os.makedirs(csv_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)
    # 기존 파일 제거 (잔여물 방지)
    for d in (csv_dir, fig_dir):
        for f in os.listdir(d):
            fp = os.path.join(d, f)
            if os.path.isfile(fp):
                os.remove(fp)
    # 새로 복사
    for f in os.listdir(run_dir):
        src = os.path.join(run_dir, f)
        if not os.path.isfile(src):
            continue
        if f.endswith('.csv'):
            shutil.copy2(src, os.path.join(csv_dir, f))
        elif f.endswith('.png'):
            shutil.copy2(src, os.path.join(fig_dir, f))


def main():
    print('=' * 60)
    print(' MZM 4-wafer 분석 시작')
    print('=' * 60)

    # 1) XML 수집
    xmls = find_all_xmls(DATA_ROOT)
    print(f'\n[1/6] 발견된 다이 XML: {len(xmls)}개')

    # 2) 다이별 추출 — 멀티코어 병렬
    print('[2/6] ER, IL, V_π 추출 중... (멀티코어)')
    with multiprocessing.Pool() as pool:
        results = pool.map(process_die, xmls)
    rows = [r for r in results if r is not None]
    df = pd.DataFrame(rows).sort_values(['Band', 'Wafer', 'Row', 'Col'])
    df = df.reset_index(drop=True)
    print(f'        → {len(df)}개 다이 처리 완료')

    # 3) Outlier 검출
    print('[3/6] Outlier 검출 (물리 바운드 + Robust Z)...')
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

    # 5) 플롯 12개 동시 생성
    print('[5/6] 플롯 생성 (4종 × 3 metric = 12작업, 멀티코어 병렬)...')
    tasks = [
        (kind, col, label, df, run_dir)
        for kind in ('wafer', '1d', '1d_mad', 'zscore', 'decompose')
        for col, label in METRICS
    ]
    with ProcessPoolExecutor(max_workers=None) as ex:
        futures = [ex.submit(_run_plot, t) for t in tasks]
        for f in futures:
            f.result()
    print('        → 모든 플롯 완료')

    # 6a) res/csv, res/figures 동기화 (먼저 — 기존 파일 지우고 새로 복사)
    print('[6/6] res/csv, res/figures 동기화...')
    _sync_to_res(run_dir)

    # 6b) 날짜별 분석 — sync 뒤에 실행해야 결과 보존됨
    print('       날짜별 분석 (data_by_date.csv + by_date_summary.png)...')
    analyze_by_date.export_and_plot()

    # 6c) 다이별 평가 (evaluation.csv + evaluation.xlsx)
    print('       다이별 평가 (evaluation.csv + evaluation.xlsx)...')
    evaluation.main()

    # 6c) GitHub push
    subprocess.run(['git', '-C', PROGRAM_ROOT, 'add', '-A'], check=True)
    status = subprocess.run(['git', '-C', PROGRAM_ROOT, 'diff', '--cached', '--quiet'])
    if status.returncode != 0:
        subprocess.run(['git', '-C', PROGRAM_ROOT, 'commit', '-m',
                        f'Update latest results ({os.path.basename(run_dir)})'], check=True)
        subprocess.run(['git', '-C', PROGRAM_ROOT, 'push', 'origin', 'main'], check=True)
        print('        → GitHub 업로드 완료')
    else:
        print('        → 변경사항 없음 — 스킵')

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
        VpiL=('Vpi_L_V_cm', 'median'),
    ).round(3)
    print(summary)

    # Width × Length 별 통계 (디바이스 디자인 그룹화)
    print('\n[Width × Length 그룹별 통계 (신뢰 데이터만)]')
    width_summary = trusted.groupby(['Width_nm', 'Length_um', 'Band']).agg(
        n=('Row', 'count'),
        ER_med=('ER_dB', 'median'), ER_std=('ER_dB', 'std'),
        IL_med=('IL_dB', 'median'), IL_std=('IL_dB', 'std'),
        Vpi_med=('Vpi_V', 'median'), Vpi_std=('Vpi_V', 'std'),
        VpiL_med=('Vpi_L_V_cm', 'median'),
    ).round(3)
    print(width_summary)


if __name__ == '__main__':
    main()
