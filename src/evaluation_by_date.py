"""날짜별 다이별 종합 평가 — evaluation.py 의 by-date 버전.

evaluation.py 가 dedup 된 data.csv (다이 1행) 를 평가하는 반면,
이 모듈은 data_by_date.csv (모든 측정 날짜 보존) 를 평가한다.

차이점:
    - Date 컬럼 추가
    - Robust Z 그룹화 단위 = (Date, Wafer, Band)
      (그룹 안 표본이 14개 그대로 유지되어 의미 있음)
    - 그 외 평가 카테고리/임계값 등은 evaluation.py 와 동일

출력:
    res/csv/evaluation_by_date.csv
    res/csv/evaluation_by_date.xlsx
"""
import sys, os
PROGRAM_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROGRAM_ROOT, 'src'))

import pandas as pd
import numpy as np

from outlier_detect import PHYSICAL_BOUNDS
from evaluation import (
    _evaluate_row, _find_figure_dir, write_xlsx as _write_xlsx_base,
)


CSV_PATH = os.path.join(PROGRAM_ROOT, 'res', 'csv', 'data_by_date.csv')
OUT_CSV  = os.path.join(PROGRAM_ROOT, 'res', 'csv', 'evaluation_by_date.csv')
OUT_XLSX = os.path.join(PROGRAM_ROOT, 'res', 'csv', 'evaluation_by_date.xlsx')


# ──────────────────────────────────────────────────────────────────────
# Group 단위 Robust Z (Date, Wafer, Band)
# ──────────────────────────────────────────────────────────────────────
def add_robust_z_per_group(df):
    """data_by_date.csv 에는 robust_z_* 가 없으므로 여기서 계산.

    그룹: (Date, Wafer, Band).  각 그룹 안에서 modified Z-score.
    """
    for col in ['ER_dB', 'IL_dB', 'Vpi_V']:
        z_col = f'robust_z_{col}'
        df[z_col] = np.nan
        for (date, wafer, band), grp in df.groupby(['Date', 'Wafer', 'Band']):
            vals = grp[col].dropna()
            if len(vals) < 3:
                continue
            m = vals.median()
            mad = (vals - m).abs().median()
            sigma = 1.4826 * mad
            if sigma == 0:
                continue
            for idx in grp.index:
                x = grp.loc[idx, col]
                if pd.isna(x):
                    continue
                df.at[idx, z_col] = round((x - m) / sigma, 2)
    return df


def _find_figure_dir_by_date(date, wafer, band, row, col):
    """date 가 주어진 경우 정확한 폴더로."""
    path = os.path.join(PROGRAM_ROOT, 'res', 'figures_per_die', date,
                        f'{band}-band', wafer, f'({row},{col})')
    return path if os.path.isdir(path) else ''


def build_evaluation_by_date(df):
    """data_by_date.csv DataFrame → 평가 결과 DataFrame."""
    df = add_robust_z_per_group(df.copy())

    base_cols = ['Date', 'Wafer', 'Band', 'Row', 'Col', 'Width_nm', 'Length_um',
                 'ER_dB', 'IL_dB', 'Vpi_V', 'Vpi_L_V_cm',
                 'R2_dlam_vs_V', 'quality_grade',
                 'imbalance_dB', 'mzm_loss_dB']
    base_cols = [c for c in base_cols if c in df.columns]

    rows = []
    for _, r in df.iterrows():
        ev = _evaluate_row(r)
        rec = {c: r.get(c) for c in base_cols}
        rec.update(ev)
        rec['figures'] = _find_figure_dir_by_date(
            r['Date'], r['Wafer'], r['Band'], int(r['Row']), int(r['Col']))
        rows.append(rec)

    out_df = pd.DataFrame(rows)
    col_order = (
        base_cols +
        ['overall', 'vpi_extraction', 'er_phys', 'il_phys', 'vpi_phys',
         'er_z', 'il_z', 'vpi_z', 'linearity', 'imbalance', 'loss',
         'notes', 'figures']
    )
    col_order = [c for c in col_order if c in out_df.columns]
    return out_df[col_order]


def main():
    df = pd.read_csv(CSV_PATH)
    print('=' * 70)
    print(f' 날짜별 평가 (evaluation_by_date) — 측정 {len(df)} 개')
    print('=' * 70)

    eval_df = build_evaluation_by_date(df)
    eval_df.to_csv(OUT_CSV, index=False, encoding='utf-8-sig')
    _write_xlsx_base(eval_df, OUT_XLSX)

    counts = eval_df['overall'].value_counts()
    print('\n[종합 분포]')
    for k in ('PASS', 'WARN', 'FAIL'):
        print(f'  {k:5s}: {counts.get(k, 0):3d} 측정')

    print('\n[날짜별 종합]')
    for date, grp in eval_df.groupby('Date'):
        n = len(grp)
        f = (grp['overall'] == 'FAIL').sum()
        w = (grp['overall'] == 'WARN').sum()
        p = (grp['overall'] == 'PASS').sum()
        print(f'  {date}:  PASS {p}/{n}, WARN {w}/{n}, FAIL {f}/{n}')

    fails = eval_df[eval_df['overall'] == 'FAIL']
    if not fails.empty:
        print(f'\n[FAIL {len(fails)}개 — 앞 10개]')
        for _, r in fails.head(10).iterrows():
            print(f'  {r["Date"]} {r["Wafer"]} {r["Band"]} '
                  f'({int(r["Row"])},{int(r["Col"])}): {r["notes"]}')

    print(f'\nCSV  저장: {OUT_CSV}')
    print(f'XLSX 저장: {OUT_XLSX}')


if __name__ == '__main__':
    main()
