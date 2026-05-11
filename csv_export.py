"""실행 시각으로 폴더 만들고 CSV 저장.

폴더명: results/YYYY-MM-DD_HH-MM-SS/
CSV 컬럼: Wafer, Band, Row, Col, Width_nm, ER_dB, IL_dB, Vpi_V,
         FSR_nm, dlam_dV_pm_per_V, is_trusted, is_outlier_*
"""
import os
import datetime
import pandas as pd


PROGRAM_ROOT = '/Users/gimhanseo/Desktop/공프/이거 개쩌는 프로그램'


def make_run_dir():
    """실행 시각 기반 결과 폴더 생성, 경로 반환."""
    ts = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    out = os.path.join(PROGRAM_ROOT, 'results', ts)
    os.makedirs(out, exist_ok=True)
    return out


def export_csv(df, run_dir, filename='data.csv'):
    """DataFrame을 run_dir 안에 CSV로 저장."""
    path = os.path.join(run_dir, filename)
    df.to_csv(path, index=False, encoding='utf-8-sig')
    print(f'CSV saved: {path}')
    return path
