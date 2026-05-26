"""C-band 전용 다이별 6 그림 생성.

C-band 디바이스 디자인 (MZMCTE_LULAB_450_500):
    - width  450 nm
    - length 500 um
    - 중심 파장 1550 nm, FSR ≈ 14.3 nm
    - ER 윈도우 1545–1561 nm (16 nm 폭)

실행:
    python3 src/plot_per_die_figures_C.py
"""
import sys, os, multiprocessing
PROGRAM_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROGRAM_ROOT, 'src'))

from plot_per_die_common import (
    DATA_ROOT, OUT_ROOT, find_all_xmls_with_meta, process_one,
)
from functools import partial


# ──────────────────────────────────────────────────────────────────────
# C-band 파라미터
# ──────────────────────────────────────────────────────────────────────
C_CONFIG = {
    'band':                  'C',
    'er_window':             (1545.0, 1561.0),
    'envelope_poly_deg':     3,
    'envelope_peak_window':  800,
    'fsr_expected':          14.3,
    'trim_frac':             0.05,
    'ref_poly_deg':          6,
    'mzi_target_bias':       -1.0,
}


def _worker(args):
    return process_one(args, C_CONFIG)


def main():
    print('=' * 70)
    print(' C-band 다이별 6 그림 생성')
    print('=' * 70)
    print(f'config: {C_CONFIG}')
    items = [it for it in find_all_xmls_with_meta(DATA_ROOT) if it[5] == 'C']
    print(f'\n발견된 C-band 측정: {len(items)}개')
    print(f'출력 위치: {OUT_ROOT}\n')

    with multiprocessing.Pool() as pool:
        results = pool.map(_worker, items)
    ok = sum(1 for r in results if r.startswith('OK'))
    fail = len(results) - ok
    print(f'\n완료: {ok}개 OK, {fail}개 실패')
    if fail > 0:
        for r in results:
            if not r.startswith('OK'):
                print('  ' + r)


if __name__ == '__main__':
    main()
