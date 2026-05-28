"""O-band 전용 다이별 6 그림 생성.

O-band 디바이스 디자인 (MZMOTE_LULAB_380_500):
    - width  380 nm  (C-band 450 nm 보다 좁음)
    - length 500 um
    - 중심 파장 1310 nm, FSR ≈ 9.87 nm  (C-band 14.3 nm 보다 짧음)
    - ER 윈도우 1306–1320 nm (14 nm 폭)

O-band 가 C-band 와 다르게 fit 이 잘 안 맞는 경향이 있어 별도 파라미터 적용:
    - envelope_poly_deg 4 (C 의 3 보다 1차 높음) — O-band grating coupler
      응답이 더 좁고 가팔라 부드러운 envelope 을 잡으려면 한 차수 더 필요
    - trim_frac 0.08 (C 의 0.05 보다 큼) — O-band 양 끝단의 outlier spike 가
      더 두드러져 더 많이 잘라내야 안정
    - fsr_expected 9.87 (조교 피드백; 우리 실측 평균은 9.81)

실행:
    python3 src/plot_per_die_figures_O.py
"""
import sys, os, multiprocessing
PROGRAM_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROGRAM_ROOT, 'src'))

from plot_per_die_common import (
    DATA_ROOT, OUT_ROOT, find_all_xmls_with_meta, process_one,
)


# ──────────────────────────────────────────────────────────────────────
# O-band 파라미터
# ──────────────────────────────────────────────────────────────────────
O_CONFIG = {
    'band':                  'O',
    'er_window':             (1306.0, 1320.0),
    'envelope_poly_deg':     4,         # O-band 의 가파른 grating coupler 응답 잡으려면 한 차수 더
    'envelope_peak_window':  800,
    'fsr_expected':          9.87,
    'trim_frac':             0.08,      # 양 끝단 spike 더 두드러짐 → 더 많이 자름
    'ref_poly_deg':          6,
    'mzi_target_bias':       -1.0,
}


def _worker(args):
    return process_one(args, O_CONFIG)


def main():
    print('=' * 70)
    print(' O-band 다이별 6 그림 생성')
    print('=' * 70)
    print(f'config: {O_CONFIG}')
    items = [it for it in find_all_xmls_with_meta(DATA_ROOT) if it[5] == 'O']
    print(f'\n발견된 O-band 측정: {len(items)}개')
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
