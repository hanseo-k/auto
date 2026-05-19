"""제안된 두 수정사항의 sensitivity 테스트.

테스트 A: V_π 추출 — 점프 임계값 / 최소 slope 필터
테스트 B: ER 윈도우 폭 변경

각각 여러 파라미터로 다시 추출해보고 결과 분포를 비교.
원본 코드는 건드리지 않고, 파라미터화된 사본을 임시로 만들어서 실행.

실행:
    python3 src/sensitivity_test.py
"""
import sys, os
PROGRAM_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROGRAM_ROOT, 'src'))

import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from xml_loader import load_die, t_dev
from extract_vpi import _parabolic_null
from analyze_by_date import find_all_xmls_with_dates


DATA_ROOT = '/Users/gimhanseo/Desktop/공프/HY202103'
OUT = os.path.join(PROGRAM_ROOT, 'doc', 'investigation')
os.makedirs(OUT, exist_ok=True)


# ──────────────────────────────────────────────────────────────────────
# 파라미터화된 V_π 추출
# ──────────────────────────────────────────────────────────────────────
def extract_vpi_param(die, jump_mult=1.5, min_slope_pm_per_V=0.0):
    """V_π 추출 파라미터 버전.

    jump_mult:           total_shift > half * jump_mult 이면 reject
    min_slope_pm_per_V:  |dλ/dV| 이 이 값 미만이면 reject (NaN 반환)
    """
    sweeps = die['sweeps']
    biases = sorted(sweeps.keys())
    if len(biases) < 3:
        return float('nan'), float('nan'), float('nan')

    L0, I0 = sweeps[biases[0]]
    distance = max(1, int(1.0 / (L0[1] - L0[0])))
    peaks, _ = find_peaks(-I0, prominence=10, distance=distance)
    deep = sorted([float(L0[p]) for p in peaks if I0[p] < -25])
    if len(deep) < 2:
        return float('nan'), float('nan'), float('nan')
    fsr = float(np.median(np.diff(deep)))

    rev_biases = [v for v in biases if v <= 0.0]
    if len(rev_biases) < 3:
        return fsr, float('nan'), float('nan')

    half = min(0.4, fsr * 0.35)
    slopes = []
    for lam0 in deep:
        positions = [_parabolic_null(*sweeps[v], lam0, half) for v in rev_biases]
        positions = np.array(positions)
        valid = ~np.isnan(positions)
        if valid.sum() < 3:
            continue
        total_shift = abs(positions[valid][-1] - positions[valid][0])
        if total_shift > half * jump_mult:
            continue
        s, _ = np.polyfit(np.array(rev_biases)[valid], positions[valid], 1)
        slopes.append(s)

    if not slopes:
        return fsr, float('nan'), float('nan')
    slopes = np.array(slopes)
    if len(slopes) >= 3:
        med = np.median(slopes); mad = np.median(np.abs(slopes - med))
        slopes = slopes[np.abs(slopes - med) <= 3 * mad + 1e-6]
    dlam_dV = float(np.mean(slopes))
    dlam_pm = dlam_dV * 1000

    # 최소 slope 필터
    if abs(dlam_pm) < min_slope_pm_per_V:
        return fsr, dlam_pm, float('nan')
    if dlam_dV == 0:
        return fsr, dlam_pm, float('nan')
    return fsr, dlam_pm, fsr / (2 * abs(dlam_dV))


# ──────────────────────────────────────────────────────────────────────
# 파라미터화된 ER 추출
# ──────────────────────────────────────────────────────────────────────
def extract_er_param(die, win_lo, win_hi):
    """주어진 윈도우 [win_lo, win_hi] nm 안에서 모든 바이어스의 peak−null."""
    ref_L, ref_IL = die['ref_L'], die['ref_IL']
    peak_t = -np.inf; null_t = np.inf
    for V, (L, IL_mzm) in die['sweeps'].items():
        T = t_dev(L, IL_mzm, ref_L, ref_IL)
        mask = (L >= win_lo) & (L <= win_hi)
        if mask.sum() == 0:
            continue
        peak_t = max(peak_t, float(T[mask].max()))
        null_t = min(null_t, float(T[mask].min()))
    if peak_t == -np.inf:
        return float('nan')
    return peak_t - null_t


# ──────────────────────────────────────────────────────────────────────
# 데이터 로드 (한 번만)
# ──────────────────────────────────────────────────────────────────────
def load_all_dies():
    items = find_all_xmls_with_dates(DATA_ROOT)
    dies = []
    for fp, date in items:
        d = load_die(fp)
        if d is None: continue
        d['_date'] = date
        d['_fp'] = fp
        dies.append(d)
    return dies


# ──────────────────────────────────────────────────────────────────────
# 테스트 A — V_π 추출 파라미터 스윕
# ──────────────────────────────────────────────────────────────────────
def test_vpi(dies):
    print('\n' + '=' * 70)
    print(' 테스트 A — V_π 점프 임계값 / 최소 slope 필터 sensitivity')
    print('=' * 70)

    jump_mults = [0.5, 1.0, 1.5, 2.0, 3.0]
    min_slopes = [0, 10, 30, 50, 100]   # pm/V

    # 1) 점프 임계값 단독 스윕 (min_slope=0)
    print('\n[A-1] 점프 임계값만 변경 (min_slope=0)')
    print(f'{"jump_mult":>10} {"n_total":>8} {"n_nan":>7} {"n_over_60V":>11} {"median_vpi":>11} {"max_vpi":>10}')
    for jm in jump_mults:
        vpis = [extract_vpi_param(d, jump_mult=jm)[2] for d in dies]
        vpis = np.array(vpis)
        n_nan = int(np.isnan(vpis).sum())
        ok = vpis[~np.isnan(vpis)]
        n_over = int((ok > 60).sum())
        med = float(np.median(ok[ok <= 60])) if (ok <= 60).any() else float('nan')
        mx = float(ok.max()) if len(ok) else float('nan')
        print(f'{jm:>10.1f} {len(vpis):>8} {n_nan:>7} {n_over:>11} {med:>11.2f} {mx:>10.2f}')

    # 2) 최소 slope 필터 추가 (jump_mult=1.5 고정)
    print('\n[A-2] 최소 |dλ/dV| 필터 추가 (jump_mult=1.5)')
    print(f'{"min_slope":>10} {"n_nan":>7} {"n_over_60V":>11} {"median_vpi":>11} {"max_vpi":>10}')
    for ms in min_slopes:
        results = [extract_vpi_param(d, min_slope_pm_per_V=ms) for d in dies]
        vpis = np.array([r[2] for r in results])
        n_nan = int(np.isnan(vpis).sum())
        ok = vpis[~np.isnan(vpis)]
        n_over = int((ok > 60).sum())
        med = float(np.median(ok[ok <= 60])) if (ok <= 60).any() else float('nan')
        mx = float(ok.max()) if len(ok) else float('nan')
        print(f'{ms:>10.0f} {n_nan:>7} {n_over:>11} {med:>11.2f} {mx:>10.2f}')

    # 3) 정상 vs 망가진 데이터 별도 분석
    print('\n[A-3] 정상 측정(06-03 D23/D24) vs 망가진(05-31 D23/D24) 의 dλ/dV 비교')
    for date_label, date_filter in [
        ('  정상 (06-03)', '2019-06-03'),
        ('망가짐 (05-31)', '2019-05-31'),
    ]:
        sub = [d for d in dies if d['_date'] == date_filter and d['wafer'] in ('D23', 'D24')]
        dlams = []
        for d in sub:
            _, dl, _ = extract_vpi_param(d)
            if not np.isnan(dl):
                dlams.append(dl)
        if dlams:
            arr = np.array(dlams)
            print(f'{date_label}: n={len(arr)}, dλ/dV (pm/V) → '
                  f'median {np.median(arr):+.1f}, '
                  f'min {arr.min():+.1f}, max {arr.max():+.1f}, '
                  f'|min| {np.abs(arr).min():.2f}')

    # 4) 추천 셋팅으로 망가진 데이터 처리 결과
    print('\n[A-4] 추천: jump_mult=1.5, min_slope=30 pm/V — 효과 확인')
    bad = [d for d in dies if d['_date'] == '2019-05-31']
    saved = sum(1 for d in bad
                if np.isnan(extract_vpi_param(d, min_slope_pm_per_V=30)[2]))
    print(f'  05-31 망가진 28개 중 NaN 처리됨: {saved}/{len(bad)}')

    # 정상 데이터(06-03) 에서 NaN 으로 잘못 처리되는 다이 수
    clean = [d for d in dies if d['_date'] in ('2019-06-03', '2019-07-12', '2019-07-15', '2019-05-26')]
    false_pos = sum(1 for d in clean
                    if np.isnan(extract_vpi_param(d, min_slope_pm_per_V=30)[2]))
    print(f'  정상 측정({len(clean)}개) 중 NaN 잘못 처리됨: {false_pos}/{len(clean)}')

    # 그림: 추천 셋팅 적용 전후
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), dpi=120)
    for ax, jm, ms, title in [
        (axes[0], 1.5, 0,  'Before  (current: jump=1.5, min_slope=0)'),
        (axes[1], 1.5, 30, 'After   (proposed: jump=1.5, min_slope=30 pm/V)'),
    ]:
        vpis = [extract_vpi_param(d, jump_mult=jm, min_slope_pm_per_V=ms)[2]
                for d in dies]
        vpis = np.array(vpis)
        ok = vpis[~np.isnan(vpis)]
        ax.hist(np.clip(ok, 0, 200), bins=np.linspace(0, 200, 41),
                color='steelblue', alpha=0.7, edgecolor='black')
        ax.axvspan(60, 200, color='red', alpha=0.15, label='> 60 V (physical)')
        ax.set_xlabel('V_π (V)')
        ax.set_ylabel('count')
        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.set_xlim(0, 200)
        ax.legend()
        n_over = int((ok > 60).sum())
        n_nan = int(np.isnan(vpis).sum())
        ax.text(0.98, 0.95, f'NaN: {n_nan}\nover 60V: {n_over}',
                transform=ax.transAxes, ha='right', va='top',
                bbox=dict(facecolor='white', alpha=0.85), fontsize=10)
    plt.tight_layout()
    out = os.path.join(OUT, 'sens_A_vpi.png')
    plt.savefig(out, bbox_inches='tight'); plt.close(fig)
    print(f'  그림 저장: {out}')


# ──────────────────────────────────────────────────────────────────────
# 테스트 B — ER 윈도우 폭 sensitivity
# ──────────────────────────────────────────────────────────────────────
def test_er(dies):
    print('\n' + '=' * 70)
    print(' 테스트 B — ER 윈도우 폭 sensitivity')
    print('=' * 70)

    # C-band 윈도우 변형들 (중심 1553)
    c_windows = [
        ('현재   1546-1560 (14nm)', 1546.0, 1560.0),
        ('1nm   1545-1561 (16nm)', 1545.0, 1561.0),
        ('2nm   1544-1562 (18nm)', 1544.0, 1562.0),
        ('4nm   1542-1564 (22nm)', 1542.0, 1564.0),
        ('극단  1535-1571 (36nm)', 1535.0, 1571.0),
    ]
    o_windows = [
        ('현재   1306-1320 (14nm)', 1306.0, 1320.0),
        ('1nm   1305-1321 (16nm)', 1305.0, 1321.0),
        ('2nm   1304-1322 (18nm)', 1304.0, 1322.0),
        ('4nm   1302-1324 (22nm)', 1302.0, 1324.0),
    ]

    for band, windows in [('C', c_windows), ('O', o_windows)]:
        sub = [d for d in dies if d['band'] == band]
        # 같은 다이가 여러 날 있으면 dedup (최신만)
        seen = {}
        for d in sub:
            key = (d['wafer'], d['row'], d['col'])
            seen[key] = d
        sub = list(seen.values())
        print(f'\n[B-{band}]  {band}-band, n={len(sub)} 다이')
        print(f'{"label":<28} {"median":>8} {"std":>7} {"min":>7} {"max":>7} {"n>45":>5} {"n<10":>5}')
        for label, lo, hi in windows:
            ers = np.array([extract_er_param(d, lo, hi) for d in sub])
            ers = ers[~np.isnan(ers)]
            n_over = int((ers > 45).sum())
            n_under = int((ers < 10).sum())
            print(f'{label:<28} {np.median(ers):>8.2f} {np.std(ers):>7.2f} '
                  f'{ers.min():>7.2f} {ers.max():>7.2f} {n_over:>5d} {n_under:>5d}')

    # 그림: C-band ER 분포 변화
    fig, axes = plt.subplots(1, len(c_windows), figsize=(3.5 * len(c_windows), 4.5), dpi=120)
    sub = [d for d in dies if d['band'] == 'C']
    seen = {}
    for d in sub:
        key = (d['wafer'], d['row'], d['col']); seen[key] = d
    sub = list(seen.values())
    for ax, (label, lo, hi) in zip(axes, c_windows):
        ers = np.array([extract_er_param(d, lo, hi) for d in sub])
        ers = ers[~np.isnan(ers)]
        ax.hist(ers, bins=np.linspace(25, 50, 26), color='steelblue',
                alpha=0.7, edgecolor='black')
        ax.axvspan(45, 50, color='red', alpha=0.15, label='> 45 (above bound)')
        ax.axvline(np.median(ers), color='black', lw=2,
                   label=f'median {np.median(ers):.1f}')
        ax.set_xlim(25, 50)
        ax.set_title(label.split('(')[0].strip(), fontsize=10, fontweight='bold')
        ax.set_xlabel('ER (dB)'); ax.set_ylabel('count')
        ax.legend(fontsize=8)
    fig.suptitle('Test B  —  C-band ER distribution vs window width',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    out = os.path.join(OUT, 'sens_B_er_window.png')
    plt.savefig(out, bbox_inches='tight'); plt.close(fig)
    print(f'\n그림 저장: {out}')


def main():
    print('데이터 로딩 (전체 다이)...')
    dies = load_all_dies()
    print(f'로드 완료: {len(dies)}개')

    test_vpi(dies)
    test_er(dies)

    print('\n' + '=' * 70)
    print(' 모든 sensitivity 테스트 완료')
    print('=' * 70)


if __name__ == '__main__':
    main()
