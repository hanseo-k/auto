"""9개 데이터 품질/방법론 진단 한번에 실행.

각 항목별 결과를 stdout 으로 출력하고 진단 그림을 doc/investigation/ 에 저장.

실행:
    python3 src/investigate.py
"""
import sys, os, re, glob
PROGRAM_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROGRAM_ROOT, 'src'))

import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from xml_loader import load_die, BAND_OF_FILE
from extract_er import extract_er
from extract_il import extract_il
from extract_vpi import extract_vpi, _parabolic_null
from analyze_by_date import find_all_xmls_with_dates


DATA_ROOT = '/Users/gimhanseo/Desktop/공프/HY202103'
OUT = os.path.join(PROGRAM_ROOT, 'doc', 'investigation')
os.makedirs(OUT, exist_ok=True)


def hdr(title):
    print('\n' + '=' * 70)
    print(f' {title}')
    print('=' * 70)


# ──────────────────────────────────────────────────────────────────────
#  #1 — 2019-05-31 V_π 폭주 원인 진단
# ──────────────────────────────────────────────────────────────────────
def investigate_05_31():
    hdr('#1  2019-05-31 V_π 폭주 원인 진단')
    items = find_all_xmls_with_dates(DATA_ROOT)
    # 05-31 측정 파일들
    bad_files = [(fp, d) for fp, d in items if d == '2019-05-31']
    print(f'2019-05-31 측정 파일: {len(bad_files)}개')

    # 한 다이 골라서 raw 스펙트럼 + null 트래킹 visualize
    target = bad_files[0][0]
    print(f'대표 파일: {os.path.basename(target)}')
    die = load_die(target)
    if die is None:
        print('파싱 실패')
        return
    sweeps = die['sweeps']
    biases = sorted(sweeps.keys())
    rev_biases = [v for v in biases if v <= 0.0]

    fig, axes = plt.subplots(2, 1, figsize=(12, 9), dpi=120)

    # 패널 1: V=가장 음의 바이어스에서의 raw + null peaks
    L0, I0 = sweeps[biases[0]]
    distance = max(1, int(1.0 / (L0[1] - L0[0])))
    peaks, _ = find_peaks(-I0, prominence=10, distance=distance)
    deep = [int(p) for p in peaks if I0[p] < -25]

    ax = axes[0]
    ax.plot(L0, I0, lw=0.7, color='steelblue', label=f'V={biases[0]:+.1f}V')
    for p in peaks:
        ax.axvline(L0[p], color='gray', alpha=0.3, lw=0.5)
    for p in deep:
        ax.axvline(L0[p], color='red', alpha=0.7, lw=1.0)
    ax.set_xlabel('λ (nm)'); ax.set_ylabel('IL (dB)')
    ax.set_title(f'#1  Raw spectrum @ V={biases[0]:+.1f}V — '
                 f'detected nulls (red = "deep", I<-25dB)')
    ax.legend(); ax.grid(alpha=0.3)

    # FSR 계산
    if len(deep) >= 2:
        deep_lams = sorted([float(L0[p]) for p in deep])
        fsr = float(np.median(np.diff(deep_lams)))
        ax.text(0.02, 0.95, f'deep nulls: {len(deep)}\nFSR ≈ {fsr:.3f} nm',
                transform=ax.transAxes, va='top', fontsize=10,
                bbox=dict(facecolor='white', alpha=0.85))
    else:
        fsr = np.nan
        ax.text(0.02, 0.95, f'⚠ deep null < 2 — FSR 추출 불가',
                transform=ax.transAxes, va='top', color='red', fontsize=11,
                fontweight='bold')

    # 패널 2: 각 null 의 바이어스별 추적
    ax = axes[1]
    if len(deep) >= 2 and not np.isnan(fsr):
        half = min(0.4, fsr * 0.35)
        for lam0 in sorted([float(L0[p]) for p in deep])[:6]:  # 최대 6개
            positions = []
            for v in rev_biases:
                L, I = sweeps[v]
                positions.append(_parabolic_null(L, I, lam0, half))
            positions = np.array(positions)
            valid = ~np.isnan(positions)
            if valid.sum() < 2:
                continue
            ax.plot(np.array(rev_biases)[valid], positions[valid],
                    marker='o', lw=1.2, ms=5, label=f'null₀={lam0:.2f}nm')
            # 점프 감지
            total_shift = abs(positions[valid][-1] - positions[valid][0])
            if total_shift > half * 1.5:
                ax.annotate('JUMP!', xy=(rev_biases[-1], positions[valid][-1]),
                            color='red', fontweight='bold', fontsize=11)

    ax.set_xlabel('V_bias (V)'); ax.set_ylabel('null wavelength (nm)')
    ax.set_title('Per-null tracking across reverse-bias sweep '
                 '(slope = dλ/dV → V_π = FSR / 2·|slope|)')
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    plt.tight_layout()
    out_path = os.path.join(OUT, '01_2019-05-31_diagnosis.png')
    plt.savefig(out_path, bbox_inches='tight')
    plt.close(fig)
    print(f'그림 저장: {out_path}')

    # 통계: 05-31 의 모든 다이에 대해 deep null 개수
    print('\n[05-31 다이별 deep null 개수]')
    deep_counts = []
    for fp, _ in bad_files:
        die = load_die(fp)
        if die is None:
            continue
        L0, I0 = die['sweeps'][sorted(die['sweeps'].keys())[0]]
        peaks, _ = find_peaks(-I0, prominence=10,
                              distance=max(1, int(1.0 / (L0[1] - L0[0]))))
        n_deep = sum(1 for p in peaks if I0[p] < -25)
        deep_counts.append(n_deep)
    if deep_counts:
        print(f'  평균 {np.mean(deep_counts):.1f}개, '
              f'min {min(deep_counts)}, max {max(deep_counts)}')
        n_too_few = sum(1 for c in deep_counts if c < 2)
        print(f'  → FSR 추출 불가 (deep < 2): {n_too_few}개')

    # 다른 정상 날짜와 비교
    print('\n[다른 날짜 비교 (deep null 개수)]')
    for date in ['2019-06-03', '2019-07-12']:
        date_files = [(fp, d) for fp, d in items if d == date][:5]
        counts = []
        for fp, _ in date_files:
            die = load_die(fp)
            if die is None:
                continue
            L0, I0 = die['sweeps'][sorted(die['sweeps'].keys())[0]]
            peaks, _ = find_peaks(-I0, prominence=10,
                                  distance=max(1, int(1.0 / (L0[1] - L0[0]))))
            n_deep = sum(1 for p in peaks if I0[p] < -25)
            counts.append(n_deep)
        if counts:
            print(f'  {date}: 평균 {np.mean(counts):.1f}, '
                  f'min {min(counts)}, max {max(counts)}')


# ──────────────────────────────────────────────────────────────────────
#  #2 — D08 반복측정 재현성
# ──────────────────────────────────────────────────────────────────────
def investigate_repeatability():
    hdr('#2  D08 반복측정 재현성 검증')
    csv_path = os.path.join(PROGRAM_ROOT, 'res', 'csv', 'data_by_date.csv')
    df = pd.read_csv(csv_path)

    # D08만, 정상 데이터만 (is_problematic=False)
    d08 = df[(df['Wafer'] == 'D08') & (~df['is_problematic'])].copy()
    if d08.empty:
        print('정상 D08 측정 없음')
        return

    # (Row, Col, Band) 별로 측정 날짜가 2개 이상인 다이만
    g = d08.groupby(['Band', 'Row', 'Col']).size().reset_index(name='n_meas')
    multi = g[g['n_meas'] >= 2]
    print(f'D08 정상측정 중 2회 이상 측정된 다이: {len(multi)}개')

    if len(multi) == 0:
        print('재측정된 다이 없음 — 다른 wafer 확인 필요')
        # 모든 wafer 에 대해 시도
        all_clean = df[~df['is_problematic']]
        g_all = all_clean.groupby(['Wafer', 'Band', 'Row', 'Col']).size().reset_index(name='n_meas')
        multi_all = g_all[g_all['n_meas'] >= 2]
        print(f'모든 wafer 통틀어 2회 이상: {len(multi_all)}개')
        if len(multi_all) > 0:
            print(multi_all.head(10).to_string(index=False))
        return

    # variation 계산
    rows = []
    for _, r in multi.iterrows():
        sub = d08[(d08['Band'] == r['Band']) &
                  (d08['Row'] == r['Row']) &
                  (d08['Col'] == r['Col'])]
        rows.append({
            'Band': r['Band'], 'Row': r['Row'], 'Col': r['Col'],
            'n': len(sub),
            'ER_range': sub['ER_dB'].max() - sub['ER_dB'].min(),
            'IL_range': sub['IL_dB'].max() - sub['IL_dB'].min(),
            'Vpi_range': sub['Vpi_V'].max() - sub['Vpi_V'].min(),
        })
    rep = pd.DataFrame(rows)
    print('\n[다이별 측정값 range (max-min)]')
    print(rep.describe().round(2))
    print('\n  → ER range 평균 {:.2f} dB / max {:.2f} dB'.format(
        rep['ER_range'].mean(), rep['ER_range'].max()))
    print('  → IL range 평균 {:.2f} dB / max {:.2f} dB'.format(
        rep['IL_range'].mean(), rep['IL_range'].max()))
    print('  → Vpi range 평균 {:.2f} V / max {:.2f} V'.format(
        rep['Vpi_range'].mean(), rep['Vpi_range'].max()))


# ──────────────────────────────────────────────────────────────────────
#  #3 — Width_nm 분포
# ──────────────────────────────────────────────────────────────────────
def investigate_width():
    hdr('#3  Width_nm 분포 점검')
    df = pd.read_csv(os.path.join(PROGRAM_ROOT, 'res', 'csv', 'data.csv'))
    print('Width_nm value_counts:')
    print(df['Width_nm'].value_counts(dropna=False).to_string())
    print('\n[Wafer × Band × Width]')
    pivot = df.groupby(['Wafer', 'Band', 'Width_nm']).size().unstack(fill_value=0)
    print(pivot)


# ──────────────────────────────────────────────────────────────────────
#  #4 — dedup 가정 (최신=최선)
# ──────────────────────────────────────────────────────────────────────
def investigate_dedup():
    hdr('#4  dedup "최신=최선" 가정 검증')
    df = pd.read_csv(os.path.join(PROGRAM_ROOT, 'res', 'csv', 'data_by_date.csv'))
    # 같은 (Wafer, Band, Row, Col) 이 여러 날 측정된 경우
    g = df.groupby(['Wafer', 'Band', 'Row', 'Col'])
    rows = []
    for key, sub in g:
        if len(sub) < 2:
            continue
        latest_idx = sub['Date'].idxmax()
        latest_bad = sub.loc[latest_idx, 'is_problematic']
        n_clean_total = (~sub['is_problematic']).sum()
        rows.append({
            'Wafer': key[0], 'Band': key[1], 'Row': key[2], 'Col': key[3],
            'n_meas': len(sub),
            'latest_date': sub.loc[latest_idx, 'Date'],
            'latest_is_problematic': latest_bad,
            'n_clean': n_clean_total,
        })
    rep = pd.DataFrame(rows)
    print(f'다중 측정 다이: {len(rep)}개')
    print(f'  → 최신이 problematic 인 다이: {rep["latest_is_problematic"].sum()}개')
    print(f'  → 최신은 OK 인데 다른 측정이 problematic 인 다이: '
          f'{((~rep["latest_is_problematic"]) & (rep["n_clean"] < rep["n_meas"])).sum()}개')
    print('\n[dedup 가정이 위험한 케이스 (최신이 problematic)]')
    risk = rep[rep['latest_is_problematic']]
    if not risk.empty:
        print(risk.to_string(index=False))
    else:
        print('  없음 — 다행히 dedup 가정 안전')


# ──────────────────────────────────────────────────────────────────────
#  #5 — D08 C/O 양 밴드 관계
# ──────────────────────────────────────────────────────────────────────
def investigate_d08_bands():
    hdr('#5  D08 C/O 양밴드 관계 확인')
    df = pd.read_csv(os.path.join(PROGRAM_ROOT, 'res', 'csv', 'data.csv'))
    d08 = df[df['Wafer'] == 'D08']
    if d08.empty:
        print('D08 데이터 없음')
        return
    print(f'D08 총 측정: {len(d08)}개 (C={len(d08[d08.Band=="C"])}, O={len(d08[d08.Band=="O"])})')
    # 같은 (Row, Col) 이 C, O 둘 다 있는지
    rc_c = set(zip(d08[d08.Band == 'C']['Row'], d08[d08.Band == 'C']['Col']))
    rc_o = set(zip(d08[d08.Band == 'O']['Row'], d08[d08.Band == 'O']['Col']))
    both = rc_c & rc_o
    print(f'(Row, Col) 분포:')
    print(f'  C-band only:   {len(rc_c - rc_o)}개')
    print(f'  O-band only:   {len(rc_o - rc_c)}개')
    print(f'  C, O 둘 다:    {len(both)}개')
    if both:
        print('\n[같은 (Row,Col) 의 C vs O 측정값]')
        for r, c in sorted(both)[:5]:
            cc = d08[(d08.Band == 'C') & (d08.Row == r) & (d08.Col == c)].iloc[0]
            oo = d08[(d08.Band == 'O') & (d08.Row == r) & (d08.Col == c)].iloc[0]
            print(f'  ({r:+d},{c:+d})  C: ER={cc.ER_dB:.1f} IL={cc.IL_dB:.1f} Vpi={cc.Vpi_V:.1f}  '
                  f'|  O: ER={oo.ER_dB:.1f} IL={oo.IL_dB:.1f} Vpi={oo.Vpi_V:.1f}')


# ──────────────────────────────────────────────────────────────────────
#  #6 — ER 윈도우 FSR 일관성
# ──────────────────────────────────────────────────────────────────────
def investigate_er_window():
    hdr('#6  ER 윈도우 / FSR 일관성')
    df = pd.read_csv(os.path.join(PROGRAM_ROOT, 'res', 'csv', 'data.csv'))
    print('[FSR 분포 (Wafer × Band)]')
    print(df.groupby(['Wafer', 'Band'])['FSR_nm'].describe().round(3))
    # 윈도우 폭 14 nm 에 들어가는 null 개수 추정
    print('\n[14nm 윈도우 안 예상 null 개수 (14/FSR)]')
    df['n_nulls_in_win'] = 14.0 / df['FSR_nm']
    print(df.groupby(['Wafer', 'Band'])['n_nulls_in_win'].describe().round(2))


# ──────────────────────────────────────────────────────────────────────
#  #7 — IV 곡선 시각화
# ──────────────────────────────────────────────────────────────────────
def investigate_iv():
    hdr('#7  IV 곡선 시각화')
    items = find_all_xmls_with_dates(DATA_ROOT)
    # 한 wafer 의 다이 4개 IV 곡선
    targets = []
    for tag_w in ['D07', 'D08', 'D23', 'D24']:
        for fp, d in items:
            if f'/{tag_w}/' in fp:
                targets.append((tag_w, fp))
                break
    fig, axes = plt.subplots(1, len(targets), figsize=(4.2 * len(targets), 4), dpi=120)
    if len(targets) == 1:
        axes = [axes]
    fig.suptitle('#7  IV curves — sanity check (leakage / breakdown)',
                 fontsize=13, fontweight='bold', y=1.02)
    for ax, (w, fp) in zip(axes, targets):
        die = load_die(fp)
        if die is None or die['iv_V'] is None:
            ax.set_title(f'{w}: no IV')
            continue
        V, I = die['iv_V'], die['iv_I']
        ax.plot(V, I * 1e3, '-o', ms=3, lw=1)
        ax.axhline(0, color='gray', lw=0.5)
        ax.axvline(0, color='gray', lw=0.5)
        ax.set_xlabel('V (V)'); ax.set_ylabel('I (mA)')
        ax.set_title(f'{w}  ({die["band"]}-band, R={die["row"]} C={die["col"]})',
                     fontsize=10)
        ax.grid(alpha=0.3)
    plt.tight_layout()
    out_path = os.path.join(OUT, '07_iv_curves.png')
    plt.savefig(out_path, bbox_inches='tight')
    plt.close(fig)
    print(f'그림 저장: {out_path}')


# ──────────────────────────────────────────────────────────────────────
#  #8 — FSR 일관성 맵
# ──────────────────────────────────────────────────────────────────────
def investigate_fsr_map():
    hdr('#8  FSR 일관성 맵')
    df = pd.read_csv(os.path.join(PROGRAM_ROOT, 'res', 'csv', 'data.csv'))
    pairs = sorted({(w, b) for w, b in zip(df['Wafer'], df['Band'])},
                   key=lambda x: (x[1] != 'C', x[0]))
    n = len(pairs)
    fig, axes = plt.subplots(1, n, figsize=(3.4 * n, 4), dpi=120)
    if n == 1: axes = [axes]
    fig.suptitle('#8  FSR map  — should be approximately constant per wafer',
                 fontsize=13, fontweight='bold', y=1.02)
    vmin = df['FSR_nm'].min(); vmax = df['FSR_nm'].max()
    last = None
    for ax, (w, b) in zip(axes, pairs):
        sub = df[(df['Wafer'] == w) & (df['Band'] == b)]
        sc = ax.scatter(sub['Col'], sub['Row'], c=sub['FSR_nm'],
                        cmap='viridis', vmin=vmin, vmax=vmax,
                        s=320, marker='s', ec='black', lw=0.6)
        last = sc
        for _, r in sub.iterrows():
            ax.text(r['Col'], r['Row'], f'{r["FSR_nm"]:.2f}',
                    ha='center', va='center', fontsize=7,
                    color='white' if r['FSR_nm'] < (vmin+vmax)/2 else 'black')
        ax.set_xlim(-5.5, 5.5); ax.set_ylim(-5.5, 5.5)
        ax.set_aspect('equal')
        ax.set_title(f'{w} [{b}-band]', fontsize=10, fontweight='bold')
        ax.grid(alpha=0.25)
        ax.set_xlabel('Col'); ax.set_ylabel('Row')
    fig.colorbar(last, ax=axes, shrink=0.85, pad=0.02, label='FSR (nm)')
    out_path = os.path.join(OUT, '08_fsr_map.png')
    plt.savefig(out_path, bbox_inches='tight')
    plt.close(fig)
    print(f'그림 저장: {out_path}')
    # 통계
    print('\n[FSR 통계 (per wafer × band)]')
    print(df.groupby(['Wafer', 'Band'])['FSR_nm'].agg(['mean', 'std', 'min', 'max']).round(3))


# ──────────────────────────────────────────────────────────────────────
#  #9 — dlam_dV 웨이퍼맵
# ──────────────────────────────────────────────────────────────────────
def investigate_dlam_map():
    hdr('#9  dλ/dV 웨이퍼맵')
    df = pd.read_csv(os.path.join(PROGRAM_ROOT, 'res', 'csv', 'data.csv'))
    pairs = sorted({(w, b) for w, b in zip(df['Wafer'], df['Band'])},
                   key=lambda x: (x[1] != 'C', x[0]))
    n = len(pairs)
    fig, axes = plt.subplots(1, n, figsize=(3.4 * n, 4), dpi=120)
    if n == 1: axes = [axes]
    fig.suptitle('#9  dλ/dV map (pm/V)  — values near 0 cause V_π runaway',
                 fontsize=13, fontweight='bold', y=1.02)
    vals = df['dlam_dV_pm_per_V'].dropna()
    vmin, vmax = float(vals.min()), float(vals.max())
    last = None
    for ax, (w, b) in zip(axes, pairs):
        sub = df[(df['Wafer'] == w) & (df['Band'] == b)]
        sc = ax.scatter(sub['Col'], sub['Row'], c=sub['dlam_dV_pm_per_V'],
                        cmap='RdBu_r', vmin=-abs(vmax), vmax=abs(vmax),
                        s=320, marker='s', ec='black', lw=0.6)
        last = sc
        for _, r in sub.iterrows():
            v = r['dlam_dV_pm_per_V']
            if pd.isna(v):
                continue
            ax.text(r['Col'], r['Row'], f'{v:+.1f}', ha='center', va='center',
                    fontsize=7, color='black')
        ax.set_xlim(-5.5, 5.5); ax.set_ylim(-5.5, 5.5)
        ax.set_aspect('equal')
        ax.set_title(f'{w} [{b}-band]', fontsize=10, fontweight='bold')
        ax.grid(alpha=0.25)
        ax.set_xlabel('Col'); ax.set_ylabel('Row')
    fig.colorbar(last, ax=axes, shrink=0.85, pad=0.02, label='dλ/dV (pm/V)')
    out_path = os.path.join(OUT, '09_dlam_dv_map.png')
    plt.savefig(out_path, bbox_inches='tight')
    plt.close(fig)
    print(f'그림 저장: {out_path}')
    print('\n[dλ/dV 통계]')
    print(df.groupby(['Wafer', 'Band'])['dlam_dV_pm_per_V'].agg(['mean', 'std', 'min', 'max']).round(2))


def main():
    investigate_05_31()
    investigate_repeatability()
    investigate_width()
    investigate_dedup()
    investigate_d08_bands()
    investigate_er_window()
    investigate_iv()
    investigate_fsr_map()
    investigate_dlam_map()
    print('\n' + '=' * 70)
    print(' 모든 진단 완료')
    print(f' 그림 저장 위치: {OUT}')
    print('=' * 70)


if __name__ == '__main__':
    main()
