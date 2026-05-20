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

    # 패널 1: raw spectrum + 검출된 null
    L0, I0 = sweeps[biases[0]]
    distance = max(1, int(1.0 / (L0[1] - L0[0])))
    peaks, _ = find_peaks(-I0, prominence=10, distance=distance)
    deep_idx = [int(p) for p in peaks if I0[p] < -25]
    deep_lams = sorted([float(L0[p]) for p in deep_idx])

    # 비교용: 정상 측정 (06-03) 의 dλ/dV 도 같이 보여주기 위해 첫 정상 다이도 로드
    normal_target = None
    for fp, d in items:
        if d == '2019-06-03':
            tdie = load_die(fp)
            if tdie is not None and tdie['wafer'] == die['wafer'] and tdie['band'] == die['band']:
                normal_target = tdie; break
    if normal_target is None:
        for fp, d in items:
            if d == '2019-06-03':
                tdie = load_die(fp)
                if tdie is not None and tdie['band'] == die['band']:
                    normal_target = tdie; break

    # FSR
    if len(deep_lams) >= 2:
        fsr = float(np.median(np.diff(deep_lams)))
    else:
        fsr = float('nan')

    # ── Figure 구성: 위 1개(raw spectrum), 아래는 null 별 분리 ──
    n_nulls_show = min(len(deep_lams), 4)
    if n_nulls_show == 0:
        n_nulls_show = 1
    fig = plt.figure(figsize=(14, 9), dpi=120)
    gs = fig.add_gridspec(2, n_nulls_show, height_ratios=[1, 1.2])

    # ── 상단: raw spectrum (가로 전체) ──
    ax_top = fig.add_subplot(gs[0, :])
    ax_top.plot(L0, I0, lw=0.7, color='steelblue', label=f'V={biases[0]:+.1f}V')
    for p in peaks:
        ax_top.axvline(L0[p], color='gray', alpha=0.3, lw=0.5)
    for p in deep_idx:
        ax_top.axvline(L0[p], color='red', alpha=0.7, lw=1.0)
    ax_top.set_xlabel('λ (nm)'); ax_top.set_ylabel('IL (dB)')
    ax_top.set_title(f'#1  Raw spectrum @ V={biases[0]:+.1f}V  —  '
                     f'detected nulls (red = "deep", I<-25 dB)   '
                     f'[FSR ≈ {fsr:.3f} nm]')
    ax_top.legend(); ax_top.grid(alpha=0.3)

    # ── 하단: null 별 개별 패널 (Δλ 로 표시, 정상 ref overlay) ──
    half = min(0.4, fsr * 0.35) if not np.isnan(fsr) else 0.4

    for i in range(n_nulls_show):
        ax = fig.add_subplot(gs[1, i])
        if i >= len(deep_lams):
            ax.axis('off'); continue
        lam0 = deep_lams[i]

        # 망가진 측정 (현재 die)
        positions = np.array([_parabolic_null(*sweeps[v], lam0, half)
                              for v in rev_biases])
        valid = ~np.isnan(positions)
        if valid.sum() >= 2:
            d_lam_pm = (positions - lam0) * 1000   # nm → pm
            vb = np.array(rev_biases)
            ax.plot(vb[valid], d_lam_pm[valid], 'o-', color='crimson',
                    ms=6, lw=1.5, label='2019-05-31 (broken)')
            # 선형 fit slope
            s, _ = np.polyfit(vb[valid], positions[valid], 1)
            ax.text(0.04, 0.95, f'slope = {s*1000:+.2f} pm/V',
                    transform=ax.transAxes, va='top',
                    color='crimson', fontsize=9, fontweight='bold')

        # 정상 측정 (가능하면 overlay)
        if normal_target is not None:
            n_sweeps = normal_target['sweeps']
            n_rev = sorted([v for v in n_sweeps if v <= 0.0])
            n_positions = np.array([_parabolic_null(*n_sweeps[v], lam0, half)
                                    for v in n_rev])
            n_valid = ~np.isnan(n_positions)
            if n_valid.sum() >= 2:
                d_pm = (n_positions - lam0) * 1000
                ax.plot(np.array(n_rev)[n_valid], d_pm[n_valid], 's--',
                        color='steelblue', ms=5, lw=1.2, alpha=0.85,
                        label='2019-06-03 (normal)')
                s2, _ = np.polyfit(np.array(n_rev)[n_valid],
                                    n_positions[n_valid], 1)
                ax.text(0.04, 0.85, f'slope = {s2*1000:+.2f} pm/V',
                        transform=ax.transAxes, va='top',
                        color='steelblue', fontsize=9, fontweight='bold')

        ax.axhline(0, color='gray', lw=0.6, ls=':')
        ax.set_xlabel('V_bias (V)'); ax.set_ylabel('Δλ_null (pm)')
        ax.set_title(f'null #{i+1}  @ {lam0:.2f} nm', fontsize=10)
        ax.grid(alpha=0.3)
        if i == 0:
            ax.legend(fontsize=8, loc='lower left')

    fig.suptitle('Null tracking across bias — broken (red) vs normal (blue)',
                 fontsize=12, fontweight='bold', y=0.99)
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
    hdr('#7  IV 곡선 시각화 (semi-log)')
    items = find_all_xmls_with_dates(DATA_ROOT)
    # 한 wafer 의 다이 4개 IV 곡선
    targets = []
    for tag_w in ['D07', 'D08', 'D23', 'D24']:
        for fp, d in items:
            if f'/{tag_w}/' in fp:
                targets.append((tag_w, fp))
                break

    # 위: 선형 / 아래: |I| log scale  — 2행 × N열
    fig, axes = plt.subplots(2, len(targets),
                             figsize=(4.0 * len(targets), 7), dpi=120,
                             squeeze=False)
    fig.suptitle('#7  IV curves — linear (top) and |I| log (bottom)\n'
                 'reverse bias V<0: near-linear / forward V>0: exponential',
                 fontsize=12, fontweight='bold', y=0.995)

    for col_idx, (w, fp) in enumerate(targets):
        die = load_die(fp)
        if die is None or die['iv_V'] is None:
            for r in range(2):
                axes[r][col_idx].set_title(f'{w}: no IV'); axes[r][col_idx].axis('off')
            continue
        V, I = die['iv_V'], die['iv_I']

        # 진단 통계
        I_at_neg2 = float(np.interp(-2.0, V, I)) if V.min() <= -2 else float('nan')
        I_at_pos05 = float(np.interp(0.5, V, I)) if V.max() >= 0.5 else float('nan')

        # 상단: 선형 (단위 자동 결정 — μA 또는 mA)
        max_abs = float(np.max(np.abs(I)))
        if max_abs < 1e-3:   # < 1 mA → μA 표시
            unit_scale, unit_label = 1e6, 'I (μA)'
        else:
            unit_scale, unit_label = 1e3, 'I (mA)'
        ax = axes[0][col_idx]
        ax.plot(V, I * unit_scale, '-o', ms=3, lw=1, color='steelblue')
        ax.axhline(0, color='gray', lw=0.5); ax.axvline(0, color='gray', lw=0.5)
        ax.set_xlabel('V (V)'); ax.set_ylabel(unit_label)
        ax.set_title(f'{w}  ({die["band"]}-band, R={die["row"]} C={die["col"]})',
                     fontsize=10, fontweight='bold')
        ax.grid(alpha=0.3)
        ax.text(0.04, 0.95,
                f'I @ -2V: {I_at_neg2*1e9:.2f} nA\n'
                f'I @ +0.5V: {I_at_pos05*unit_scale:.2f} {unit_label[-3:-1]}',
                transform=ax.transAxes, va='top', fontsize=9,
                bbox=dict(facecolor='white', alpha=0.85))

        # 하단: |I| log
        ax = axes[1][col_idx]
        I_abs = np.abs(I)
        # 0 이나 매우 작은 값 floor 처리
        floor = max(1e-12, I_abs[I_abs > 0].min() if (I_abs > 0).any() else 1e-12)
        I_clip = np.clip(I_abs, floor, None)
        # 역바이어스(V<0) 빨강, 순방향(V>=0) 파랑
        rev_mask = V < 0
        fwd_mask = V >= 0
        ax.semilogy(V[rev_mask], I_clip[rev_mask], 'o-', color='crimson',
                    ms=4, lw=1, label='reverse (V<0)')
        ax.semilogy(V[fwd_mask], I_clip[fwd_mask], 's-', color='steelblue',
                    ms=4, lw=1, label='forward (V≥0)')
        ax.axvline(0, color='gray', lw=0.5)
        ax.set_xlabel('V (V)'); ax.set_ylabel('|I| (A, log)')
        ax.grid(alpha=0.3, which='both')
        ax.legend(fontsize=8, loc='lower right')
        # forward 영역의 기대 기울기 (이상적: q/kT ≈ 38.9 /V at 300K)
        # 즉 log10(I) 가 V 당 약 +16.9 (이상계수 n=1, kT/q=25.85mV)
        # 우리 데이터에서 fit slope 계산
        if fwd_mask.sum() >= 2 and (I_clip[fwd_mask] > floor * 10).any():
            try:
                slope_log10 = np.polyfit(V[fwd_mask], np.log10(I_clip[fwd_mask]), 1)[0]
                ax.text(0.04, 0.95,
                        f'fwd slope ≈ {slope_log10:.1f} dec/V\n'
                        f'(ideal Si: ~16.9)',
                        transform=ax.transAxes, va='top', fontsize=9,
                        bbox=dict(facecolor='white', alpha=0.85))
            except Exception:
                pass

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out_path = os.path.join(OUT, '07_iv_curves.png')
    plt.savefig(out_path, bbox_inches='tight')
    plt.close(fig)
    print(f'그림 저장: {out_path}')

    # 진단: 누설전류 / forward 기울기 통계
    print('\n[다이별 IV 진단 — 처음 한 다이씩만]')
    for w, fp in targets:
        die = load_die(fp)
        if die is None or die['iv_V'] is None: continue
        V, I = die['iv_V'], die['iv_I']
        I_neg2 = float(np.interp(-2.0, V, I)) if V.min() <= -2 else float('nan')
        I_pos05 = float(np.interp(0.5, V, I)) if V.max() >= 0.5 else float('nan')
        print(f'  {w}: I@-2V={I_neg2*1e9:+8.2f} nA, '
              f'I@+0.5V={I_pos05*1e6:+8.2f} μA  '
              f'(ratio={abs(I_pos05/I_neg2):.0e} if linear scale)')


# ──────────────────────────────────────────────────────────────────────
#  #8 — FSR 일관성 맵
# ──────────────────────────────────────────────────────────────────────
def _contour_map_panel(ax, sub, value_col, vlim, cmap='viridis',
                       show_numbers=False, number_fmt='.2f'):
    """Delaunay 삼각분할 등고선 + 실측 점 (옵션 숫자)."""
    import matplotlib.tri as mtri
    x = sub['Col'].to_numpy(dtype=float)
    y = sub['Row'].to_numpy(dtype=float)
    z = sub[value_col].to_numpy(dtype=float)
    valid = ~np.isnan(z)
    x, y, z = x[valid], y[valid], z[valid]
    if len(x) < 3:
        return None
    tri = mtri.Triangulation(x, y)
    levels = np.linspace(vlim[0], vlim[1], 20)
    tcf = ax.tricontourf(tri, z, levels=levels, cmap=cmap,
                         vmin=vlim[0], vmax=vlim[1], extend='both')
    ax.tricontour(tri, z, levels=levels, colors='white',
                  linewidths=0.3, alpha=0.4)
    ax.scatter(x, y, c='black', s=22, ec='white', lw=0.5, zorder=5)
    if show_numbers:
        for xi, yi, zv in zip(x, y, z):
            ax.text(xi, yi + 0.4, f'{zv:{number_fmt}}',
                    ha='center', va='bottom', fontsize=7,
                    color='black',
                    bbox=dict(facecolor='white', alpha=0.7,
                              edgecolor='none', pad=0.5))
    ax.set_xlim(-6, 6); ax.set_ylim(-6, 6)
    ax.set_aspect('equal'); ax.grid(alpha=0.2)
    return tcf


def investigate_fsr_map():
    hdr('#8  FSR 일관성 맵 (contour)')
    df = pd.read_csv(os.path.join(PROGRAM_ROOT, 'res', 'csv', 'data.csv'))
    pairs = sorted({(w, b) for w, b in zip(df['Wafer'], df['Band'])},
                   key=lambda x: (x[1] != 'C', x[0]))
    n = len(pairs)
    fig, axes = plt.subplots(1, n, figsize=(3.6 * n, 4.5), dpi=120)
    if n == 1: axes = [axes]
    fig.suptitle('#8  FSR map (continuous surface) — should be approx. constant per wafer',
                 fontsize=12, fontweight='bold', y=1.02)
    # 같은 밴드끼리 같은 컬러 스케일 사용 (C-band 14nm vs O-band 9nm 가 섞이지 않게)
    last = None
    for ax, (w, b) in zip(axes, pairs):
        sub = df[(df['Wafer'] == w) & (df['Band'] == b)]
        band_df = df[df['Band'] == b]['FSR_nm']
        vmin = float(band_df.min()); vmax = float(band_df.max())
        tcf = _contour_map_panel(ax, sub, 'FSR_nm', (vmin, vmax),
                                 cmap='viridis', show_numbers=True,
                                 number_fmt='.2f')
        if tcf is not None: last = tcf
        ax.set_title(f'{w} [{b}-band]\nFSR ≈ {sub["FSR_nm"].mean():.2f} nm',
                     fontsize=10, fontweight='bold')
        ax.set_xlabel('Col'); ax.set_ylabel('Row')
        # 패널별 컬러바 (밴드마다 스케일 다르므로)
        plt.colorbar(tcf, ax=ax, shrink=0.85, pad=0.02, label='FSR (nm)')
    plt.tight_layout()
    out_path = os.path.join(OUT, '08_fsr_map.png')
    plt.savefig(out_path, bbox_inches='tight')
    plt.close(fig)
    print(f'그림 저장: {out_path}')
    print('\n[FSR 통계 (per wafer × band)]')
    print(df.groupby(['Wafer', 'Band'])['FSR_nm'].agg(['mean', 'std', 'min', 'max']).round(3))


# ──────────────────────────────────────────────────────────────────────
#  #9 — dlam_dV 웨이퍼맵
# ──────────────────────────────────────────────────────────────────────
def investigate_dlam_map():
    hdr('#9  dλ/dV 웨이퍼맵 (contour)')
    df = pd.read_csv(os.path.join(PROGRAM_ROOT, 'res', 'csv', 'data.csv'))
    pairs = sorted({(w, b) for w, b in zip(df['Wafer'], df['Band'])},
                   key=lambda x: (x[1] != 'C', x[0]))
    n = len(pairs)
    fig, axes = plt.subplots(1, n, figsize=(3.6 * n, 4.5), dpi=120)
    if n == 1: axes = [axes]
    fig.suptitle('#9  dλ/dV map (continuous surface, pm/V) — values near 0 → V_π runaway',
                 fontsize=12, fontweight='bold', y=1.02)
    vals = df['dlam_dV_pm_per_V'].dropna()
    # 전체 데이터 기준으로 스케일 통일 (음수 영역만 의미 있으므로 음수 범위)
    vmin, vmax = float(vals.min()), float(vals.max())
    # 0 중심 발산형 스케일이 적절하지만 우리 데이터는 모두 음수라 viridis 가 나음
    # 양수가 있으면 RdBu_r, 모두 음수면 viridis_r (어두울수록 양호)
    cmap = 'RdBu_r' if (vmin < 0 < vmax) else 'viridis_r'

    last = None
    for ax, (w, b) in zip(axes, pairs):
        sub = df[(df['Wafer'] == w) & (df['Band'] == b)]
        tcf = _contour_map_panel(ax, sub, 'dlam_dV_pm_per_V', (vmin, vmax),
                                 cmap=cmap, show_numbers=False)
        if tcf is not None: last = tcf
        ax.set_title(f'{w} [{b}-band]\n'
                     f'mean = {sub["dlam_dV_pm_per_V"].mean():+.1f} pm/V',
                     fontsize=10, fontweight='bold')
        ax.set_xlabel('Col'); ax.set_ylabel('Row')
    if last is not None:
        fig.colorbar(last, ax=axes, shrink=0.85, pad=0.02,
                     label='dλ/dV (pm/V)')
    plt.tight_layout()
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
