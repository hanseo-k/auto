"""Raw spectrum + ALIGN reference 시각화.

각 (Wafer, Band) 그룹에 대해 두 종류 그림을 생성한다.

1. raw_spectra_<wafer>_<band>.png
   해당 그룹의 모든 다이 스펙트럼을 한 패널에 overlay.
   V = -2V 의 raw transmission 만 보여서 다이 간 산포를 한눈에 비교.
   ALIGN 참조 (커플러 + 도파로 손실) 도 같은 그리드에서 두꺼운 회색선으로 overlay.

2. raw_spectra_detail_<wafer>_<band>_R<r>_C<c>.png  (대표 다이 1 개)
   해당 그룹에서 R^2 가 가장 높은 다이를 골라 모든 바이어스 (6 점) 의
   스펙트럼을 색 변화로 그림. 같이 ALIGN 참조와 T_dev (ALIGN 차감) 도 표시.

3. align_reference_summary.png
   5 개 (Wafer, Band) 그룹의 ALIGN 참조 스펙트럼을 한 패널에 overlay
   하여 커플러/도파로 손실의 그룹 간 차이를 보여준다.

출력: res/diagnostics/
"""
import sys, os
PROGRAM_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROGRAM_ROOT, 'src'))

import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap

from xml_loader import load_die, find_all_xmls, t_dev
from extract_er import ER_WINDOW_NM
from extract_vpi import extract_vpi


DATA_ROOT = os.path.join(PROGRAM_ROOT, 'data', 'HY202103')
OUT = os.path.join(PROGRAM_ROOT, 'res', 'diagnostics')
os.makedirs(OUT, exist_ok=True)


def _group_band_color(band):
    return '#1f77b4' if band == 'C' else '#d62728'


def _load_group_dies(items):
    """전체 XML 을 (Wafer, Band) -> [die,...] 로 정리."""
    groups = {}
    for fp in items:
        d = load_die(fp)
        if d is None:
            continue
        key = (d['wafer'], d['band'])
        groups.setdefault(key, []).append(d)
    return groups


# ──────────────────────────────────────────────────────────────────────
# Figure 1: 그룹별 raw spectra overlay (다이별 변동)
# ──────────────────────────────────────────────────────────────────────
def plot_group_overlay(group_dies, group_key, save_path):
    wafer, band = group_key
    win_lo, win_hi = ER_WINDOW_NM[band]

    fig, ax = plt.subplots(figsize=(11, 5.5), dpi=140)

    # 다이별 V=-2V 스펙트럼 overlay
    cmap = get_cmap('viridis')
    n = len(group_dies)
    for i, die in enumerate(sorted(group_dies, key=lambda d: (d['row'], d['col']))):
        bias = sorted(die['sweeps'].keys())[0]   # 가장 negative bias
        L, IL = die['sweeps'][bias]
        ax.plot(L, IL, lw=0.5, alpha=0.6,
                color=cmap(i / max(n - 1, 1)),
                label=f'R={die["row"]:+d} C={die["col"]:+d}')

    # ALIGN 참조 — 그룹의 첫 다이의 ALIGN 을 대표로 (다이 간 ALIGN 은 거의 동일)
    ref_die = group_dies[0]
    if ref_die['ref_L'] is not None:
        ax.plot(ref_die['ref_L'], ref_die['ref_IL'],
                color='black', lw=2.2, linestyle='--',
                label='ALIGN reference (coupler + waveguide)', zorder=10)

    # ER window 음영
    ax.axvspan(win_lo, win_hi, color='lightgreen', alpha=0.12,
               label=f'ER window [{win_lo:.0f}-{win_hi:.0f}] nm')

    ax.set_xlabel('λ (nm)')
    ax.set_ylabel('Measured Transmission (dB)')
    ax.set_title(f'{wafer} [{band}-band] — Raw spectra (V={bias:+.1f}V) + ALIGN reference\n'
                 f'(다이 {n}개 overlay; ALIGN = 점선 검정)',
                 fontsize=11, fontweight='bold')
    ax.grid(alpha=0.3)
    # legend: ALIGN 과 ER window 만 메인으로 표시, 다이 라벨은 너무 많음
    handles, labels = ax.get_legend_handles_labels()
    main_h, main_l = [], []
    for h, l in zip(handles, labels):
        if 'ALIGN' in l or 'ER window' in l:
            main_h.append(h); main_l.append(l)
    ax.legend(main_h, main_l, fontsize=9, loc='lower left',
              framealpha=0.85)
    plt.tight_layout()
    plt.savefig(save_path, bbox_inches='tight')
    plt.close(fig)
    print(f'Group overlay saved: {save_path}')


# ──────────────────────────────────────────────────────────────────────
# Figure 2: 대표 다이 상세 — 모든 바이어스 + ALIGN ref + T_dev
# ──────────────────────────────────────────────────────────────────────
def plot_detail(die, save_path):
    wafer, band = die['wafer'], die['band']
    row, col = die['row'], die['col']
    win_lo, win_hi = ER_WINDOW_NM[band]
    sweeps = die['sweeps']
    biases = sorted(sweeps.keys())

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 9), dpi=140, sharex=True)
    fig.suptitle(f'{wafer} [{band}-band] R={row:+d} C={col:+d}  '
                 f'— raw spectra (top) and ALIGN-subtracted T_dev (bottom)',
                 fontsize=12, fontweight='bold', y=0.995)

    # ── 상단: raw spectra + ALIGN reference ──
    cmap = get_cmap('coolwarm')
    n_b = len(biases)
    for i, V in enumerate(biases):
        L, IL = sweeps[V]
        color = cmap(i / max(n_b - 1, 1))
        ax1.plot(L, IL, lw=1.0, color=color, label=f'V={V:+.1f} V')

    if die['ref_L'] is not None:
        ax1.plot(die['ref_L'], die['ref_IL'], color='black', lw=2.0,
                 linestyle='--',
                 label='ALIGN reference', zorder=10)

    ax1.axvspan(win_lo, win_hi, color='lightgreen', alpha=0.12)
    ax1.set_ylabel('Measured Transmission (dB)')
    ax1.set_title('Raw spectra (전체 바이어스) + ALIGN reference', fontsize=10)
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=8, loc='lower left', ncol=2, framealpha=0.85)

    # ── 하단: T_dev = MZM - ALIGN (커플러 영향 제거된 device transfer) ──
    if die['ref_L'] is not None:
        for i, V in enumerate(biases):
            L, IL = sweeps[V]
            T = t_dev(L, IL, die['ref_L'], die['ref_IL'])
            color = cmap(i / max(n_b - 1, 1))
            ax2.plot(L, T, lw=1.0, color=color, label=f'V={V:+.1f} V')
    ax2.axhline(0, color='gray', lw=0.5)
    ax2.axvspan(win_lo, win_hi, color='lightgreen', alpha=0.12)
    ax2.set_xlabel('λ (nm)')
    ax2.set_ylabel('T_dev = IL_mzm - IL_ref (dB)')
    ax2.set_title('ALIGN 차감 후 device transfer function (커플러 영향 제거)', fontsize=10)
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=8, loc='lower left', ncol=2, framealpha=0.85)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(save_path, bbox_inches='tight')
    plt.close(fig)
    print(f'Detail plot saved: {save_path}')


# ──────────────────────────────────────────────────────────────────────
# Figure 3: ALIGN reference 그룹 간 비교
# ──────────────────────────────────────────────────────────────────────
def plot_align_summary(groups, save_path):
    """5 개 (Wafer, Band) 의 ALIGN 참조를 한 패널에 overlay.
    커플러 + 도파로 손실의 그룹 간 차이를 보여준다."""
    fig, ax = plt.subplots(figsize=(11, 5.5), dpi=140)
    colors = {'C': '#1f77b4', 'O': '#d62728'}
    linestyle_map = {}  # wafer 별로 점선/실선/...
    line_styles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1))]

    keys = sorted(groups.keys(), key=lambda k: (k[1] != 'C', k[0]))
    for i, key in enumerate(keys):
        wafer, band = key
        rep = groups[key][0]
        if rep['ref_L'] is None: continue
        ls = line_styles[i % len(line_styles)]
        ax.plot(rep['ref_L'], rep['ref_IL'],
                lw=1.5, ls=ls, color=colors[band],
                label=f'{wafer} [{band}]')

    ax.set_xlabel('λ (nm)')
    ax.set_ylabel('ALIGN Transmission (dB)')
    ax.set_title('ALIGN reference 비교 — 그룹별 커플러 + 도파로 손실 차이',
                 fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=10, loc='lower center', ncol=5, framealpha=0.85)
    plt.tight_layout()
    plt.savefig(save_path, bbox_inches='tight')
    plt.close(fig)
    print(f'ALIGN summary saved: {save_path}')


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────
def main():
    print('=' * 70)
    print(' Raw spectra + ALIGN reference 그래프 생성')
    print('=' * 70)

    items = find_all_xmls(DATA_ROOT)
    print(f'발견된 XML: {len(items)}개')
    groups = _load_group_dies(items)
    print(f'(Wafer, Band) 그룹: {len(groups)}개')

    # 그룹별 raw spectra overlay
    for key in sorted(groups.keys(), key=lambda k: (k[1] != 'C', k[0])):
        wafer, band = key
        save_path = os.path.join(OUT, f'raw_spectra_{wafer}_{band}.png')
        plot_group_overlay(groups[key], key, save_path)

    # 대표 다이 상세 — 각 그룹마다 가장 깔끔한 다이 (Row, Col 가운데에 가까운)
    for key in sorted(groups.keys(), key=lambda k: (k[1] != 'C', k[0])):
        wafer, band = key
        rep = min(groups[key], key=lambda d: abs(d['row']) + abs(d['col']))
        save_path = os.path.join(OUT,
            f'raw_spectra_detail_{wafer}_{band}_R{rep["row"]:+d}_C{rep["col"]:+d}.png')
        plot_detail(rep, save_path)

    # ALIGN reference 비교
    save_path = os.path.join(OUT, 'align_reference_summary.png')
    plot_align_summary(groups, save_path)

    print('\n저장 완료:', OUT)


if __name__ == '__main__':
    main()
