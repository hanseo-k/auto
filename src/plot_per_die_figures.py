"""다이별 6 종류 분석 그림 생성 — PE2_Lec06 Objective 6 (page 30).

각 다이 (Date, Band, Wafer, Row, Col) 마다 다음 6 개 그림을 생성한다.

    01_mzm_ref_spectra.png   : 모든 바이어스 raw spectra + ALIGN reference
    02_ref_polyfit.png       : ALIGN spectrum 의 2~6 차 polynomial fit
    03_flat_transmission.png : 2-step flatten (envelope subtraction) 후 spectra
    04_mzi_fit.png           : 한 바이어스에서 MZI cos² 모델 fit
    05_iv_semilog.png        : IV semilog plot (|I|)
    06_iv_fit.png            : IV reverse polynomial + forward Shockley fit

폴더 구조:
    res/figures_per_die/
        <Date>/
            <Band>-band/
                <Wafer>/
                    R<row>_C<col>/
                        01_*.png ~ 06_*.png

실행:
    python3 src/plot_per_die_figures.py
"""
import sys, os, re, glob
PROGRAM_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROGRAM_ROOT, 'src'))

import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from scipy.optimize import curve_fit
import multiprocessing

from xml_loader import load_die, BAND_OF_FILE
from extract_er import ER_WINDOW_NM


DATA_ROOT = os.path.join(PROGRAM_ROOT, 'data', 'HY202103')
OUT_ROOT  = os.path.join(PROGRAM_ROOT, 'res', 'figures_per_die')


# ──────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────
def _r2(y, y_fit):
    """R² 계산."""
    ss_res = float(np.sum((y - y_fit) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')


def _find_peaks(y, window=14):
    """국소 윈도우 기반 peak 검출."""
    y = np.asarray(y)
    n = len(y)
    peaks = []
    for i in range(window, n - window):
        if y[i] == np.max(y[i - window: i + window + 1]):
            peaks.append(i)
    return np.array(peaks, dtype=int)


def _mzi_model(wl, A, B, fsr, phi, wl0):
    """T(λ) = A + B · cos²(π(λ - λ₀)/FSR + φ)"""
    return A + B * np.cos(np.pi * (wl - wl0) / fsr + phi) ** 2


def _fit_AB(wl, T, fsr, phi, wl0):
    """주어진 (fsr, phi, wl0) 에서 A, B 는 closed-form (선형) 으로 fit."""
    c2 = np.cos(np.pi * (wl - wl0) / fsr + phi) ** 2
    M = np.column_stack([np.ones_like(wl), c2])
    coef, *_ = np.linalg.lstsq(M, T, rcond=None)
    return float(coef[0]), float(coef[1])


def _shockley(V, Is, n):
    """이상 다이오드: I = I_s · (exp(V/(n·V_T)) − 1)"""
    VT = 0.02585
    return Is * (np.exp(np.clip(V / (n * VT), -100, 100)) - 1)


def _find_all_xmls_with_meta(data_root):
    """경로에서 (xml_path, date, wafer, row, col, band) 를 수집."""
    items = []
    for tag, band in BAND_OF_FILE.items():
        pattern = os.path.join(data_root, '**', f'*_DCM_{tag}.xml')
        for fp in sorted(glob.glob(pattern, recursive=True)):
            fname = os.path.basename(fp)
            m_date = re.search(r'(\d{4})(\d{2})(\d{2})_\d{6}', fp)
            m_rc = re.search(r'(D\d+)_\((-?\d+),(-?\d+)\)', fname)
            if not (m_date and m_rc):
                continue
            date_str = f'{m_date.group(1)}-{m_date.group(2)}-{m_date.group(3)}'
            wafer = m_rc.group(1)
            row = int(m_rc.group(2)); col = int(m_rc.group(3))
            items.append((fp, date_str, wafer, row, col, band))
    return items


# ──────────────────────────────────────────────────────────────────────
# Plot 01 — MZM raw spectra + ALIGN reference
# ──────────────────────────────────────────────────────────────────────
def plot_01_mzm_ref_spectra(die, save_path):
    sweeps = die['sweeps']
    biases = sorted(sweeps.keys())
    cmap = get_cmap('coolwarm')
    fig, ax = plt.subplots(figsize=(8.5, 5), dpi=120)
    for i, V in enumerate(biases):
        L, IL = sweeps[V]
        ax.plot(L, IL, lw=1.0, color=cmap(i / max(len(biases) - 1, 1)),
                label=f'MZM {V:+.1f} V')
    if die['ref_L'] is not None:
        ax.plot(die['ref_L'], die['ref_IL'], color='black', lw=2.0,
                label='Reference')
    ax.set_xlabel('Wavelength [nm]')
    ax.set_ylabel('Measured transmission [dB]')
    ax.set_title('MZM and Reference spectra', fontsize=11, fontweight='bold')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc='lower center', ncol=4, framealpha=0.85)
    plt.tight_layout()
    plt.savefig(save_path, bbox_inches='tight')
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────
# Plot 02 — ALIGN reference polynomial fit (2~6 차)
# ──────────────────────────────────────────────────────────────────────
def plot_02_ref_polyfit(die, save_path):
    if die['ref_L'] is None:
        return
    ref_L, ref_IL = die['ref_L'], die['ref_IL']
    wl_mean = float(np.mean(ref_L))

    fig, ax = plt.subplots(figsize=(8.5, 5), dpi=120)
    ax.scatter(ref_L, ref_IL, s=6, color='C0', alpha=0.4, label='Reference data')

    cmap = get_cmap('viridis')
    degrees = [2, 3, 4, 5, 6]
    for i, deg in enumerate(degrees):
        coef = np.polyfit(ref_L - wl_mean, ref_IL, deg)
        fit = np.polyval(coef, ref_L - wl_mean)
        r2 = _r2(ref_IL, fit)
        ax.plot(ref_L, fit, lw=1.5, color=cmap(i / (len(degrees) - 1)),
                label=f'{deg}th polyfit, R²={r2:.4f}')

    ax.set_xlabel('Wavelength [nm]')
    ax.set_ylabel('Transmission [dB]')
    ax.set_title('Reference Spectrum and Polynomial Fitting',
                 fontsize=11, fontweight='bold')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc='lower center', ncol=2, framealpha=0.85)
    plt.tight_layout()
    plt.savefig(save_path, bbox_inches='tight')
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────
# Helpers for envelope (2-step flatten)
# ──────────────────────────────────────────────────────────────────────
def _envelope_from_peaks(bias_records, poly_deg=4):
    """모든 바이어스 spectrum 의 peak 들을 모아 글로벌 envelope 다항식 fit."""
    x_all, y_all = [], []
    wl_ref = None
    for V, (wl, y) in bias_records:
        pk = _find_peaks(y, window=14)
        if len(pk) == 0:
            continue
        x_all.append(wl[pk]); y_all.append(y[pk])
        if wl_ref is None:
            wl_ref = wl
    if not x_all:
        return None, None
    x_all = np.concatenate(x_all)
    y_all = np.concatenate(y_all)
    wl_mean = float(np.mean(wl_ref))
    coef = np.polyfit(x_all - wl_mean, y_all, poly_deg)

    def env(wl):
        return np.polyval(coef, wl - wl_mean)
    return env, wl_mean


# ──────────────────────────────────────────────────────────────────────
# Plot 03 — Flat transmission (2-step flatten)
# ──────────────────────────────────────────────────────────────────────
def plot_03_flat_transmission(die, save_path):
    sweeps = die['sweeps']
    biases = sorted(sweeps.keys())
    bias_records = [(V, sweeps[V]) for V in biases]
    env, wl_mean = _envelope_from_peaks(bias_records, poly_deg=4)
    if env is None:
        return

    cmap = get_cmap('coolwarm')
    fig, ax = plt.subplots(figsize=(8.5, 5), dpi=120)
    for i, (V, (wl, y)) in enumerate(bias_records):
        flat = y - env(wl)
        flat -= np.max(flat)   # max → 0 으로 정규화
        color = cmap(i / max(len(biases) - 1, 1))
        ax.plot(wl, flat, lw=0.9, color=color, label=f'MZM {V:+.1f} V')

    # reference 도 같이 (envelope 차감)
    if die['ref_L'] is not None:
        ref_flat = die['ref_IL'] - env(die['ref_L'])
        ref_flat -= np.max(ref_flat)
        ax.plot(die['ref_L'], ref_flat, color='magenta', lw=1.0, alpha=0.7,
                label='Reference (flat)')

    ax.set_xlabel('Wavelength [nm]')
    ax.set_ylabel('Flat transmission [dB]')
    ax.set_title('Flat transmission as measured (2-step flatten)',
                 fontsize=11, fontweight='bold')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc='lower center', ncol=4, framealpha=0.85)
    plt.tight_layout()
    plt.savefig(save_path, bbox_inches='tight')
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────
# Plot 04 — MZI cos² model fit (single bias)
# ──────────────────────────────────────────────────────────────────────
def plot_04_mzi_fit(die, save_path, target_bias=-1.0):
    sweeps = die['sweeps']
    biases = sorted(sweeps.keys())
    bias_records = [(V, sweeps[V]) for V in biases]

    # target_bias 에 가장 가까운 bias 선택
    V_sel = min(biases, key=lambda v: abs(v - target_bias))
    wl_sel, y_sel = sweeps[V_sel]

    env, wl_mean = _envelope_from_peaks(bias_records, poly_deg=4)
    if env is None:
        return

    # 2-step flatten
    flat_dB = y_sel - env(wl_sel)
    flat_dB -= np.max(flat_dB)
    T_lin = 10 ** (flat_dB / 10)
    T_lin = T_lin / np.max(T_lin)

    # FSR 초기값 — valley-to-valley 거리
    peak_idx = _find_peaks(flat_dB, window=14)
    if len(peak_idx) < 3:
        fsr_init = 10.0   # fallback
    else:
        mid = peak_idx[len(peak_idx) // 2]
        left = peak_idx[len(peak_idx) // 2 - 1]
        right = peak_idx[len(peak_idx) // 2 + 1]
        lv = left + int(np.argmin(flat_dB[left:mid]))
        rv = mid  + int(np.argmin(flat_dB[mid:right]))
        fsr_init = float(wl_sel[rv] - wl_sel[lv])
        if fsr_init <= 0:
            fsr_init = 10.0

    # Grid search (A, B 는 closed-form)
    wl0_guess = float(np.mean(wl_sel))
    fsr_grid = np.linspace(max(0.5, fsr_init * 0.7), fsr_init * 1.3, 21)
    phi_grid = np.linspace(0, np.pi, 18, endpoint=False)
    wl0_grid = wl0_guess + np.linspace(-fsr_init / 2, fsr_init / 2, 11)

    best = (np.inf, None)
    for f in fsr_grid:
        for p in phi_grid:
            for w in wl0_grid:
                A, B = _fit_AB(wl_sel, T_lin, f, p, w)
                r = T_lin - _mzi_model(wl_sel, A, B, f, p, w)
                ss = float(r @ r)
                if ss < best[0]:
                    best = (ss, (A, B, f, p, w))

    # Local refinement
    A, B, f, p, w = best[1]
    sf = (fsr_grid[1] - fsr_grid[0]) / 2
    sp = (phi_grid[1] - phi_grid[0]) / 2
    sw = (wl0_grid[1] - wl0_grid[0]) / 2
    for _ in range(40):
        improved = False
        for df in (-sf, 0, sf):
            for dp in (-sp, 0, sp):
                for dw in (-sw, 0, sw):
                    ff, pp, ww = f + df, p + dp, w + dw
                    Ai, Bi = _fit_AB(wl_sel, T_lin, ff, pp, ww)
                    r = T_lin - _mzi_model(wl_sel, Ai, Bi, ff, pp, ww)
                    ss = float(r @ r)
                    if ss < best[0] - 1e-15:
                        best = (ss, (Ai, Bi, ff, pp, ww))
                        A, B, f, p, w = Ai, Bi, ff, pp, ww
                        improved = True
        if not improved:
            sf *= 0.5; sp *= 0.5; sw *= 0.5
            if sf < 1e-7:
                break

    T_fit = _mzi_model(wl_sel, A, B, f, p, w)
    # 극솟값 보정
    T_fit -= (np.min(T_fit) - np.min(T_lin))
    r2 = _r2(T_lin, T_fit)

    fig, ax = plt.subplots(figsize=(8.5, 5), dpi=120)
    ax.plot(wl_sel, T_lin, 'b-', lw=1.0,
            label=f'Flat MZM raw ({V_sel:+.1f} V)')
    ax.plot(wl_sel, T_fit, 'k--', lw=2.0,
            label=f'MZI fit, R²={r2:.4f}')
    ax.set_xlabel('Wavelength [nm]')
    ax.set_ylabel('Normalized transmission')
    ax.set_title(f'MZM fitting after 2-step flatten ({V_sel:+.1f} V)\n'
                 f'A={A:.3f}, B={B:.3f}, FSR={f:.3f} nm, φ={p:.3f}, λ₀={w:.3f}',
                 fontsize=10, fontweight='bold')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9, framealpha=0.85)
    plt.tight_layout()
    plt.savefig(save_path, bbox_inches='tight')
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────
# Plot 05 — IV semilog
# ──────────────────────────────────────────────────────────────────────
def plot_05_iv_semilog(die, save_path):
    if die['iv_V'] is None:
        return
    V = die['iv_V']
    # 부호 컨벤션: forward I > 0 으로 통일 (HY202103 raw 는 음수로 출력)
    I_signed = -die['iv_I']
    I_abs = np.clip(np.abs(I_signed), 1e-15, None)

    fig, ax = plt.subplots(figsize=(8.5, 5), dpi=120)
    ax.semilogy(V, I_abs, 'o', ms=6, mfc='steelblue', mec='black',
                lw=0, label='Measured IV')
    ax.axvline(0, color='gray', lw=0.5)
    ax.set_xlabel('Voltage [V]')
    ax.set_ylabel('Current [A]')
    ax.set_title('IV analysis', fontsize=11, fontweight='bold')
    ax.grid(alpha=0.3, which='both')
    ax.legend(fontsize=9, framealpha=0.85)
    plt.tight_layout()
    plt.savefig(save_path, bbox_inches='tight')
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────
# Plot 06 — IV fitting (reverse poly + forward Shockley)
# ──────────────────────────────────────────────────────────────────────
def plot_06_iv_fit(die, save_path):
    if die['iv_V'] is None:
        return
    V = np.asarray(die['iv_V'], dtype=float)
    I_signed = -np.asarray(die['iv_I'], dtype=float)
    I_abs = np.clip(np.abs(I_signed), 1e-15, None)

    # 분할
    V_rev_mask = V <= 0.25
    V_fwd_mask = V >= 0.5

    V_rev, I_rev = V[V_rev_mask], I_abs[V_rev_mask]
    V_fwd, I_fwd = V[V_fwd_mask], I_abs[V_fwd_mask]

    fig, ax = plt.subplots(figsize=(8.5, 5), dpi=120)
    ax.semilogy(V, I_abs, 'o', ms=6, mfc='gray', mec='black',
                lw=0, label='Measured IV')

    info_lines = []
    # 역방향 polynomial fit (log10|I|)
    if len(V_rev) >= 4:
        deg_rev = 6
        logI_rev = np.log10(I_rev)
        coef_rev = np.polyfit(V_rev, logI_rev, deg_rev)
        V_rev_dense = np.linspace(V_rev.min(), V_rev.max(), 200)
        I_rev_fit = 10 ** np.polyval(coef_rev, V_rev_dense)
        r2_rev = _r2(logI_rev, np.polyval(coef_rev, V_rev))
        ax.plot(V_rev_dense, I_rev_fit, '-', color='orange', lw=2,
                label=f'Reverse polynomial fit (R²={r2_rev:.3f})')
        info_lines.append(f'R²_rev = {r2_rev:.4f}')

    # 순방향 Shockley fit
    if len(V_fwd) >= 3:
        try:
            popt, _ = curve_fit(_shockley, V_fwd, I_fwd,
                                p0=[1e-12, 1.5],
                                bounds=([1e-20, 0.5], [1e-3, 5]),
                                maxfev=5000)
            Is, n_ideal = popt
            V_fwd_dense = np.linspace(V_fwd.min(), V_fwd.max(), 200)
            I_fwd_fit = _shockley(V_fwd_dense, *popt)
            r2_fwd = _r2(I_fwd, _shockley(V_fwd, *popt))
            ax.plot(V_fwd_dense, I_fwd_fit, '-', color='green', lw=2,
                    label=f'Forward diode fit (n={n_ideal:.2f})')
            info_lines.append(f'Is = {Is:.2e} A')
            info_lines.append(f'n  = {n_ideal:.3f}')
            info_lines.append(f'R²_fwd = {r2_fwd:.4f}')
        except Exception as e:
            info_lines.append(f'Shockley fit fail: {e}')

    if info_lines:
        ax.text(0.02, 0.97, '\n'.join(info_lines), transform=ax.transAxes,
                va='top', fontsize=9, family='monospace',
                bbox=dict(facecolor='white', alpha=0.85))

    ax.axvline(0, color='gray', lw=0.5)
    ax.set_xlabel('Voltage [V]')
    ax.set_ylabel('Current [A]')
    ax.set_title('IV analysis — reverse polynomial + forward diode fit',
                 fontsize=11, fontweight='bold')
    ax.grid(alpha=0.3, which='both')
    ax.legend(fontsize=9, loc='lower right', framealpha=0.85)
    plt.tight_layout()
    plt.savefig(save_path, bbox_inches='tight')
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────
# Per-die worker
# ──────────────────────────────────────────────────────────────────────
def _process_one(args):
    fp, date, wafer, row, col, band = args
    die = load_die(fp)
    if die is None:
        return f'{fp}: parse fail'

    out_dir = os.path.join(OUT_ROOT, date, f'{band}-band', wafer,
                            f'R{row:+d}_C{col:+d}')
    os.makedirs(out_dir, exist_ok=True)

    try:
        plot_01_mzm_ref_spectra(die, os.path.join(out_dir, '01_mzm_ref_spectra.png'))
        plot_02_ref_polyfit   (die, os.path.join(out_dir, '02_ref_polyfit.png'))
        plot_03_flat_transmission(die, os.path.join(out_dir, '03_flat_transmission.png'))
        plot_04_mzi_fit       (die, os.path.join(out_dir, '04_mzi_fit.png'))
        plot_05_iv_semilog    (die, os.path.join(out_dir, '05_iv_semilog.png'))
        plot_06_iv_fit        (die, os.path.join(out_dir, '06_iv_fit.png'))
        return f'OK  {date}/{band}-band/{wafer}/R{row:+d}_C{col:+d}'
    except Exception as e:
        return f'FAIL {date}/{band}-band/{wafer}/R{row:+d}_C{col:+d}: {e}'


def main():
    print('=' * 70)
    print(' 다이별 6 종류 분석 그림 생성 (PE2_Lec06 Objective 6)')
    print('=' * 70)
    items = _find_all_xmls_with_meta(DATA_ROOT)
    print(f'\n발견된 측정: {len(items)}개')
    print(f'출력 위치: {OUT_ROOT}\n')

    with multiprocessing.Pool() as pool:
        results = pool.map(_process_one, items)

    ok = sum(1 for r in results if r.startswith('OK'))
    fail = len(results) - ok
    print(f'\n완료: {ok}개 OK, {fail}개 실패')
    if fail > 0:
        print('실패 목록 (앞 10개):')
        for r in results:
            if not r.startswith('OK'):
                print('  ' + r)
                fail -= 1
                if fail <= 0: break


if __name__ == '__main__':
    main()
