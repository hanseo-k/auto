"""다이별 6 종류 분석 그림 — 공통 helper + plot 함수.

band 별 파라미터는 호출자가 dict 로 전달.  실제 main entry 는
    src/plot_per_die_figures_C.py  (C-band 전용)
    src/plot_per_die_figures_O.py  (O-band 전용)

config 구조:
    {
      'band':                 'C' or 'O',
      'er_window':            (lo_nm, hi_nm),
      'envelope_poly_deg':    int,
      'envelope_peak_window': int,
      'fsr_expected':         float (nm),
      'trim_frac':            float (양 끝단 trim 비율),
      'ref_poly_deg':         int (03 ref-fit 차수),
      'mzi_target_bias':      float (04 fit 대상 bias),
    }
"""
import sys, os, re, glob
PROGRAM_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROGRAM_ROOT, 'src'))

import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from scipy.signal import find_peaks as _scipy_find_peaks

from xml_loader import load_die, BAND_OF_FILE


DATA_ROOT = os.path.join(PROGRAM_ROOT, 'data', 'HY202103')
OUT_ROOT  = os.path.join(PROGRAM_ROOT, 'res', 'figures_per_die')


# ──────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────
def r2(y, y_fit):
    ss_res = float(np.sum((y - y_fit) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')


def find_peaks_local(y, window=14):
    y = np.asarray(y); n = len(y); peaks = []
    for i in range(window, n - window):
        if y[i] == np.max(y[i - window: i + window + 1]):
            peaks.append(i)
    return np.array(peaks, dtype=int)


def find_peaks_global(y, window=800):
    y = np.asarray(y); n = len(y); peaks = []
    for i in range(0, n):
        lo = max(0, i - window)
        hi = min(n, i + window + 1)
        if y[i] == np.max(y[lo:hi]):
            peaks.append(i)
    return np.array(peaks, dtype=int)


def mzi_model(wl, A, B, fsr, phi, wl0):
    return A + B * np.cos(np.pi * (wl - wl0) / fsr + phi) ** 2


def fit_AB(wl, T, fsr, phi, wl0):
    c2 = np.cos(np.pi * (wl - wl0) / fsr + phi) ** 2
    M = np.column_stack([np.ones_like(wl), c2])
    coef, *_ = np.linalg.lstsq(M, T, rcond=None)
    return float(coef[0]), float(coef[1])


def envelope_from_peaks(bias_records, poly_deg, peak_window):
    """모든 bias 의 큰 global peak 을 모아 다항식 envelope fit."""
    x_all, y_all = [], []
    wl_ref = None
    for V, (wl, y) in bias_records:
        pk = find_peaks_global(y, window=peak_window)
        if len(pk) == 0:
            continue
        x_all.append(wl[pk]); y_all.append(y[pk])
        if wl_ref is None:
            wl_ref = wl
    if not x_all:
        return None, None
    x_all = np.concatenate(x_all); y_all = np.concatenate(y_all)
    wl_mean = float(np.mean(wl_ref))
    coef = np.polyfit(x_all - wl_mean, y_all, poly_deg)
    def env(wl): return np.polyval(coef, wl - wl_mean)
    return env, wl_mean


def find_all_xmls_with_meta(data_root):
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
            items.append((fp, date_str, m_rc.group(1),
                          int(m_rc.group(2)), int(m_rc.group(3)), band))
    return items


# ──────────────────────────────────────────────────────────────────────
# Plot 01 — MZM raw spectra + ALIGN reference
# ──────────────────────────────────────────────────────────────────────
def plot_01_mzm_ref_spectra(die, save_path, config):
    sweeps = die['sweeps']; biases = sorted(sweeps.keys())
    cmap = get_cmap('coolwarm')
    fig, ax = plt.subplots(figsize=(8.5, 5), dpi=120)
    for i, V in enumerate(biases):
        L, IL = sweeps[V]
        ax.plot(L, IL, lw=1.0, color=cmap(i / max(len(biases) - 1, 1)),
                label=f'MZM {V:+.1f} V')
    if die['ref_L'] is not None:
        ax.plot(die['ref_L'], die['ref_IL'], color='black', lw=2.0,
                label='Reference')
    ax.set_xlabel('Wavelength [nm]'); ax.set_ylabel('Measured transmission [dB]')
    ax.set_title(f'MZM and Reference spectra ({config["band"]}-band)',
                 fontsize=11, fontweight='bold')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc='lower center', ncol=4, framealpha=0.85)
    plt.tight_layout(); plt.savefig(save_path, bbox_inches='tight'); plt.close(fig)


# ──────────────────────────────────────────────────────────────────────
# Plot 02 — ALIGN reference polynomial fit (2~6 차)
# ──────────────────────────────────────────────────────────────────────
def plot_02_ref_polyfit(die, save_path, config):
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
        r2v = r2(ref_IL, fit)
        ax.plot(ref_L, fit, lw=1.5, color=cmap(i / (len(degrees) - 1)),
                label=f'{deg}th polyfit, R²={r2v:.4f}')
    ax.set_xlabel('Wavelength [nm]'); ax.set_ylabel('Transmission [dB]')
    ax.set_title(f'Reference Spectrum and Polynomial Fitting ({config["band"]}-band)',
                 fontsize=11, fontweight='bold')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc='lower center', ncol=2, framealpha=0.85)
    plt.tight_layout(); plt.savefig(save_path, bbox_inches='tight'); plt.close(fig)


# ──────────────────────────────────────────────────────────────────────
# Plot 03 — Flat transmission (ref-fit subtraction)
# ──────────────────────────────────────────────────────────────────────
def plot_03_flat_transmission(die, save_path, config):
    if die['ref_L'] is None:
        return
    sweeps = die['sweeps']; biases = sorted(sweeps.keys())
    ref_L, ref_IL = die['ref_L'], die['ref_IL']
    wl_mean = float(np.mean(ref_L))
    ref_coef = np.polyfit(ref_L - wl_mean, ref_IL, config['ref_poly_deg'])
    def ref_fit_at(wl): return np.polyval(ref_coef, wl - wl_mean)

    cmap = get_cmap('coolwarm')
    fig, ax = plt.subplots(figsize=(8.5, 5), dpi=120)
    for i, V in enumerate(biases):
        wl, y = sweeps[V]
        flat = y - ref_fit_at(wl); flat -= np.max(flat)
        ax.plot(wl, flat, lw=1.0,
                color=cmap(i / max(len(biases) - 1, 1)),
                label=f'MZM {V:+.1f} V')
    ref_flat = ref_IL - ref_fit_at(ref_L); ref_flat -= np.max(ref_flat)
    ax.plot(ref_L, ref_flat, color='magenta', lw=1.4,
            label=f'Reference (flat after {config["ref_poly_deg"]}-th polyfit)')
    ax.set_xlabel('Wavelength [nm]'); ax.set_ylabel('Flat transmission [dB]')
    ax.set_title(f'Flat transmission — reference polynomial subtracted ({config["band"]}-band)',
                 fontsize=11, fontweight='bold')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc='lower center', ncol=4, framealpha=0.85)
    plt.tight_layout(); plt.savefig(save_path, bbox_inches='tight'); plt.close(fig)


# ──────────────────────────────────────────────────────────────────────
# Plot 04 — MZI cos² model fit
# ──────────────────────────────────────────────────────────────────────
def plot_04_mzi_fit(die, save_path, config):
    sweeps = die['sweeps']; biases = sorted(sweeps.keys())
    bias_records = [(V, sweeps[V]) for V in biases]

    target_bias = config['mzi_target_bias']
    V_sel = min(biases, key=lambda v: abs(v - target_bias))
    wl_full, y_full = sweeps[V_sel]

    env, wl_mean = envelope_from_peaks(
        bias_records,
        poly_deg=config['envelope_poly_deg'],
        peak_window=config['envelope_peak_window'])
    if env is None:
        return

    # 양 끝단 trim
    n_pts = len(wl_full)
    trim = max(50, int(n_pts * config['trim_frac']))
    wl_sel = wl_full[trim:-trim]
    y_sel  = y_full[trim:-trim]

    flat_dB = y_sel - env(wl_sel); flat_dB -= np.max(flat_dB)
    T_lin = 10 ** (flat_dB / 10); T_lin = T_lin / np.max(T_lin)

    # FSR 초기값 — global peak 으로
    peak_idx = find_peaks_global(flat_dB, window=config['envelope_peak_window'])
    if len(peak_idx) < 3:
        peak_idx = find_peaks_local(flat_dB, window=14)
    if len(peak_idx) < 3:
        fsr_init = config['fsr_expected']
    else:
        mid = peak_idx[len(peak_idx) // 2]
        left = peak_idx[len(peak_idx) // 2 - 1]
        right = peak_idx[len(peak_idx) // 2 + 1]
        lv = left + int(np.argmin(flat_dB[left:mid]))
        rv = mid  + int(np.argmin(flat_dB[mid:right]))
        fsr_init = float(wl_sel[rv] - wl_sel[lv])
        if fsr_init <= 0:
            fsr_init = config['fsr_expected']

    wl0_guess = float(np.mean(wl_sel))
    fsr_grid = np.linspace(max(0.5, fsr_init * 0.7), fsr_init * 1.3, 41)
    phi_grid = np.linspace(0, np.pi, 36, endpoint=False)
    wl0_grid = wl0_guess + np.linspace(-fsr_init / 2, fsr_init / 2, 21)

    best = (np.inf, None)
    for f in fsr_grid:
        for p in phi_grid:
            for w in wl0_grid:
                A, B = fit_AB(wl_sel, T_lin, f, p, w)
                rr = T_lin - mzi_model(wl_sel, A, B, f, p, w)
                ss = float(rr @ rr)
                if ss < best[0]:
                    best = (ss, (A, B, f, p, w))
    A, B, f, p, w = best[1]
    sf = (fsr_grid[1] - fsr_grid[0]) / 2
    sp = (phi_grid[1] - phi_grid[0]) / 2
    sw = (wl0_grid[1] - wl0_grid[0]) / 2
    for _ in range(80):
        improved = False
        for df in (-sf, 0, sf):
            for dp in (-sp, 0, sp):
                for dw in (-sw, 0, sw):
                    ff, pp, ww = f + df, p + dp, w + dw
                    Ai, Bi = fit_AB(wl_sel, T_lin, ff, pp, ww)
                    rr = T_lin - mzi_model(wl_sel, Ai, Bi, ff, pp, ww)
                    ss = float(rr @ rr)
                    if ss < best[0] - 1e-15:
                        best = (ss, (Ai, Bi, ff, pp, ww))
                        A, B, f, p, w = Ai, Bi, ff, pp, ww
                        improved = True
        if not improved:
            sf *= 0.5; sp *= 0.5; sw *= 0.5
            if sf < 1e-7:
                break

    T_fit = mzi_model(wl_sel, A, B, f, p, w)
    r2v = r2(T_lin, T_fit)
    fit_span = float(T_fit.max() - T_fit.min())
    T_fit_plot = ((T_fit - T_fit.min()) / fit_span) if fit_span > 1e-12 else T_fit

    fig, ax = plt.subplots(figsize=(8.5, 5), dpi=120)
    ax.plot(wl_sel, T_lin, 'b-', lw=1.0,
            label=f'Flat MZM raw ({V_sel:+.1f} V)')
    ax.plot(wl_sel, T_fit_plot, 'k--', lw=2.0,
            label=f'MZI fit (normalized), R²={r2v:.4f}')
    ax.set_xlabel('Wavelength [nm]'); ax.set_ylabel('Normalized transmission')
    ax.set_title(f'MZM fitting after flatten ({V_sel:+.1f} V, {config["band"]}-band)\n'
                 f'A={A:.3f}, B={B:.3f}, FSR={f:.3f} nm, φ={p:.3f}, λ₀={w:.3f}',
                 fontsize=10, fontweight='bold')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9, framealpha=0.85)
    plt.tight_layout(); plt.savefig(save_path, bbox_inches='tight'); plt.close(fig)


# ──────────────────────────────────────────────────────────────────────
# Plot 05 — IV semilog
# ──────────────────────────────────────────────────────────────────────
def plot_05_iv_semilog(die, save_path, config):
    if die['iv_V'] is None:
        return
    V = die['iv_V']
    I_signed = -die['iv_I']
    I_abs = np.clip(np.abs(I_signed), 1e-15, None)
    fig, ax = plt.subplots(figsize=(8.5, 5), dpi=120)
    ax.semilogy(V, I_abs, 'o', ms=6, mfc='steelblue', mec='black',
                lw=0, label='Measured IV')
    ax.axvline(0, color='gray', lw=0.5)
    ax.set_xlabel('Voltage [V]'); ax.set_ylabel('Current [A]')
    ax.set_title(f'IV analysis ({config["band"]}-band)',
                 fontsize=11, fontweight='bold')
    ax.grid(alpha=0.3, which='both')
    ax.legend(fontsize=9, framealpha=0.85)
    plt.tight_layout(); plt.savefig(save_path, bbox_inches='tight'); plt.close(fig)


# ──────────────────────────────────────────────────────────────────────
# Plot 06 — IV fit (V<=0.25 3차, V>=0.5 1차)
# ──────────────────────────────────────────────────────────────────────
def plot_06_iv_fit(die, save_path, config):
    if die['iv_V'] is None:
        return
    V = np.asarray(die['iv_V'], dtype=float)
    I_signed = -np.asarray(die['iv_I'], dtype=float)
    I_abs = np.clip(np.abs(I_signed), 1e-15, None)
    logI = np.log10(I_abs)

    V_low_mask  = V <= 0.25
    V_high_mask = V >= 0.5
    V_low,  logI_low  = V[V_low_mask],  logI[V_low_mask]
    V_high, logI_high = V[V_high_mask], logI[V_high_mask]

    fig, ax = plt.subplots(figsize=(8.5, 5), dpi=120)
    ax.semilogy(V, I_abs, 'o', ms=6, mfc='gray', mec='black',
                lw=0, label='Measured IV')

    info = []
    if len(V_low) >= 4:
        coef_low = np.polyfit(V_low, logI_low, 3)
        V_low_d = np.linspace(V_low.min(), V_low.max(), 200)
        I_low_fit = 10 ** np.polyval(coef_low, V_low_d)
        r2l = r2(logI_low, np.polyval(coef_low, V_low))
        ax.plot(V_low_d, I_low_fit, '-', color='orange', lw=2,
                label=f'V∈[-2.0, 0.25]: 3rd polyfit (R²={r2l:.3f})')
        info.append(f'V≤0.25:  3rd polyfit, R² = {r2l:.4f}')

    if len(V_high) >= 2:
        coef_high = np.polyfit(V_high, logI_high, 1)
        V_high_d = np.linspace(V_high.min(), V_high.max(), 200)
        I_high_fit = 10 ** np.polyval(coef_high, V_high_d)
        slope = coef_high[0]
        VT = 0.02585
        ideal = 1.0 / (slope * np.log(10) * VT) if slope > 0 else float('nan')
        ax.plot(V_high_d, I_high_fit, '-', color='green', lw=2,
                label=f'V∈[0.5, 1.0]: 1st polyfit (slope={slope:.1f} dec/V)')
        info.append(f'V≥0.5:   1st polyfit, slope = {slope:.2f} dec/V')
        info.append(f'         ideality n = {ideal:.2f}')

    if info:
        ax.text(0.02, 0.97, '\n'.join(info), transform=ax.transAxes,
                va='top', fontsize=9, family='monospace',
                bbox=dict(facecolor='white', alpha=0.85))
    ax.axvline(0, color='gray', lw=0.5)
    ax.axvspan(0.25, 0.5, color='lightgray', alpha=0.3,
               label='gap (no measurement)')
    ax.set_xlabel('Voltage [V]'); ax.set_ylabel('Current [A]')
    ax.set_title(f'IV analysis — V≤0.25: 3rd polyfit / V≥0.5: 1st polyfit ({config["band"]}-band)',
                 fontsize=11, fontweight='bold')
    ax.grid(alpha=0.3, which='both')
    ax.legend(fontsize=9, loc='lower right', framealpha=0.85)
    plt.tight_layout(); plt.savefig(save_path, bbox_inches='tight'); plt.close(fig)


# ──────────────────────────────────────────────────────────────────────
# Helper: null finding near ER window center
# ──────────────────────────────────────────────────────────────────────
def _parabolic_min(L, IL, idx, half_pts=4):
    """idx 주변 ±half_pts 에서 parabolic fit 으로 정밀 null wavelength."""
    lo, hi = max(0, idx - half_pts), min(len(L), idx + half_pts + 1)
    if hi - lo < 3:
        return float(L[idx])
    Lw, Iw = L[lo:hi], IL[lo:hi]
    a, b, _ = np.polyfit(Lw, Iw, 2)
    if a <= 0:
        return float(L[idx])
    return float(-b / (2 * a))


def _scipy_null_indices(IL_w, prominence=10.0, distance=20):
    """IL 의 깊은 null index — prominence 와 distance 로 작은 ripple 제외.

    prominence: 깊이 기준 (dB). 10 이면 인접 baseline 보다 ≥ 10 dB 깊은 null만.
    distance:   인접 null 사이 최소 거리 (sample 수). 노이즈 ripple 클러스터 제거.
    """
    idx, _ = _scipy_find_peaks(-IL_w, prominence=prominence, distance=distance)
    return idx.astype(int)


def find_two_nulls_in_window(wl, IL, win_lo, win_hi,
                              prominence=10.0, distance=20):
    """ER window 안에서 prominent 한 두 인접 null 의 wavelength.

    반환: (deep_wl, neighbor_wl) 또는 None.  파장 순서로 정렬 안 됨.
    """
    mask = (wl >= win_lo) & (wl <= win_hi)
    wl_w, IL_w = wl[mask], IL[mask]
    if len(wl_w) < 50:
        return None
    null_idx = _scipy_null_indices(IL_w, prominence=prominence, distance=distance)
    if len(null_idx) < 2:
        return None
    deepest = int(null_idx[np.argmin(IL_w[null_idx])])
    deepest_wl = _parabolic_min(wl_w, IL_w, deepest)
    sorted_idx = sorted(null_idx.tolist())
    pos = sorted_idx.index(deepest)
    if pos + 1 < len(sorted_idx):
        neighbor_idx = sorted_idx[pos + 1]
    elif pos - 1 >= 0:
        neighbor_idx = sorted_idx[pos - 1]
    else:
        return None
    neighbor_wl = _parabolic_min(wl_w, IL_w, neighbor_idx)
    return (deepest_wl, neighbor_wl)


def find_deepest_null_in_window(wl, IL, win_lo, win_hi,
                                 prominence=10.0, distance=20):
    """ER window 안 가장 prominent 한 (= 가장 깊은) null 의 wavelength + IL."""
    mask = (wl >= win_lo) & (wl <= win_hi)
    wl_w, IL_w = wl[mask], IL[mask]
    if len(wl_w) < 50:
        return None
    null_idx = _scipy_null_indices(IL_w, prominence=prominence, distance=distance)
    if len(null_idx) == 0:
        # fallback — 그냥 단순 최소
        idx = int(np.argmin(IL_w))
        return _parabolic_min(wl_w, IL_w, idx), float(IL_w[idx])
    deepest = int(null_idx[np.argmin(IL_w[null_idx])])
    return _parabolic_min(wl_w, IL_w, deepest), float(IL_w[deepest])


# ──────────────────────────────────────────────────────────────────────
# Plot 07 — FSR per bias (인접 두 null 거리)
# ──────────────────────────────────────────────────────────────────────
def plot_07_fsr_per_bias(die, save_path, config):
    """각 bias 마다 ER window 안에서 인접 두 깊은 null 의 거리 = FSR.
    세로 점선으로 표시하고, box 에 bias 별 FSR 표.
    """
    sweeps = die['sweeps']; biases = sorted(sweeps.keys())
    win_lo, win_hi = config['er_window']
    # 윈도우를 살짝 넓혀서 인접 null 도 포함 (FSR 만큼)
    plot_lo = win_lo - config['fsr_expected'] * 0.5
    plot_hi = win_hi + config['fsr_expected'] * 0.5

    cmap = get_cmap('coolwarm')
    fig, ax = plt.subplots(figsize=(10, 6), dpi=120)

    fsr_table = []  # (V, fsr, deep_wl, neighbor_wl)
    for i, V in enumerate(biases):
        wl, IL = sweeps[V]
        mask = (wl >= plot_lo) & (wl <= plot_hi)
        color = cmap(i / max(len(biases) - 1, 1))
        ax.plot(wl[mask], IL[mask], lw=0.9, color=color,
                label=f'{V:+.1f} V')
        nulls = find_two_nulls_in_window(wl, IL, plot_lo, plot_hi)
        if nulls is None:
            continue
        a, b = sorted(nulls)
        fsr_val = abs(b - a)
        fsr_table.append((V, fsr_val, a, b))
        # 세로 점선
        ax.axvline(a, color=color, ls='--', alpha=0.55, lw=0.9)
        ax.axvline(b, color=color, ls='--', alpha=0.55, lw=0.9)

    # FSR 표 (box)
    if fsr_table:
        rows = [f'{V:+.1f} V:  FSR = {f:.3f} nm   (nulls {a:.2f}, {b:.2f})'
                for V, f, a, b in fsr_table]
        ax.text(0.02, 0.97, 'FSR per bias\n' + '\n'.join(rows),
                transform=ax.transAxes, va='top',
                fontsize=8, family='monospace',
                bbox=dict(facecolor='white', alpha=0.9))

    ax.set_xlabel('Wavelength [nm]')
    ax.set_ylabel('Measured transmission [dB]')
    ax.set_title(f'FSR per bias — {config["band"]}-band '
                 f'(nulls in [{plot_lo:.0f}, {plot_hi:.0f}] nm)',
                 fontsize=11, fontweight='bold')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc='lower right', ncol=2, framealpha=0.85)
    plt.tight_layout(); plt.savefig(save_path, bbox_inches='tight'); plt.close(fig)


# ──────────────────────────────────────────────────────────────────────
# Plot 08 — Null zoom (가장 깊은 null 영역 확대, 전압별)
# ──────────────────────────────────────────────────────────────────────
def plot_08_null_zoom(die, save_path, config, half_window=2.0):
    """ER window 중심 부근 가장 깊은 null 의 ±half_window nm 확대.
    전압이 클수록 null 이 더 왼쪽으로 이동하는 것이 한눈에 보인다.
    """
    sweeps = die['sweeps']; biases = sorted(sweeps.keys())
    win_lo, win_hi = config['er_window']
    # baseline (V=0) 의 가장 깊은 null 위치 기준으로 zoom 영역 정함
    if 0.0 in sweeps:
        wl0, IL0 = sweeps[0.0]
    else:
        V0 = min(biases, key=lambda v: abs(v))
        wl0, IL0 = sweeps[V0]
    res = find_deepest_null_in_window(wl0, IL0, win_lo, win_hi)
    if res is None:
        return
    center_wl, _ = res
    zoom_lo = center_wl - half_window
    zoom_hi = center_wl + half_window

    cmap = get_cmap('coolwarm')
    fig, ax = plt.subplots(figsize=(9, 6), dpi=120)

    null_labels = []  # (V, wl) — 텍스트 표 만들기용
    for i, V in enumerate(biases):
        wl, IL = sweeps[V]
        mask = (wl >= zoom_lo) & (wl <= zoom_hi)
        color = cmap(i / max(len(biases) - 1, 1))
        ax.plot(wl[mask], IL[mask], lw=1.4, color=color,
                label=f'{V:+.1f} V')
        # 각 bias 의 null 위치
        res_V = find_deepest_null_in_window(wl, IL, zoom_lo, zoom_hi)
        if res_V is not None:
            null_wl, null_IL = res_V
            ax.axvline(null_wl, color=color, ls=':', lw=0.8, alpha=0.7)
            # null 위치에 작은 marker (아래쪽)
            ax.plot(null_wl, null_IL, 'v', color=color, ms=8,
                    mec='black', mew=0.5, zorder=5)
            null_labels.append((V, null_wl))

    # V=0 reference 세로선
    ax.axvline(center_wl, color='black', ls='-', lw=0.5, alpha=0.4)

    # 좌상단 박스에 각 bias 의 null wavelength 표
    if null_labels:
        lines = ['Null wavelength per bias:']
        for V, w in sorted(null_labels):
            lines.append(f'  {V:+.1f} V :  {w:.3f} nm')
        ax.text(0.02, 0.97, '\n'.join(lines), transform=ax.transAxes,
                va='top', fontsize=8, family='monospace',
                bbox=dict(facecolor='white', alpha=0.9, edgecolor='gray'))

    ax.set_xlabel('Wavelength [nm]')
    ax.set_ylabel('Measured transmission [dB]')
    ax.set_title(f'Focus on spectral null near {center_wl:.2f} nm '
                 f'({config["band"]}-band, ±{half_window:.1f} nm)',
                 fontsize=11, fontweight='bold')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc='lower center', ncol=4, framealpha=0.85)
    plt.tight_layout(); plt.savefig(save_path, bbox_inches='tight'); plt.close(fig)


# ──────────────────────────────────────────────────────────────────────
# Plot 09 — Phase shift vs bias
# ──────────────────────────────────────────────────────────────────────
def plot_09_phase_shift(die, save_path, config):
    """전압 별 null wavelength 이동량 (pm) 과 phase shift (rad).
    baseline = V=0 (없으면 가장 0에 가까운 bias).
    """
    sweeps = die['sweeps']; biases = sorted(sweeps.keys())
    win_lo, win_hi = config['er_window']
    if 0.0 in sweeps:
        V0 = 0.0
    else:
        V0 = min(biases, key=lambda v: abs(v))
    wl0, IL0 = sweeps[V0]
    res0 = find_deepest_null_in_window(wl0, IL0, win_lo, win_hi)
    if res0 is None:
        return
    null0_wl = res0[0]

    points = []  # (V, dlam_pm, dphase_rad)
    fsr = config['fsr_expected']
    for V in biases:
        wl, IL = sweeps[V]
        # 같은 null 추적 — V=0 null 위치 ± half_FSR 영역 안에서
        half = min(0.4, fsr * 0.35)
        sub_lo = null0_wl - half
        sub_hi = null0_wl + half
        res = find_deepest_null_in_window(wl, IL, sub_lo, sub_hi)
        if res is None:
            continue
        dlam_nm = res[0] - null0_wl
        dlam_pm = dlam_nm * 1000.0
        dphase = 2 * np.pi * dlam_nm / fsr   # phase shift (rad)
        points.append((V, dlam_pm, dphase))

    if not points:
        return
    V_arr   = np.array([p[0] for p in points])
    dlam_arr = np.array([p[1] for p in points])
    dphi_arr = np.array([p[2] for p in points])

    fig, ax1 = plt.subplots(figsize=(9, 5.5), dpi=120)
    ax1.plot(V_arr, dlam_arr, 'o-', color='steelblue', lw=1.5, ms=6,
             label='Δλ_null (pm)')
    ax1.axhline(0, color='gray', lw=0.5)
    ax1.axvline(0, color='gray', lw=0.5)
    ax1.set_xlabel('Bias V [V]')
    ax1.set_ylabel('Δλ_null vs V=0 [pm]', color='steelblue')
    ax1.tick_params(axis='y', labelcolor='steelblue')
    ax1.grid(alpha=0.3)

    ax2 = ax1.twinx()
    ax2.plot(V_arr, dphi_arr, 's--', color='crimson', lw=1.2, ms=5,
             label='Δφ (rad)')
    ax2.axhline(np.pi, color='crimson', ls=':', lw=0.7, alpha=0.5)
    ax2.text(V_arr.max(), np.pi, ' π', color='crimson',
             va='center', fontsize=9)
    ax2.set_ylabel('Δφ [rad]  (= 2π · Δλ / FSR)', color='crimson')
    ax2.tick_params(axis='y', labelcolor='crimson')

    ax1.set_title(f'Phase shift per bias ({config["band"]}-band)\n'
                  f'null tracked from {null0_wl:.3f} nm at V={V0:+.1f} V; '
                  f'FSR = {fsr:.2f} nm',
                  fontsize=11, fontweight='bold')
    plt.tight_layout(); plt.savefig(save_path, bbox_inches='tight'); plt.close(fig)


# ──────────────────────────────────────────────────────────────────────
# Plot 10 — V_pi·L vs bias
# ──────────────────────────────────────────────────────────────────────
def plot_10_vpiL_per_bias(die, save_path, config):
    """각 bias 별 V_π·L 추정값.

    Phase shift 정의:  Δφ(V) = 2π · (λ_null(V) − λ_null(0)) / FSR
    각 bias 에서 V → π 위상 변화까지 외삽:
        V_π(V_i) = V_i · π / Δφ(V_i)
    V_π·L = V_π · L  (length 는 die['length_um'] → cm 환산)

    Δφ → 0 (= V == V0) 인 점만 제외.  나머지 모든 bias 표시.
    작은 |Δφ| 는 자연스럽게 큰 V_π 값을 주지만 그것도 디바이스의 비선형
    특성이므로 정보로 보존.
    """
    sweeps = die['sweeps']; biases = sorted(sweeps.keys())
    win_lo, win_hi = config['er_window']
    if 0.0 in sweeps:
        V0 = 0.0
    else:
        V0 = min(biases, key=lambda v: abs(v))
    res0 = find_deepest_null_in_window(sweeps[V0][0], sweeps[V0][1],
                                       win_lo, win_hi)
    if res0 is None:
        return
    null0_wl = res0[0]

    L_um = die.get('length_um') or 500
    L_cm = L_um * 1e-4

    fsr = config['fsr_expected']
    points = []
    for V in biases:
        if abs(V - V0) < 1e-9:
            continue
        wl, IL = sweeps[V]
        half = min(0.4, fsr * 0.35)
        res = find_deepest_null_in_window(wl, IL, null0_wl - half, null0_wl + half)
        if res is None:
            continue
        dlam_nm = res[0] - null0_wl
        dphi = 2 * np.pi * dlam_nm / fsr
        if abs(dphi) < 1e-6:   # 사실상 V == V0 만 제외 (분모 0 방지)
            continue
        vpi_at_V = (V - V0) * np.pi / dphi
        vpiL = vpi_at_V * L_cm
        points.append((V, vpiL, vpi_at_V, dphi))

    if not points:
        return
    V_arr    = np.array([p[0] for p in points])
    vpiL_arr = np.array([p[1] for p in points])
    vpi_arr  = np.array([p[2] for p in points])

    fig, ax = plt.subplots(figsize=(9, 5.5), dpi=120)
    rev_mask = V_arr <= 0
    fwd_mask = V_arr > 0
    ax.plot(V_arr[rev_mask], np.abs(vpiL_arr[rev_mask]), 'o-',
            color='steelblue', lw=1.5, ms=6, label='reverse')
    ax.plot(V_arr[fwd_mask], np.abs(vpiL_arr[fwd_mask]), 's-',
            color='crimson', lw=1.5, ms=6, label='forward')
    ax.axhline(np.median(np.abs(vpiL_arr)), color='gray', ls=':', lw=0.8,
               label=f'median = {np.median(np.abs(vpiL_arr)):.3f} V·cm')
    ax.set_xlabel('Bias V [V]')
    ax.set_ylabel('|V_π · L| [V·cm]')
    ax.set_title(f'V_π·L per bias ({config["band"]}-band, L = {L_um} μm)\n'
                 f'extrapolated from null shift at each V using FSR = {fsr:.2f} nm',
                 fontsize=11, fontweight='bold')
    ax.grid(alpha=0.3, which='both')
    ax.legend(fontsize=9, framealpha=0.85)
    plt.tight_layout(); plt.savefig(save_path, bbox_inches='tight'); plt.close(fig)


# ──────────────────────────────────────────────────────────────────────
# Per-die worker (config-driven)
# ──────────────────────────────────────────────────────────────────────
def process_one(args, config):
    fp, date, wafer, row, col, band = args
    die = load_die(fp)
    if die is None:
        return f'{fp}: parse fail'
    out_dir = os.path.join(OUT_ROOT, date, f'{band}-band', wafer,
                            f'({row},{col})')
    os.makedirs(out_dir, exist_ok=True)
    try:
        plot_01_mzm_ref_spectra  (die, os.path.join(out_dir, '01_mzm_ref_spectra.png'), config)
        plot_02_ref_polyfit      (die, os.path.join(out_dir, '02_ref_polyfit.png'), config)
        plot_03_flat_transmission(die, os.path.join(out_dir, '03_flat_transmission.png'), config)
        plot_04_mzi_fit          (die, os.path.join(out_dir, '04_mzi_fit.png'), config)
        plot_05_iv_semilog       (die, os.path.join(out_dir, '05_iv_semilog.png'), config)
        plot_06_iv_fit           (die, os.path.join(out_dir, '06_iv_fit.png'), config)
        plot_07_fsr_per_bias     (die, os.path.join(out_dir, '07_fsr_per_bias.png'), config)
        plot_08_null_zoom        (die, os.path.join(out_dir, '08_null_zoom.png'), config)
        plot_09_phase_shift      (die, os.path.join(out_dir, '09_phase_shift.png'), config)
        plot_10_vpiL_per_bias    (die, os.path.join(out_dir, '10_vpiL_per_bias.png'), config)
        return f'OK  {date}/{band}-band/{wafer}/({row},{col})'
    except Exception as e:
        return f'FAIL {date}/{band}-band/{wafer}/({row},{col}): {e}'
