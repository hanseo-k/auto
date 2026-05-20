# 자동분석폴더 — MZM 4-wafer 분석

HY202103 MZM 측정 XML 을 자동으로 파싱해서 ER, IL, V_π 를 추출하고,
웨이퍼맵 / 1D 분포 / Robust Z 신뢰도 맵 / 날짜별 분석을 생성하는 파이프라인.

## 폴더 구조

```
.
├── run.py                # 메인 실행 파일
├── requirements.txt      # 패키지 의존성
├── src/                  # 분석 모듈
│   ├── xml_loader.py             # XML 파싱
│   ├── extract_er.py             # ER 추출 (sensitivity test 검증된 16 nm 윈도우)
│   ├── extract_il.py             # IL 추출
│   ├── extract_vpi.py            # V_π 추출 (slope filter + vpi_status 꼬리표)
│   ├── extract_passive_params.py # coupler split ratio, MZM section loss
│   ├── outlier_detect.py         # 물리바운드 + Robust Z outlier
│   ├── csv_export.py             # 결과 폴더 / CSV 저장
│   ├── plot_common.py            # 플롯 공통 헬퍼
│   ├── wafer_map.py              # 웨이퍼맵 (Delaunay surface)
│   ├── plot_1d.py                # 1D 분포 (IQR 박스)
│   ├── plot_1d_mad.py            # 1D 분포 (MAD 박스)
│   ├── zscore_map.py             # Robust Z-score 격자맵
│   ├── decompose_variation.py    # systematic / random 분해 (Xing 2023)
│   ├── analyze_by_date.py        # 날짜별 분석 + 물리바운드 위반 강조
│   ├── investigate.py            # 9개 진단 한번에 실행 (개발/검증용)
│   └── sensitivity_test.py       # ER 윈도우 / V_π slope filter sensitivity
├── data/                 # 입력 데이터 (HY202103 은 외부 경로 참조)
├── doc/
│   ├── fig_*.png              # 방법론 그림
│   └── investigation/         # 진단 결과 그림
└── res/
    ├── csv/                   # 최신 CSV/XLSX (git 추적)
    ├── figures/               # 최신 PNG (git 추적)
    └── <YYYY-MM-DD_HH-MM-SS>/ # 매 실행 결과 (gitignored)
```

## 실행

```bash
python3 run.py
```

실행 시:
1. `DATA_ROOT` 경로의 모든 XML 자동 탐색 (dedup: 같은 다이는 최신만)
2. 멀티코어 병렬로 다이 추출 + 플롯 12개 동시 생성
3. 결과를 `res/<timestamp>/` 에 저장
4. 최신 결과를 `res/csv/`, `res/figures/` 로 복사
5. 날짜별 분석 (`data_by_date.csv` + `.xlsx` + `by_date_summary.png`)
6. GitHub `main` 브랜치로 자동 push

## 의존성

`requirements.txt`. 주요: numpy, pandas, scipy, matplotlib, openpyxl

---

# Physical Bounds — 기준과 근거

`outlier_detect.PHYSICAL_BOUNDS` 는 "물리적으로 가능한 값" 의 hard limit.
이 범위 밖이면 **장비/추출 알고리즘의 오류**로 간주하고 `is_problematic`
플래그를 띄움. z-score (상대비교) 와 별개로 동작.

## 디바이스 가정

HY202103 는 표준 **TE-mode Si depletion MZM** 로 가정:
- single-arm drive (push-pull 아님)
- MMI splitter/combiner, single-stage
- polarization filter, cascaded MZI, 추가 DC bias section 없음

## 1. ER (Extinction Ratio) — **10 ~ 45 dB**

| | 값 | 근거 |
|--|---|----|
| 하한 | 10 dB | working device 의 최소선. 그 이하는 MMI 균형 깨짐/도파로 손상 |
| 상한 | 45 dB | 우리 ER 정의(peak−null across biases) + HY202103 실측 분포 + 측정 아티팩트 margin |

문헌상 standard fixed-bias ER 한계는 ~30 dB (Witzens 2018) 지만, 우리 ER 정의는
"모든 바이어스 중 best peak − 모든 바이어스 중 best null" 이라 더 크게 나옴.
HY202103 C-band 실측 32~40 dB 를 정상으로 포함하기 위해 상한 45 dB 채택.

## 2. IL (Insertion Loss, ON-state) — **−15 ~ −1 dB**

| | 값 | 근거 |
|--|---|----|
| 상한 | −1 dB | 비물리적. 커플링/MMI/도파로 손실 합 ≥ 1 dB |
| 하한 | −15 dB | 표준 working device 의 하한 |

## 3. V_π (Half-wave Voltage) — **2 ~ 60 V**

| | 값 | 근거 |
|--|---|----|
| 하한 | 2 V | V_π·L = 1 V·cm × L_max = 5 mm |
| 상한 | 60 V | V_π·L = 3 V·cm × L_min = 0.5 mm |

60 V 초과는 산술적 폭주 (`dλ/dV → 0`) → **slope filter 가 먼저 NaN 처리**.

## Physical Bounds — Empirical Validation

처음에 ER 상한을 35 dB 로 잡았을 때 HY202103 C-band 다이 28 개가 모두
outlier 처리됨 (실측 median 37 dB). 이는 HY202103 의 실제 성능이 문헌 평균보다
높음 + 우리 ER 정의가 fixed-bias ER 보다 큰 값을 줌을 의미.  → 45 dB 로 완화.

**일반 권장:** 새 디바이스 데이터에 적용할 때는 문헌 bound 를 1차 추정으로
쓰되, 실측 분포를 보고 검증/조정 필수.

## 인용

**디바이스 / 물리 한계:**
1. **Soref & Bennett**, "Electrooptical effects in silicon", *IEEE J. Quantum Electron.* **23**, 123–129 (1987).
2. **Reed et al.**, "Silicon optical modulators", *Nature Photonics* **4**, 518–526 (2010).
3. **Patel et al.**, "Design, analysis, and transmission system performance of a 41 GHz silicon photonic modulator", *Opt. Express* **23**, 14263 (2015).
4. **Witzens**, "High-Speed Silicon Photonics Modulators", *Proc. IEEE* **106**, 2158–2182 (2018).

**Wafer-scale 분석 / 공정 변동:**
5. **Selvaraja et al.**, "Process variation in silicon photonic devices", *Appl. Opt.* **52**, 7638 (2013) — classic reference for Si photonic process variation magnitudes.
6. **Xing, Dong, Khan, Bogaerts**, "Capturing the Effects of Spatial Process Variations in Silicon Photonic Circuits", *ACS Photonics* **10**, 928 (2023) — 본 프로젝트의 variation decomposition 방법론의 출처.

**Multi-parameter 추출:**
7. **Xu et al.**, "Optical and geometric parameter extraction across 300-mm photonic integrated circuit wafers", *APL Photonics* **9**, 016104 (2024) — 본 프로젝트의 passive 파라미터 (split ratio, loss) 추출 영감.

**Outlier 검출 통계:**
8. **Iglewicz & Hoaglin**, "How to Detect and Handle Outliers", *ASQC Quality Press* (1993) — Modified Z-score / MAD 의 원전.

---

# Why Reverse Bias? — MZM 작동 원리의 핵심

광통신용 Si MZM 은 **거의 모두 reverse bias (depletion mode) 만 사용**.
우리 측정 데이터도 V = −2, −1.5, −1, −0.5, 0, +0.5 V (5 reverse + 1 forward).
Forward 가 1 개만 있는 건 **검증용**, 운용은 reverse 영역.

## Depletion (reverse) vs Injection (forward) 비교

| | **Reverse (depletion)** | **Forward (injection)** |
|---|----|----|
| 메커니즘 | depletion region 너비 조절 | carrier 직접 주입 |
| 응답 속도 | RC 회로 시정수 (~ps) → **≥ 50 GHz** | carrier lifetime (~ns) → ≤ 1 GHz |
| 전력 소비 | ~ pJ/bit (capacitor 충방전만) | mW 단위 (정상상태 전류) |
| 발열 | 거의 없음 | 큼 |
| 선형성 | √V 비례, 거의 선형 | exp(V/V_T) 강한 비선형 |
| 통신용 사용 | ✅ **표준** | ❌ 너무 느림 + 발열 |

## 왜 reverse 가 빠른가
**Reverse bias = 전기장 인가만으로 carrier 가 즉시 재분포** (charge redistribution, no actual current).
시정수가 RC 회로 의 그것: τ = RC ≈ 5 ps → bandwidth ~ 30 GHz.

반면 **forward = carrier 가 실제로 도파로로 이동해야 함** (diffusion).
Minority carrier lifetime ~ ns → bandwidth 수십 MHz.  100G 광통신엔 못 씀.

## 왜 reverse 가 선형인가
Depletion width 가 `W ∝ √(V_built-in − V)` 의 부드러운 함수.
도파로 안 평균 carrier density 변화도 부드러움 → **dλ/dV ≈ const** 의 선형 거동.

Forward 의 carrier injection 은 Shockley 식대로 **지수 함수** → 비선형 distortion.

→ 우리 `extract_vpi.py` 가 reverse 만 사용 (`rev_biases = [v for v in biases if v <= 0]`).

## Linearity R² — 선형성 정량화

V_π 추출의 신뢰성은 **dλ/dV 곡선이 얼마나 직선인가** 에 달림.
완벽한 depletion 동작이면 R² ≈ 1, 비선형/노이즈가 끼면 R² 가 떨어짐.

`extract_vpi` 가 각 다이마다 `linearity_R2` 컬럼을 반환:
- **R² > 0.99**: 깨끗한 선형 동작 (이상적)
- **R² 0.95 ~ 0.99**: 정상 (약간의 측정 노이즈)
- **R² < 0.9**: 비선형 또는 측정 품질 저하 → V_π 신뢰도 낮음
- **R² < 0.5**: 측정 거의 망가짐 (slope_filter 가 잡지 못한 케이스라도 의심)

이 컬럼으로 **단순 "측정됨/망가짐" 이 아니라 신뢰도 등급** 매김 가능.

---

# Slope Filter — V_π 폭주 방지

`extract_vpi.MIN_SLOPE_PM_PER_V = 10 pm/V`

## 왜 필요한가

V_π = FSR / (2·|dλ/dV|).  바이어스 인가했는데 null 파장이 거의 안 움직이면
(dλ/dV ≈ 0) V_π 가 산술적으로 폭주.  이건 **디바이스 자체의 V_π 가 무한대**가
아니라 **측정이 망가졌다는 신호** (probe contact 불량, 케이블 단선, 측정
소프트웨어 버그 등).

## 임계값 10 pm/V 의 근거

HY202103 의 정상/망가짐 데이터에서 검증 (`src/sensitivity_test.py`):

| 측정 상태 | dλ/dV (pm/V) median | |min| |
|----------|--------------------|------|
| 정상 (06-03 D23/D24) | -119  | 107 |
| 망가짐 (05-31 D23/D24) | +0.3 | 0.06 |

정상의 최저 107 과 망가짐의 최고 4.7 사이 안전한 분기점 → **10 pm/V** 채택.

## 적용 효과

```
변경 전: 망가진 28개 중 28개가 V_π = 1062 ~ 78633 V 로 폭주
변경 후: 망가진 28개 모두 V_π = NaN (정상 70개는 100% 영향 없음)
```

dλ/dV 값은 그대로 보고함 (CSV 의 `dlam_dV_pm_per_V` 컬럼) — 진단 가능.

## `vpi_status` — 명시적 상태 꼬리표

`extract_vpi` 는 V_π 추출 결과와 함께 **`vpi_status`** 라는 명시적 꼬리표를 같이
반환함. 추출이 왜 실패/성공했는지를 추론하지 않고도 알 수 있게 하는 장치.

| `vpi_status` | 의미 | `vpi_V` |
|--------------|------|---------|
| `ok` | 정상 추출 | 실제 값 |
| `slope_filter` | \|dλ/dV\| < 10 pm/V (측정 망가짐) | NaN |
| `no_nulls` | deep null < 2 개 — FSR 추출 실패 | NaN |
| `no_slopes` | 모든 null tracking 이 점프로 reject | NaN |
| `few_biases` | reverse-bias 데이터 < 3 개 | NaN |
| `no_sweeps` | sweep 데이터 자체 없음 | NaN |

CSV/XLSX 모듈은 이 꼬리표를 받아서 `extract_vpi.status_to_reason()` 을 호출,
사람 친화적 reason 문자열을 만들어 `reason_Vpi_V` 컬럼에 기록:

```
broken: |dλ/dV| < 10 pm/V (slope filter, |dλ/dV|=0.06)
```

**왜 꼬리표가 필요한가:**
이전엔 `Vpi_V` 가 NaN 이면 reason 컬럼이 그냥 `missing` 으로 표시됐는데,
그게 (a) slope filter 발동 (b) FSR 추출 실패 (c) 데이터 부족 중 어느 건지
모호했음. 명시적 꼬리표로 한 번에 분리 가능.

코드 흐름:
```
extract_vpi.extract_vpi(die)
    → { fsr_nm, dlam_dV_pm_per_V, vpi_V, vpi_status }
                                          ↓
run.py / analyze_by_date.py
    → row dict 에 그대로 보존
                                          ↓
analyze_by_date._vpi_reason
    → status_to_reason(vpi_status, dlam_dV_pm_per_V)
    → 'broken: |dλ/dV| < 10 pm/V (slope filter, |dλ/dV|=0.06)'
                                          ↓
xlsx 의 reason_Vpi_V 셀에 표시 (빨간 배경)
```

---

# ER 윈도우 — Sensitivity 검증

## 왜 16 nm 인가 (C-band)

C-band FSR ≈ 14.3 nm.  과거 14 nm 윈도우는 FSR 보다 좁아 다이별로 null 이
0.97~0.99 개만 포함 → 가끔 0 개 포함하는 다이 발생 → ER 불안정.

`src/sensitivity_test.py` 결과:

| 윈도우 폭 | median ER | std | n>45 |
|----------|-----------|-----|------|
| 14 nm (이전) | 36.87 | **1.47** | 0 |
| **16 nm (현재)** | 37.11 | **1.24** ✅ | 0 |
| 18 nm | 37.31 | 1.40 | 0 |
| 22 nm | 38.03 | 1.36 | 0 |
| 36 nm (극단) | 39.53 | **2.48** ❌ | 2 |

→ 16 nm 가 sweet spot.  너무 넓히면 ALIGN reference 가 신뢰성 떨어지는 밴드
   엣지를 포함해서 노이즈 + outlier 발생.

O-band 는 FSR 9.8 nm 이라 14 nm 윈도우에 이미 1.4 개 null 안정 포함 → 변경 없음.

---

# Passive 파라미터 — Coupler Split Ratio & MZM Section Loss

APL Photonics 2024 ("Optical and geometric parameter extraction across 300mm
PIC wafers") 의 방법론을 우리 단순화 버전으로 적용. 추가 측정 없이 기존 spectrum
에서 추출.

## 물리 모델

비균형 MZI 의 transfer function:
```
T(λ) = a² + b² + 2ab·cos(Δφ)
```
여기서 `a`, `b` 는 두 팔의 amplitude.  ER 의 linear 비로부터 amplitude 비
`k = b/a` 를 닫힌형식으로 유도 가능:
```
√ER_linear = (1+k)/(1-k)
k          = (√ER_linear − 1) / (√ER_linear + 1)
```

## 출력 컬럼

| 컬럼 | 의미 | 이상적 |
|------|------|--------|
| `amplitude_ratio_k` | b/a (두 팔 amplitude 비) | 1.0 |
| `power_split_ratio` | k² (power) | 1.0 |
| `imbalance_dB` | −20·log₁₀(k) | 0 dB |
| `mzm_loss_dB` | −T_dev_peak (양수) | 0 dB |

## HY202103 결과

| Wafer-Band | Imbalance (dB) | MZM loss (dB) | 해석 |
|------------|----------------|---------------|------|
| D07-C, D08-C | ~0.24 | 0.2~0.3 | 매우 균형 잡힌 splitter, 저손실 |
| D08-O, D23-O, D24-O | ~0.44 | ~1.0 | 약간 불균형, 손실 약간 큼 |

→ C-band 가 ER 이 더 높은 이유 부분 설명됨 (splitter 균형이 더 좋음).

---

# Variation Decomposition — Spatial Trend vs Random Noise

Xing et al. ACS Photonics 2023 의 방법론을 우리 데이터에 적용.
각 (Wafer × Band) 그룹의 다이 metric 을 두 컴포넌트로 분해:
```
observed(R, C)  =  systematic(R, C)  +  random(R, C)
                   (2D 다항식 fit)      (잔차)
```

## 왜 분해하는가

기존 `zscore_map.py` 는 두 컴포넌트를 통째로 보여줘서 "공정 gradient" 와
"random 변동" 을 구분 못 함.  분해하면:

- **systematic 패턴이 크면** (R² 높음) → **공정 분포 문제** (lithography focus,
   etch uniformity, deposition gradient).  공정 엔지니어가 손볼 수 있음.
- **random 이 크면** (R² 낮음) → **다이별 random 변동** 또는 측정 노이즈.
   공정 control 보다는 redundancy / statistical yield 관리 영역.

## HY202103 결과 (R² = systematic 이 분산을 얼마나 설명하는가)

| Wafer-Band | ER R² | IL R² | Vpi R² | 해석 |
|------------|-------|-------|--------|------|
| D08-O | 0.48 | **0.90** | 0.72 | IL/Vpi 의 90% 가 공정 gradient |
| D23-O | 0.64 | 0.83 | **0.85** | Vpi 강한 systematic 패턴 |
| D07-C | 0.38 | 0.55 | **0.10** | Vpi 거의 random (gradient 없음) |
| D24-O | 0.31 | 0.81 | 0.40 | IL 만 systematic |
| D08-C | 0.51 | 0.67 | 0.29 | 모두 적당히 systematic |

**핵심 시사점:**
- O-band wafer 들이 IL, Vpi 에서 **강한 공간 gradient** 보임 → 공정 uniformity 문제
- C-band 의 Vpi 는 **거의 random** → 다이별 random 변동이 dominant
- IL 이 가장 systematic 한 패턴 → 손실은 위치에 의존적인 공정 변동에 가장 민감

## 출력

`res/figures/decompose_<ER|IL|Vpi>.png`: 각 (Wafer × Band) 행마다 3 패널
[observed | systematic | random], R² 동시 표기.

## Reference

Y. Xing, J. Dong, U. Khan, W. Bogaerts, "Capturing the Effects of Spatial
Process Variations in Silicon Photonic Circuits", *ACS Photonics* **10**, 928 (2023).

---

# Outlier Detection 방법론

## Robust Z-score (per Wafer × Band)

```
σ_robust = 1.4826 × MAD
z'       = (x − median) / σ_robust
outlier  if |z'| > 3
```

## 왜 3-sigma 가 아닌가
- 3-sigma 는 정규성 가정 + 충분한 표본 요구. 웨이퍼당 14 다이로는 부적합
- median/MAD 는 outlier 에 의해 자체가 오염되지 않음 (robust)
- σ_robust = 1.4826 × MAD 는 정규분포에서 σ 와 같으므로 3-sigma 와 직접 비교 가능

참고: Iglewicz & Hoaglin, "How to Detect and Handle Outliers", ASQC, 1993.

## 왜 그룹 단위 (per Wafer-Band) 인가 — Hampel(이웃 비교) 가 아닌 이유
- 웨이퍼당 14 다이로는 8-이웃 Hampel 의 표본이 너무 작음 (3~8 개, 가장자리 더 작음)
- 그룹 전체 (~14 개) 의 median/MAD 가 훨씬 안정

---

# 데이터 품질 진단 결과 (1~9 항목)

전체 9개 항목 진단 완료 (`src/investigate.py` + `src/sensitivity_test.py`).
주요 발견 요약:

## 🚨 즉시 조치

| # | 항목 | 발견 | 조치 |
|---|------|------|------|
| 1 | 2019-05-31 V_π 폭주 28개 | **dλ/dV ≈ 0** (점프 아님, 측정 자체 망가짐) | slope filter 10 pm/V 추가 |
| 6 | ER 윈도우 좁음 | C-band 14nm 가 FSR=14.3nm 보다 좁아 null 0.98개만 포함 | 16 nm 로 확장 |

## ⚠️ 한계점 (코드 변경 없음, 미래 작업)

| # | 항목 | 발견 |
|---|------|------|
| 2 | 재현성 데이터 없음 | 같은 다이 정상 재측정 0 개 → 측정 신뢰도 정량화 불가. **추후 측정 캠페인에 포함 필요** |
| 4 | dedup 가정 안전 | 운 좋게 망가진 측정은 항상 옛날 → 최신=최선 가정 통과. 하지만 항상 그렇진 않음 |

## ✅ 확인된 사실

| # | 항목 | 결과 |
|---|------|------|
| 3 | Width 분포 | 380nm (O-band) / 450nm (C-band) — Band 와 1:1 매핑 → 현재 그룹화 (Wafer × Band) 가 자동으로 width 분리 |
| 5 | D08 C/O 양 밴드 | 같은 (Row,Col) 14 개에 다른 디바이스 (380/450nm) 가 함께 있음 → 별개 행 취급 정당 |
| 8 | FSR 일관성 | C ≈ 14.3±0.1, O ≈ 9.8±0.07 — 매우 일관적 |
| 9 | dλ/dV 정상값 | C: ≈ -210, O(D08): ≈ -174, O(D23/24): ≈ -120 pm/V — 모두 0 에서 멀리 안전 |

진단 그림: `doc/investigation/01_2019-05-31_diagnosis.png` ~ `09_dlam_dv_map.png`,
`sens_A_vpi.png`, `sens_B_er_window.png`

---

# 출력 컬럼

## `res/csv/data.csv` (dedup 된 main 결과)

| 컬럼 | 설명 |
|------|------|
| Wafer, Band, Row, Col, Width_nm | 다이 식별자 |
| ER_dB | Extinction Ratio (peak − null, 고정 윈도우 16/14 nm) |
| IL_dB | Insertion Loss (V=−1V, ±5nm peak transmission) |
| Vpi_V | V_π (FSR / 2·\|dλ/dV\|), slope filter 적용됨 |
| FSR_nm, dlam_dV_pm_per_V | V_π 계산 중간값 (진단용) |
| is_outlier_* | 항목별 outlier 플래그 (physical bound OR Robust Z) |
| robust_z_* | Robust Z-score |
| is_trusted | 세 항목 모두 trusted 인가 |

## `res/csv/data_by_date.csv` / `.xlsx` (날짜별, dedup 없음)

| 컬럼 | 설명 |
|------|------|
| Date | 측정 날짜 (YYYY-MM-DD) |
| Wafer, Band, Row, Col, Width_nm | 다이 식별자 |
| ER_dB, IL_dB, Vpi_V | 측정값 |
| reason_X | 물리바운드 위반 사유 (예: `over: 5478.34 > 80.0`) |
| out_of_bound_X | bool, reason 비어있지 않으면 True |
| is_problematic | 세 항목 중 하나라도 위반하면 True |

xlsx 버전은 문제 셀에 **빨간 배경** 적용.

---

# 디자인 결정 노트

## 왜 dedup 하는가?
같은 `(Wafer, Row, Col, Band)` 에 대해 가장 최신 측정만 남김.  근거: 재측정은
보통 이전 측정에 문제가 있어서 수행됨 → 최신이 가장 신뢰할 만함.
**한계:** 이 가정이 깨지면 main 분석이 망가짐.  현재 데이터에선 안전한 것 확인됨 (#4).

## 왜 dλ/dV 값을 NaN 으로 안 바꾸고 그대로 두는가?
Slope filter 발동했을 때 V_π 만 NaN, `dlam_dV_pm_per_V` 는 측정값 그대로 보고.
이유: 진단을 위해서.  dλ/dV 값을 보면 *왜* V_π 가 NaN 인지 즉시 알 수 있음
(0 근처면 측정 망가짐, 정상값이면 다른 이유).
