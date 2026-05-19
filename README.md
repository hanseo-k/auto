# 자동분석폴더 — MZM 4-wafer 분석

HY202103 MZM 측정 XML을 자동으로 파싱해서 ER, IL, V_π 를 추출하고,
웨이퍼맵 / 1D 분포 / Robust Z 신뢰도 맵을 생성하는 파이프라인.

## 폴더 구조

```
.
├── run.py              # 메인 실행 파일
├── requirements.txt    # 패키지 의존성
├── src/                # 분석 모듈
│   ├── xml_loader.py        # XML 파싱
│   ├── extract_er.py        # ER 추출
│   ├── extract_il.py        # IL 추출
│   ├── extract_vpi.py       # V_π 추출
│   ├── outlier_detect.py    # 물리바운드 + Robust Z outlier
│   ├── csv_export.py        # 결과 폴더 / CSV 저장
│   ├── wafer_map.py         # 웨이퍼맵 (Delaunay surface)
│   ├── plot_1d.py           # 1D 분포 (IQR 박스)
│   ├── plot_1d_mad.py       # 1D 분포 (MAD 박스)
│   ├── zscore_map.py        # Robust Z-score 격자맵
│   └── analyze_by_date.py   # 날짜별 분석 + 물리바운드 위반 강조
├── data/               # 입력 데이터 (HY202103은 외부 경로 참조)
├── doc/                # 방법론 그림
└── res/
    ├── csv/                 # 최신 CSV/XLSX (git 추적)
    ├── figures/             # 최신 PNG (git 추적)
    └── <YYYY-MM-DD_HH-MM-SS>/  # 매 실행 결과 (gitignored)
```

## 실행

```bash
python3 run.py
```

실행 시:
1. `DATA_ROOT` 경로의 모든 XML 자동 탐색 (dedup: 같은 다이는 최신만)
2. 멀티코어 병렬로 다이 추출 + 플롯 생성
3. 결과를 `res/<timestamp>/` 에 저장
4. 최신 결과를 `res/csv/`, `res/figures/` 로 복사
5. **날짜별 분석** (`data_by_date.csv` + `.xlsx` + `by_date_summary.png`)
6. GitHub `main` 브랜치로 자동 push

## 의존성

`requirements.txt` 참고. 주요 패키지:
- numpy, pandas, scipy, matplotlib
- openpyxl (xlsx 빨간 셀 강조용)

---

# Physical Bounds — 기준과 근거

`outlier_detect.PHYSICAL_BOUNDS` 는 "물리적으로 가능한 값" 의 hard limit 입니다.
이 범위를 벗어난 값은 **장비/추출 알고리즘의 오류**로 간주하고 `is_problematic`
플래그를 띄웁니다. z-score (상대비교) 와 별개로 동작합니다.

## 디바이스 가정

HY202103 는 표준 **TE-mode Si depletion MZM** 로 가정:
- single-arm drive (push-pull 아님)
- MMI splitter/combiner, single-stage
- polarization filter, cascaded MZI, 추가 DC bias section 없음

이 가정이 깨지면 (예: push-pull 이거나 cascaded) 일부 상한이 더 커질 수 있음.

## 1. ER (Extinction Ratio) — **10 ~ 45 dB**

| | 값 | 근거 |
|--|---|----|
| 하한 | 10 dB | working device 의 최소선. 그 이하는 MMI 균형이 깨졌거나 도파로 손상으로 판단 |
| 상한 | 45 dB | 우리 ER 정의(peak−null across biases) + HY202103 실측 분포 + 측정 아티팩트 margin |

**왜 이 값인가:**
- 표준 fixed-bias ER 의 산업적 상한 ≈ 30 dB (Witzens 2018, Table II) — single-arm MMI 의 splitter imbalance 한계
- Patel 2015 의 standard MMI Si MZM: ~25 dB
- **우리 ER 정의가 다름:** `extract_er.py` 는 *모든 바이어스의 peak − 모든 바이어스의 null* 를 fixed window 안에서 계산. 즉 "ever-achieved peak" 와 "ever-achieved null" 의 차이이므로 fixed-bias ER 보다 자연스럽게 크게 측정됨
- **HY202103 실측 분포 (`res/csv/data.csv` 기준):**
  - C-band: 32 ~ 40 dB (median 37)
  - O-band: 29 ~ 35 dB (median 32)
- 상한을 35 dB로 잡으면 C-band 정상 측정 전체가 outlier 처리됨 → 비합리
- 상한 45 dB 는 실측 max(40.34) 위 ~5 dB margin. 그 이상은 노이즈 floor 인공값으로 판단
- 하한 10 dB 는 working device 의 의미있는 minimum (10 dB 미만은 modulator 로 사용 불가)

**검증 절차:** ER bound 는 문헌만으로는 부족하다는 점에 주의. 새로운 디바이스
세대를 분석할 때는 이 코드를 한 번 돌리고 `res/csv/data.csv` 의 ER 분포를 보고
실측 max + 5 dB 정도로 상한을 재조정하는 것이 권장됨.

## 2. IL (Insertion Loss, ON-state) — **−15 ~ −1 dB**

| | 값 | 근거 |
|--|---|----|
| 상한 | −1 dB | 비물리적. 커플링(1~3 dB) + MMI(0.5 dB) + 도파로 ≥ 1 dB. 0 dB 부근은 측정/캘리브레이션 오류 |
| 하한 | −15 dB | 표준 working device 의 하한. 그 이하는 도파로 단절/커플링 실패 |

**왜 이 값인가:**
- 표준 Si MZM ON-state IL ≈ 4 ~ 10 dB (Reed 2010 review, Witzens 2018)
- 커플러 + MMI + phase shifter 합치면 최소 4 dB
- ALIGN 레퍼런스 차감 후 device IL 만 보더라도 ≥ 1 dB
- 0 dB 에 가까우면 → 레퍼런스 정규화 실패 가능성

**측정 정의:** V = −1 V 바이어스에서 ± 5 nm 윈도우 내 peak transmission.
ALIGN 차감 안 함 (XML 의 IL 값이 이미 device IL 로 기록된 것으로 가정).

## 3. V_π (Half-wave Voltage) — **2 ~ 60 V**

| | 값 | 근거 |
|--|---|----|
| 하한 | 2 V | V_π·L = 1 V·cm × L_max = 5 mm 의 이론적 하한 |
| 상한 | 60 V | V_π·L = 3 V·cm × L_min = 0.5 mm 의 이론적 상한 |

**왜 이 값인가:**
- Si depletion MZM 의 V_π·L 은 **plasma dispersion 효과**로 결정되며, Soref & Bennett (1987) 의 이론과 후속 측정 결과 **1 ~ 3 V·cm** 가 산업적 표준 (Witzens 2018)
- 우리 디바이스의 phase shifter 길이는 0.5 ~ 5 mm 범위로 가정 (보통 1~3 mm)
- 60 V 초과는 **산술적 폭주**: V_π = FSR / (2 · |dλ/dV|) 에서 dλ/dV → 0 이면 V_π → ∞. 실제 V_π 가 아니라 **null 트래킹/fit 실패**의 신호
- 우리 `data_by_date.csv` 의 2019-05-31 측정 28 개에서 Vpi 1062 ~ 78633 V 가 발견됐는데, 이게 바로 이 케이스

**측정 정의:** V_π = FSR / (2 · |dλ_null / dV|).
V ≤ 0 구간의 null 위치를 parabolic fit 으로 추적 → 선형 slope = dλ/dV.

## Physical Bounds — Empirical Validation

처음에 문헌만 보고 ER 상한을 35 dB 로 잡았을 때, HY202103 의 C-band 다이 28 개가
**모두 outlier 처리**되는 문제가 발생함 (실측 ER median = 37 dB).

이는 두 가지를 시사함:
1. **HY202103 의 실제 ER 성능이 문헌 평균보다 높음** — high-uniformity MMI / 긴 phase shifter 등 high-ER 친화적 설계로 보임
2. **우리 ER 정의가 fixed-bias ER 보다 큰 값을 줌** — peak/null 둘 다 모든 바이어스 중 best 를 골라서

→ ER 상한을 45 dB 로 완화. 일반적으로 새 디바이스 데이터에 적용할 때는 문헌
   bound 를 1차 추정으로 쓰되, 실측 분포를 보고 검증 / 조정해야 함.

## 인용

1. **Soref & Bennett**, "Electrooptical effects in silicon", *IEEE J. Quantum Electron.* **23**, 123–129 (1987).
   - Si plasma dispersion 의 원전. V_π·L 한계의 출발점.
2. **Reed et al.**, "Silicon optical modulators", *Nature Photonics* **4**, 518–526 (2010).
   - Si MZM IL / ER 의 표준 review.
3. **Patel et al.**, "Design, analysis, and transmission system performance of a 41 GHz silicon photonic modulator", *Opt. Express* **23**, 14263 (2015).
   - 표준 MMI single-arm Si MZM 의 측정값 reference.
4. **Witzens**, "High-Speed Silicon Photonics Modulators", *Proc. IEEE* **106**, 2158–2182 (2018).
   - 종합 review. Table II 에 ER / IL / V_π 의 산업 평균치 정리.

---

# 출력 컬럼 정의

## `res/csv/data.csv` (dedup된 main 결과)

| 컬럼 | 설명 |
|------|------|
| Wafer, Band, Row, Col, Width_nm | 다이 식별자 |
| ER_dB | Extinction Ratio (peak − null, 고정 윈도우) |
| IL_dB | Insertion Loss (V=−1V, ±5nm peak transmission) |
| Vpi_V | V_π (FSR / 2·\|dλ/dV\|) |
| FSR_nm, dlam_dV_pm_per_V | V_π 계산 중간값 |
| is_outlier_* | 항목별 outlier 플래그 (physical bound OR Robust Z) |
| robust_z_* | Robust Z-score (median/MAD, per Wafer-Band 그룹) |
| is_trusted | 세 항목 모두 trusted 인가 |

## `res/csv/data_by_date.csv` / `.xlsx` (날짜별, dedup 없음)

| 컬럼 | 설명 |
|------|------|
| Date | 측정 날짜 (YYYY-MM-DD) |
| Wafer, Band, Row, Col, Width_nm | 다이 식별자 |
| ER_dB, IL_dB, Vpi_V | 측정값 |
| reason_X | 물리바운드 위반 사유 (예: `over: 5478.34 > 80.0`) |
| out_of_bound_X | bool, reason 이 비어있지 않으면 True |
| is_problematic | 세 항목 중 하나라도 위반하면 True |

xlsx 버전은 `is_problematic=True` 셀과 `reason_X` 비어있지 않은 셀에 **빨간 배경** 적용.

---

# 분석 로드맵 (To Do)

지금까지 한 deep-thinking 결과 우선순위:

## 🚨 즉시 조치 (데이터 품질)

1. **2019-05-31 망가진 측정 원인 진단**
   - V_π 폭주 (1062 ~ 78633 V) 의 근본 원인을 raw 스펙트럼에서 직접 확인
   - 후보: null 트래킹이 옆 null 로 점프 / parabolic fit 실패 / 측정 SNR
   - 그 날 데이터를 명시적으로 제외하거나, 추출 알고리즘을 수정

2. **D08 반복측정 재현성 검증**
   - D08 은 05-26, 07-12 등 여러 날짜에 측정됨 → 동일 다이를 다른 날 측정한 값의 일관성 확인
   - ER / IL / V_π 의 measurement-to-measurement variation 정량화
   - 장비 신뢰도 / 디바이스 안정성 둘 다 확인 가능

3. **Width_nm 분포 점검**
   - `MZMCTE_LULAB_450_500` 의 450 이 width — 모든 다이가 450nm 인지, 아니면 여러 width 가 섞여있는지 확인
   - 다른 width 면 outlier 검출 그룹화 (Wafer × Band) 가 부적절. (Wafer × Band × Width) 로 분할 필요

## 💭 방법론 점검

4. **dedup "최신 = 최선" 가정 재고**
   - 운 나쁘게 망가진 측정이 최신이면 main 분석이 그대로 망가짐
   - 옵션: 다중 측정의 median 사용, 또는 "검증된" 측정만 사용

5. **D08 의 C/O 양 밴드 처리**
   - 같은 물리 다이가 2 개 행으로 분리되어 있음 — 의도된 구조인지 확인

6. **ER 윈도우 검증**
   - C-band 1546-1560 / O-band 1306-1320 의 14 nm 폭이 모든 다이에서 같은 수의 null 을 포함하는지 확인 (FSR 일관성 sanity check)

## 🔧 빠진 시각화 (한 줄 추가)

7. **IV 곡선 시각화**
   - `xml_loader` 가 iv_V/iv_I 를 추출하지만 한 번도 안 씀
   - depletion MZM 의 IV 는 측정 건전성의 1차 단서 (누설전류 큰 다이는 ER 측정 자체가 깨짐)

8. **FSR 일관성 맵**
   - FSR 이 다이/웨이퍼별로 얼마나 일관적인가? 튀는 다이는 구조/측정 문제

9. **`dlam_dV` 웨이퍼맵**
   - V_π 폭주의 직접 원인 지표 (0 에 가까우면 V_π → ∞)
   - 망가진 측정 진단의 핵심 도구

---

## 디자인 결정 노트

### 왜 dedup 하는가?
`xml_loader.find_all_xmls` 는 같은 `(Wafer, Row, Col, Band)` 에 대해 **가장 최신** 측정만 남김.
근거: 동일 다이의 재측정은 보통 이전 측정에 문제가 있어서 수행됨 → 최신이 가장 신뢰할 만함.
**한계:** 이 가정이 깨지는 케이스가 있을 수 있음 (위 로드맵 #4 참고).

### 왜 Robust Z (MAD) 인가? — 3-sigma 가 아닌 이유
- 3-sigma 는 **정규성 가정** 과 충분한 표본을 요구. 우리는 웨이퍼당 ~14 다이로 표본이 작음
- median 과 MAD 는 outlier 에 의해 자체가 오염되지 않음 (robust)
- σ_robust = 1.4826 × MAD 는 정규분포에서 σ 와 같으므로 3-sigma 와 직접 비교 가능
- 참고: Iglewicz & Hoaglin, "How to Detect and Handle Outliers", ASQC Quality Press, 1993

### 왜 그룹 단위 (per Wafer-Band) 인가? — 이웃 비교(Hampel) 가 아닌 이유
- 웨이퍼당 14 다이로는 이웃 8 개로 MAD 계산하면 표본이 너무 작아 노이즈에 휘둘림
- 가장자리 다이는 이웃이 3-4 개밖에 없어 더 심함
- 그룹 전체 (~14 개) 의 median/MAD 가 더 안정적
