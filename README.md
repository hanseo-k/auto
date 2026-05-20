# Wafer-scale MZM 분석 파이프라인

웨이퍼 단위 광·전기 측정 데이터로부터 Mach-Zehnder Modulator (MZM) 의
주요 파라미터를 자동 추출하고 통계 분석을 수행한다. HY202103 디바이스
세트를 기준으로 설계되었으며, 같은 XML 스키마를 가진 임의의 depletion-mode
실리콘 MZM 데이터에 동일하게 적용된다.

---

## 1. 서론

### 1.1 목적

웨이퍼 단위 실리콘 광공정은 수백 개의 동일 설계 디바이스를 동시 생산하지만
공정 편차, 측정 장비 드리프트, 프로브 접촉 품질의 차이로 인해 다이마다
값의 산포가 발생한다. 이 산포 데이터로부터 다음을 수행할 필요가 있다.

1. 운용에 유의미한 파라미터 (소광비, 삽입손실, 반파전압) 를 모든 다이에
   대해 추출한다.
2. 망가지거나 신뢰성이 낮은 측정과 실제 디바이스 편차를 구분한다.
3. 남은 편차를 공간적으로 체계적인 추세 (공정 엔지니어가 대응할 수 있는
   부분) 와 랜덤 산포 (다중 측정·redundancy 로 대응해야 하는 부분) 로
   분해한다.

이 코드베이스는 위 세 작업을 하나의 재현 가능한 파이프라인으로 처리한다.

### 1.2 범위

입력은 HY202103 스키마를 따르는 XML 파일들이며, 각 파일은 한 다이의 광
스펙트럼 (여섯 개의 바이어스에서 측정) 과 전류-전압 특성을 포함한다.

출력은 다이별 결과 표 (CSV, XLSX), metric 별 웨이퍼 맵, 1D 분포 그림,
Robust Z-score 맵, 측정 날짜별 분석 CSV (위반 항목이 빨간색으로 강조된
Excel 동봉), systematic/random 공간 변동 분해 그림이다.

### 1.3 파이프라인 개요

각 다이마다 세 개의 1차 metric (ER, IL, V_pi) 과 여러 파생 파라미터
(FSR, dlambda/dV, 선형성 R^2, 커플러 imbalance, MZM 구간 손실) 가
계산된다. 다이별 행을 모아 집계한 뒤, 물리적 한계와 (Wafer, Band)
그룹 단위 Robust Z-score 로 outlier 를 플래깅한다. 결과는 디스크에
기록되고 GitHub 원격 저장소에 자동으로 푸시된다.

```
XML 파일
  -> 다이별 추출 (ER, IL, V_pi, passive 파라미터)
  -> outlier 플래깅 (물리바운드 + Robust Z)
  -> 집계
  -> 웨이퍼맵 / 1D 분포 / Z-score 맵 / decomposition / by-date
  -> CSV + XLSX + PNG 를 res/ 에 기록
  -> GitHub 로 자동 commit/push
```

---

## 2. 이론

### 2.1 MZM 동작 원리

MZM 은 입력광을 두 갈래로 분기시키고, 한쪽 팔의 위상을 전압으로 조절한
뒤 두 광을 재결합시킨다. 출력 투과율은

```
T(lambda, V) = |a * exp(j*phi1) + b * exp(j*phi2(V))|^2
             = a^2 + b^2 + 2*a*b*cos(Delta_phi(V))
```

이때 a, b 는 splitter 통과 후 각 팔의 amplitude 이고, 두 팔의 위상차
`Delta_phi` 는 바이어스 V 에 의존한다. 두 팔이 동위상이면 최대 투과,
180 도 위상차이면 최소 투과 (이상적으로 0) 가 된다.

파장에 따라 위상차가 변하므로 투과율은 파장 축에서 주기적인 null 패턴을
보이며, 이 주기 (FSR, Free Spectral Range) 와 바이어스 인가 시 각
null 이 이동하는 속도 (dlambda/dV) 가 반파전압을 결정한다.

```
V_pi = FSR / (2 * |dlambda/dV|)
```

### 2.2 Plasma Dispersion

실리콘 depletion-mode MZM 의 위상 변조는 plasma dispersion 효과
(Soref & Bennett, 1987) 에 의해 발생한다.

```
Delta_n = -8.8e-22 * Delta_N_e - 8.5e-18 * (Delta_N_h)^0.8
```

여기서 `Delta_N_e`, `Delta_N_h` 는 자유 전자/정공 농도 변화량이다. PN
접합에 역방향 바이어스를 인가하면 공핍 영역이 넓어지면서 광 모드 영역의
자유 carrier 가 *감소* 하고 유효 굴절률이 *증가* 한다. 순방향 바이어스
에서는 carrier 가 주입되어 굴절률이 *감소* 한다.

### 2.3 역방향 바이어스가 운용 영역인 이유

상용 고속 실리콘 MZM 은 거의 모두 역방향 바이어스 (depletion mode) 로만
동작한다. 두 영역의 비교를 Table 1 에 정리한다.

**Table 1.** 역방향 (depletion) 과 순방향 (injection) 의 비교.

| 특성 | 역방향 (depletion) | 순방향 (injection) |
|---|---|---|
| 메커니즘 | 공핍 영역 폭 변조 | carrier 직접 주입 |
| 응답 시간 | RC 시정수 (~5 ps) | minority carrier lifetime (~1-100 ns) |
| 대역폭 | 30 GHz 이상 | 10-100 MHz |
| 정상상태 전력 | sub-pJ/bit (capacitive 만) | mW (정상 전류) |
| 선형성 | `W ~ sqrt(V_bi - V)`, 거의 선형 | Shockley 방정식, 지수 비선형 |
| 통신 적합성 | 표준 | ~100 MHz 이상 부적합 |

역방향 바이어스는 정상상태 전류 없이 carrier 의 *공간 분포만* 조절한다.
이 재분포는 capacitive 시정수로 동작하므로 현대 광통신이 요구하는
대역폭을 만족한다. 순방향 바이어스는 실제 carrier 수송을 요구하므로
느리고, 또한 정상 전류가 흘러 큰 전력을 소비한다.

따라서 추출 코드는 `dlambda/dV` 피팅에서 역방향 (V <= 0) 만 사용한다.
순방향 데이터는 raw XML 에 보존되지만 slope 추정에는 포함되지 않는다.

### 2.4 품질 지표로서의 선형성

이상적 depletion-mode MZM 은 `Delta_lambda` 대 `V_bias` 가 엄밀히
선형이다. 선형성 이탈은 다음 중 하나를 의미한다.

- 기대보다 낮은 전압에서 순방향 도통이 시작됨 (직렬 저항 큼 또는
  의도되지 않은 순방향 동작).
- 적당한 역방향 바이어스에서 이미 carrier injection 이 일어남
  (공정 결함).
- 측정 노이즈 또는 null 추적 실패.

따라서 linear fit 의 결정계수 R^2 는 다이별 신호 품질 지표 역할을 한다.
파이프라인은 각 다이에서 성공적으로 추적된 모든 null 의 R^2 중앙값을
`R2_dlam_vs_V` 컬럼에 보고한다.

### 2.5 ER 로부터 커플러 분배비 도출

두 팔의 amplitude 비를 `k = b/a (k <= 1)` 로 두면 선형 ER 은

```
ER_linear = ((1 + k) / (1 - k))^2
```

이를 풀면

```
k = (sqrt(ER_linear) - 1) / (sqrt(ER_linear) + 1)
```

즉 직접 측정되는 ER 한 값에서 splitter 의 amplitude 불균형이 닫힌
형식으로 유도된다. 파이프라인은 각 다이마다 `amplitude_ratio_k`,
`power_split_ratio` (k^2), `imbalance_dB` (-20 log10 k) 를 보고하며,
각 개별 바이어스에서의 imbalance 도 따로 기록한다.

### 2.6 Outlier 검출 전략

서로 독립적인 두 층을 결합한다.

1. **물리바운드.** 각 metric 을 디바이스 물리와 관측 분포로부터 유도된
   hard range 와 비교한다 (5.1 절). 범위 밖 값은 디바이스 편차가 아닌
   *추출 실패* 로 간주하여 by-date CSV 에서 `is_problematic = True`
   로 플래깅한다.

2. **Robust Z-score.** 각 (Wafer, Band) 그룹 안에서 modified Z-score
   `z' = (x - median) / (1.4826 * MAD)` 를 계산한다. `|z'| > 3` 인
   다이를 플래깅한다. Robust Z 는 물리바운드와 독립적이며 둘 중 하나만
   걸려도 outlier 로 표시된다.

Z-score 의 비교 단위를 (Wafer, Band) 그룹으로 잡은 이유는 웨이퍼당
14 개 라는 표본 크기로는 공간 이웃 (Hampel) 기반 비교가 불안정하기
때문이다.

### 2.7 변동 분해

Xing et al. (ACS Photonics 2023) 의 방법에 따라 다이별 metric 을
부드러운 공간 성분과 잔차로 분해한다.

```
observed(R, C) = systematic(R, C) + random(R, C)
```

`systematic` 은 (Row, Col) 의 2차 다항식이며 least squares 로 적합한다.

```
R^2 = Var(systematic) / Var(observed)
```

는 다이 간 변동 중 위치에 의존적인 부분의 비율을 나타낸다. R^2 가
크면 공정 균일성 개선이 효과적이고, 작으면 남은 변동은 본질적으로
랜덤이므로 통계적 대응 (다이 수 증가, redundancy) 이 필요하다.

---

## 3. 방법론

### 3.1 입력 스키마

XML 파일 하나는 (Wafer, DieRow, DieColumn, Band) 한 항목을 담는다.
사용되는 필드는 다음과 같다.

- `TestSiteInfo` 속성 Wafer, DieRow, DieColumn.
- `Modulator` 의 이름이 `ALIGN...` 인 항목: passive 레퍼런스 스펙트럼
  (파장 L, 투과 전력 IL).
- `Modulator` 의 이름이 `MZM...` 인 항목: 측정 대상. 바이어스 별로
  최대 6 개의 `WavelengthSweep` 와 1 개의 `IVMeasurement` 를 제공한다.

파일명 접미사가 `_DCM_LMZC` 이면 C-band (lambda_c = 1550 nm),
`_DCM_LMZO` 이면 O-band (lambda_c = 1310 nm) 이다. Width 는 modulator
이름 `MZMCTE_LULAB_<width>_<length>` 에서 파싱된다.

### 3.2 다이별 파라미터 추출

#### 3.2.1 소광비 (ER_dB)

각 바이어스에서 ALIGN 참조를 차감한 transfer function 은

```
T_dev(lambda) = IL_mzm(lambda) - IL_ref(lambda)
```

소광비는 고정 파장 윈도우 내에서

```
ER_dB = max over biases ( max_lambda T_dev )
      - min over biases ( min_lambda T_dev )
```

윈도우 경계는 밴드당 최소 한 개의 FSR 을 안정적으로 포함하도록
선택된다 (5.2 절). 전체 바이어스 sweep 의 best peak 과 deepest null
을 사용하는 정의는 단일 바이어스 ER 보다 자연스럽게 크며, 디바이스의
최대 modulation depth 를 반영한다.

#### 3.2.2 삽입손실 (IL_dB)

`IL_dB` 는 V = -1 V 의 스펙트럼에서 밴드 중심 ± 5 nm 윈도우 안의
peak `IL_mzm` 값이다. XML 의 IL 필드가 이미 device-level 투과율이라고
가정하므로 ALIGN 차감을 수행하지 않는다.

#### 3.2.3 반파전압 (Vpi_V)

`extract_vpi` 모듈은 다음 절차로 동작한다.

1. V = -2 V 스펙트럼에서 `scipy.signal.find_peaks` 로 null 들을
   식별한다. -25 dB 보다 깊은 null 만 retain.
2. FSR 은 deep null 들의 간격 median 으로 한다.
3. 각 deep null 에 대해 모든 역바이어스 점에서 parabolic fit 으로
   null 파장을 정밀화한다.
4. (V, null wavelength) 에 직선 피팅을 적용해 `dlambda/dV` 를 얻는다.
5. 총 시프트가 추적 윈도우의 1.5 배를 초과하는 fit 은 null 추적
   실패로 reject 한다.
6. 살아남은 slope 들에 대해 3*MAD outlier 트림 후 평균을 취한다.
7. Slope filter 적용: `|dlambda/dV| < MIN_SLOPE_PM_PER_V` (기본
   10 pm/V) 이면 측정 망가짐으로 판정하고 `Vpi_V` 를 NaN 으로
   설정한다. slope 값 자체는 진단을 위해 보존한다.

반환 dict 에는 명시적 상태 필드 `vpi_status` 가 포함되어 다운스트림
코드가 `ok / slope_filter / no_nulls / no_slopes / few_biases /
no_sweeps` 를 구분할 수 있다.

#### 3.2.4 선형성 R^2 (R2_dlam_vs_V)

추적된 각 null 의 linear fit R^2 를 계산하고, slope 평균에 사용한
3*MAD 트림 후 살아남은 fit 들의 median R^2 를 보고한다. V_pi 추출의
1차 신호 품질 지표이다.

#### 3.2.5 품질 등급 (quality_grade)

`R2_dlam_vs_V` 와 `vpi_status` 로부터 등급을 부여한다.

| 등급 | 조건 |
|---|---|
| A | R^2 >= 0.99 |
| B | 0.95 <= R^2 < 0.99 |
| C | 0.90 <= R^2 < 0.95 |
| D | 0.50 <= R^2 < 0.90 |
| F | R^2 < 0.50 또는 vpi_status != 'ok' |

A-B 는 고신뢰 신호 추출에 적합한 다이이며, C-D 는 측정 불확실성이
증가하는 영역, F 는 다운스트림 통계에서 제외해야 하는 다이를 의미한다.

#### 3.2.6 Passive 파라미터

측정된 ER 로부터 splitter amplitude 비를 2.5 절의 닫힌 형식으로
계산한다. MZM 구간 propagation loss 는 ER 윈도우 안에서 모든
바이어스의 T_dev 중 최댓값의 음수로 정의한다 (즉 ALIGN 참조 대비
가장 잘 통과하는 시점의 손실).

바이어스별 splitter imbalance 는 `imbalance_V<V>_dB` 컬럼으로
기록되어, 운용점에 따라 splitter 가 거의 일정하게 유지되는지 (이상적
MMI 라면 그러해야 한다) 를 확인할 수 있다.

### 3.3 물리바운드

bound 는 문헌과 HY202103 의 실측 분포 (5.1 절) 로부터 유도한다.

| Metric | 하한 | 상한 | 근거 |
|---|---|---|---|
| ER_dB | 10 | 45 | Witzens 2018 및 실측 margin |
| IL_dB | -15 | -1 | Reed et al. 2010 |
| Vpi_V | 2 | 60 | V_pi*L = 1 ~ 3 V*cm, L = 0.5 ~ 5 mm |

범위를 벗어난 값은 `reason_<metric>` 에 위반 내용 (예
`over: 5478.34 > 80.0`) 으로 기록되며, `is_problematic` 은 세
metric 의 위반 OR 이다.

### 3.4 Robust Z-score

(Wafer, Band) 그룹마다

```
sigma_robust = 1.4826 * MAD( 그룹 값들 )
z'           = (x - median) / sigma_robust
outlier      if |z'| > 3
```

NaN 값을 가진 다이는 outlier 로 강제 플래깅한다. 상수 1.4826 은
정규분포에서 동일 MAD 를 가질 때의 표준편차 값으로, 3 이라는 임계가
표준 3-sigma 와 정합하면서도 outlier 자신에 오염되지 않는 robust
성질을 가진다 (Iglewicz & Hoaglin, 1993).

### 3.5 변동 분해

`decompose_variation.fit_systematic` 은 각 (Wafer, Band) 그룹의
(Row, Col) 에 대한 2차 이변수 다항식을 fit 한다. 자유 파라미터 6 개,
표본 14 개로 자유도가 8 남으므로 안정적인 fit 이 가능하면서 과적합
위험은 통제된다.

분해 그림은 그룹별로 세 패널 (observed, systematic, random) 을
출력하며, fit 패널에 R^2 을 표기한다.

---

## 4. 구현

### 4.1 폴더 구조

```
.
|-- run.py                            진입점
|-- requirements.txt
|-- src/
|   |-- xml_loader.py                 XML 파싱
|   |-- extract_er.py                 다이별 ER
|   |-- extract_il.py                 다이별 IL
|   |-- extract_vpi.py                V_pi, FSR, dlambda/dV, R^2, 등급
|   |-- extract_passive_params.py     커플러 imbalance, MZM loss
|   |-- outlier_detect.py             물리바운드 + Robust Z
|   |-- csv_export.py                 실행 폴더 생성
|   |-- plot_common.py                플롯 공통 헬퍼
|   |-- wafer_map.py                  연속 surface 웨이퍼맵
|   |-- plot_1d.py                    IQR 박스플롯
|   |-- plot_1d_mad.py                MAD 박스플롯
|   |-- zscore_map.py                 Robust Z-score 격자
|   |-- decompose_variation.py        systematic/random 분해
|   |-- analyze_by_date.py            날짜별 분석 (XLSX 포함)
|   |-- investigate.py                9 항목 진단
|   `-- sensitivity_test.py           ER 윈도우/slope filter 민감도
|-- data/                             (현재 비어있음; HY202103 는
|                                      run.py 의 DATA_ROOT 가 외부 참조)
|-- doc/                              방법론 그림, 진단 결과
`-- res/
    |-- csv/                          최신 실행의 추적 CSV/XLSX
    |-- figures/                      최신 실행의 추적 PNG
    `-- <timestamp>/                  매 실행 raw 결과 (gitignored)
```

### 4.2 모듈 역할

`run.py` 는 파이프라인을 조율한다. `xml_loader.find_all_xmls` 로
XML 을 수집 (각 `(Wafer, Row, Col, Band)` 키에 대해 최신 측정만
보존), `multiprocessing.Pool` 로 다이별 추출을 코어 수만큼 병렬화,
집계된 DataFrame 에 `outlier_detect.mark_outliers` 를 적용, 결과를
실행 폴더에 기록한 후 최신 결과를 `res/csv` 와 `res/figures` 로
mirror 하고, `analyze_by_date.export_and_plot` 으로 날짜별 분석을
생성한 뒤, 추적 대상 출력 파일을 GitHub 원격 저장소에 자동 commit
및 push 한다.

플롯 생성은 `concurrent.futures.ProcessPoolExecutor` 로 4 종 플롯 ×
3 metric = 12 개의 병렬 작업으로 fan-out 된다.

### 4.3 실행

```
python3 run.py
```

standalone 진단 스크립트:

```
python3 src/investigate.py        9 항목 진단 보고
python3 src/sensitivity_test.py   ER 윈도우/slope filter 민감도 스윕
```

### 4.4 출력

매 실행마다 두 CSV (와 동반 PNG) 가 생성된다.

`res/csv/data.csv` 는 dedup 된 다이당 한 행을 가진다. 컬럼은
식별자, 1차 metric, V_pi 진단치, splitter 파라미터, outlier 플래그
순으로 묶여 있다.

`res/csv/data_by_date.csv` (와 `data_by_date.xlsx`) 는 dedup 없이
모든 측정 날짜를 보존한다. XLSX 버전은 `reason_<metric>` 이 비어있지
않은 셀과 `is_problematic` 컬럼의 True 셀에 빨간 배경을 적용한다.

매 실행 스냅샷은 `res/<YYYY-MM-DD_HH-MM-SS>/` 에 저장되며
gitignored. 최신 스냅샷만 `res/csv` 와 `res/figures` 로 mirror 되어
git 으로 추적된다.

### 4.5 CSV 컬럼 정의

`data.csv`:

| 컬럼 | 설명 |
|---|---|
| Wafer, Band, Row, Col, Width_nm | 식별자 |
| ER_dB | 소광비 (3.2.1) |
| IL_dB | 삽입손실 (3.2.2) |
| Vpi_V | 반파전압 (3.2.3) |
| FSR_nm | V=-2 V 스펙트럼의 free spectral range |
| dlam_dV_pm_per_V | 역바이어스에서의 null 이동률 평균 |
| R2_dlam_vs_V | Delta_lambda vs V 직선 피팅의 median R^2 |
| quality_grade | R2_dlam_vs_V 와 vpi_status 로부터 결정되는 A-F 등급 |
| vpi_status | ok / slope_filter / no_nulls / no_slopes / few_biases / no_sweeps 중 하나 |
| amplitude_ratio_k, power_split_ratio, imbalance_dB | 도출된 splitter 파라미터 |
| mzm_loss_dB | MZM 구간 propagation loss (양수 dB) |
| imbalance_V<bias>_dB | 바이어스별 splitter imbalance |
| is_outlier_ER_dB, is_outlier_IL_dB, is_outlier_Vpi_V | 물리바운드 또는 Robust Z outlier |
| robust_z_ER_dB, robust_z_IL_dB, robust_z_Vpi_V | Robust Z-score 값 |
| is_trusted | 세 outlier 플래그가 모두 False 인 경우 True |

`data_by_date.csv` 는 위에 `Date` 컬럼을 추가하고, trust 컬럼 대신
`reason_<metric>`, `out_of_bound_<metric>`, `is_problematic` 을
제공한다.

---

## 5. 결과 (HY202103)

### 5.1 물리바운드의 경험적 검증

ER 상한을 처음에 Witzens (2018) Table II 의 30 dB 기준으로 35 dB 로
설정하였으나, 이는 C-band 다이 28 개 전체를 outlier 로 잘못 분류하였다
(실측 median 37 dB). 두 가지 요인이 차이를 만든다.

1. HY202103 의 splitter (MMI) 가 일반 평균보다 균형이 잘 잡혀 있어
   더 높은 ER 이 가능.
2. 본 ER 정의가 모든 바이어스에 걸친 peak-null 차이 (3.2.1) 라서
   단일 바이어스 ER 보다 자연스럽게 큼.

따라서 상한을 45 dB 로 완화하였다. 이 값은 모든 실측 디바이스 데이터를
포함하면서도 2019-05-31 망가진 측정의 비물리적 아티팩트는 여전히
거른다.

### 5.2 ER 윈도우 선택

`src/sensitivity_test.py` 로 윈도우 폭을 sweep 한 C-band 결과는
다음과 같다.

| 윈도우 | Median ER | std | n > 45 |
|---|---|---|---|
| 14 nm | 36.87 | 1.47 | 0 |
| 16 nm | 37.11 | 1.24 | 0 |
| 18 nm | 37.31 | 1.40 | 0 |
| 22 nm | 38.03 | 1.36 | 0 |
| 36 nm | 39.53 | 2.48 | 2 |

16 nm 윈도우가 다이 간 표준편차를 최소화한다 (한 개의 완전 FSR 이
윈도우에 안정적으로 들어와 peak/null 통계가 개선됨). 동시에 밴드
가장자리 영역 (ALIGN 참조가 덜 신뢰성 있는 영역) 을 포함시키지
않아 아티팩트가 발생하지 않는다. 따라서 C-band 윈도우를 [1545,
1561] nm 로 설정하였다. O-band 의 FSR 은 9.8 nm 이므로 14 nm 윈도우가
이미 1.4 FSR 을 포함하고 있어 변경하지 않았다.

### 5.3 진단 결과

2019-05-31 측정 세트 (D23, D24 wafer 의 O-band 다이 28 개) 는 모두
`|dlambda/dV| < 5 pm/V` 를 보였다. 같은 다이를 2019-06-03 에 다시
측정한 데이터는 `|dlambda/dV| > 100 pm/V` 인 정상 값이다. 10 pm/V
의 slope filter 가 두 분포를 정확히 분리한다: 망가진 28 개 모두
`vpi_status = slope_filter` 로 표시되고, 정상 측정에서는 false
positive 가 발생하지 않는다.

dedup 으로 `data.csv` 에는 06-03 의 정상 측정만 남으므로 main
분석에서는 망가진 데이터의 영향이 없다.

### 5.4 선형성 분포

`R2_dlam_vs_V` 의 그룹별 통계 (평균과 최솟값):

| 그룹 | Mean R^2 | Min R^2 |
|---|---|---|
| D08, C-band | 0.987 | 0.974 |
| D08, O-band | 0.989 | 0.978 |
| D07, C-band | 0.981 | 0.965 |
| D23, O-band | 0.977 | 0.932 |
| D24, O-band | 0.958 | 0.862 |

C-band 측정이 가장 깨끗한 선형 depletion 거동을 보인다. D24 의
O-band 다이 중 일부는 R^2 < 0.90 (등급 C-D) 으로 추출 자체는
성공하지만 통합 통계에서는 가중치를 낮춰야 한다.

### 5.5 변동 분해

위치 의존 성분이 metric 분산 중 차지하는 비율 (degree-2 다항식 fit
의 R^2):

| 그룹 | R^2 (ER) | R^2 (IL) | R^2 (V_pi) |
|---|---|---|---|
| D08, O-band | 0.48 | 0.90 | 0.72 |
| D23, O-band | 0.64 | 0.83 | 0.85 |
| D24, O-band | 0.31 | 0.81 | 0.40 |
| D07, C-band | 0.38 | 0.55 | 0.10 |
| D08, C-band | 0.51 | 0.67 | 0.29 |

모든 그룹에서 IL 이 공간적으로 가장 강하게 종속되어 있는데, 이는
두께와 폭의 gradient 가 도파로 손실에 영향을 주는 현상과 일치한다.
V_pi 는 C-band wafer 에서 거의 random 변동이지만 O-band wafer
에서는 강한 위치 의존성을 보이며, 이는 C-band 공정 run 의 phase
shifter 형상이 더 균일하게 제어되었고 O-band run 에서는 측정 가능한
wafer-scale gradient 가 존재함을 시사한다.

---

## 6. 인용

1. R. A. Soref and B. R. Bennett, "Electrooptical effects in
   silicon", *IEEE Journal of Quantum Electronics* **23**, 123-129
   (1987).
2. G. T. Reed et al., "Silicon optical modulators", *Nature
   Photonics* **4**, 518-526 (2010).
3. A. H. Patel et al., "Design, analysis, and transmission system
   performance of a 41 GHz silicon photonic modulator", *Optics
   Express* **23**, 14263 (2015).
4. J. Witzens, "High-speed silicon photonics modulators",
   *Proceedings of the IEEE* **106**, 2158-2182 (2018).
5. S. K. Selvaraja et al., "Process variation in silicon photonic
   devices", *Applied Optics* **52**, 7638 (2013).
6. Y. Xing, J. Dong, U. Khan, and W. Bogaerts, "Capturing the effects
   of spatial process variations in silicon photonic circuits",
   *ACS Photonics* **10**, 928 (2023).
7. P. Xu et al., "Optical and geometric parameter extraction across
   300-mm photonic integrated circuit wafers", *APL Photonics*
   **9**, 016104 (2024).
8. B. Iglewicz and D. Hoaglin, *How to Detect and Handle Outliers*,
   ASQC Quality Press, 1993.
