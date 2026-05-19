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
│   ├── trust_map.py         # Robust Z-score 격자맵
│   ├── make_method_figures.py
│   └── make_method_ppt.py
├── data/               # 입력 데이터 (HY202103은 외부 경로 참조)
├── doc/                # 방법론 그림
└── res/
    ├── csv/                 # 최신 CSV (git 추적)
    ├── figures/             # 최신 PNG (git 추적)
    └── <YYYY-MM-DD_HH-MM-SS>/  # 매 실행 결과 (gitignored)
```

## 실행

```bash
python3 run.py
```

실행 시:
1. `DATA_ROOT` (run.py 안에 하드코딩) 경로의 모든 XML 자동 탐색
2. 멀티코어 병렬로 다이 추출 + 플롯 생성
3. 결과를 `res/<timestamp>/` 에 저장
4. 최신 결과를 `res/csv/` 와 `res/figures/` 로 복사
5. GitHub `main` 브랜치로 자동 push

## 의존성

`requirements.txt` 참고. 주요 패키지:

- numpy, pandas, scipy
- matplotlib
- python-pptx (`make_method_ppt.py` 용)

## 출력 컬럼

`res/csv/data.csv`:

| 컬럼 | 설명 |
|------|------|
| Wafer, Band, Row, Col, Width_nm | 다이 식별자 |
| ER_dB | Extinction Ratio (peak − null, 고정 윈도우) |
| IL_dB | Insertion Loss (V=−1V, ±5nm peak transmission) |
| Vpi_V | V_π (FSR / 2·\|dλ/dV\|) |
| FSR_nm, dlam_dV_pm_per_V | V_π 계산 중간값 |
| is_outlier_* | 항목별 outlier 플래그 |
| robust_z_* | Robust Z-score |
| is_trusted | 세 항목 모두 trusted 인가 |
