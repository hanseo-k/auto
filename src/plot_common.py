"""플롯 모듈 공통 헬퍼.

`WAFER_BAND_COLOR` 와 `ordered_groups` 는 여러 플롯 모듈에서 공통으로 사용.
중복 정의 방지를 위해 여기 한 곳에 둠.
"""


# Wafer × Band → 색상 (1D 분포, by-date 등 공통 사용)
WAFER_BAND_COLOR = {
    ('D07', 'C'): '#4C72B0', ('D08', 'C'): '#7AA0CB',
    ('D08', 'O'): '#DD8452', ('D23', 'O'): '#E8A87C', ('D24', 'O'): '#F4B999',
}


def ordered_groups(df):
    """C-band 먼저, 그 다음 O-band. 각 밴드 안에서 wafer 이름 정렬."""
    return sorted({(w, b) for w, b in zip(df['Wafer'], df['Band'])},
                  key=lambda x: (x[1] != 'C', x[0]))
