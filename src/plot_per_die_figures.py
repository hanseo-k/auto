"""다이별 6 그림 생성 — dispatcher.

C-band 와 O-band 의 fit 파라미터가 달라 각각 별도 모듈로 분리됨.
이 스크립트는 두 모듈을 순차 실행한다.

실행:
    python3 src/plot_per_die_figures.py        (C + O 모두)
    python3 src/plot_per_die_figures_C.py      (C-band 만)
    python3 src/plot_per_die_figures_O.py      (O-band 만)
"""
import sys, os
PROGRAM_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROGRAM_ROOT, 'src'))

from plot_per_die_figures_C import main as main_C
from plot_per_die_figures_O import main as main_O


def main():
    main_C()
    print()
    main_O()


if __name__ == '__main__':
    main()
