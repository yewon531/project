# 뇌척수액 펩타이드의 소수성 패턴 비교

**프로그래밍기반 데이터 분석 · 2026.11.27 업데이트 · 서열 분석, 손검산 테스트 및 시각화 코드 구현**

이번 측정에서 알츠하이머 치매군(AD)에서만 검출된 펩타이드와 세 집단에서 공통으로 검출된 펩타이드는 소수성 잔기의 비율과 연속 길이에서 어떤 차이를 보일까?

## 현재 결과

| 전체 고유 서열 | AD에서만 검출 | AD·MCI·Control 공통 검출 | 나머지 검출 조합 | 제외 행 |
| ---: | ---: | ---: | ---: | ---: |
| 1,650 | 256 | 395 | 999 | 0 |

집단별로 모은 뇌척수액을 측정한 자료이므로 표의 단위는 고유 펩타이드 서열이다. 같은 기본 서열이 여러 집단에서 나오면 한 행에 검출 집단을 모두 기록한다. AD 단독 검출은 이 실험의 혼합 시료에서 관측된 결과를 뜻한다.

| 결과 파일 | 확인할 내용 |
| --- | --- |
| [peptide_summary.csv](results/peptide_summary.csv) | 한 행당 기본 서열 하나. 집단별 검출 여부(1 또는 0), 검출 조합, 길이, 소수성·친수성 개수와 비율, 소수성 잔기의 최대 연속 길이 |
| [amino_acid_counts.csv](results/amino_acid_counts.csv) | 같은 서열의 아미노산 20종별 개수. 요약표와 연결할 때는 두 표의 `sequence` 열을 사용 |
| [excluded_rows.csv](results/excluded_rows.csv) | 제외된 행의 시트명, 엑셀 행 번호, 원래 값과 사유. 이번 실행에서는 열 이름만 기록 |

## 계산 방법

프로젝트의 이분류 기준에서 소수성 잔기는 `A C F I L M V` 7종이다. 소수성 비율은 이 잔기들의 개수를 전체 길이로 나눈 백분율이다. 최대 연속 길이는 원래 서열을 왼쪽부터 읽으며 소수성 잔기가 몇 개까지 연이어 나오는지 센 값이다.

| 손검산 서열 | 길이 (aa) | 소수성 비율 | 최대 연속 길이 (aa) |
| --- | ---: | ---: | ---: |
| `LLLLSSSS` | 8 | 50% | 4 |
| `LSLSLSLS` | 8 | 50% | 1 |

## 실행하기

원 논문의 **Supplemental Table 7** 엑셀 파일을 내려받아 `data/raw/`에 둔다. 입력 시트는 `Control`, `MCI`, `AD`이며, 원본 B열의 이름은 `Annotated Sequence`다. 예를 들어 `[A].QLTRVPR.[S]`에서는 가운데 `QLTRVPR`을 분석한다. [논문 및 보충자료](https://pmc.ncbi.nlm.nih.gov/articles/PMC9240109/)

저장소 최상위 폴더에서 다음 순서로 실행한다.

```bash
python3 -m pip install -r requirements.txt
python3 analyze_peptides.py data/raw/NIHMS1771217-supplement-Supplemental_Table_7.xlsx
python3 count_amino_acids.py
python3 test_peptide_metrics.py
python3 plot_peptides.py
```

| 코드 | 역할 |
| --- | --- |
| [peptide_metrics.py](peptide_metrics.py) | 서열 한 개의 아미노산 조성, 소수성 비율과 최대 연속 길이 계산 |
| [analyze_peptides.py](analyze_peptides.py) | 시트별 미리보기, 기본 서열별 집단 통합, 요약표와 제외 기록 저장 |
| [count_amino_acids.py](count_amino_acids.py) | 요약표의 고유 서열을 읽고 아미노산별 개수를 별도 CSV로 저장 |

| [test_peptide_metrics.py](test_peptide_metrics.py) | 소수성 비율이 같은 두 서열의 최대 연속 길이 손검산 |
| [plot_peptides.py](plot_peptides.py) | 상자그림과 산점도를 results 폴더에 PNG로 저장 |
| [requirements.txt](requirements.txt) | 실행에 필요한 Python 라이브러리 |

### FASTA 입력 및 Codespaces

`python3 analyze_peptides.py example.fasta`로 FASTA도 분석할 수 있다. FASTA에는 집단 정보가 없으므로 검출 조합은 `unknown`으로 기록되며, 집단 비교 그래프는 Excel 입력 결과에 대해서만 실행한다. FASTA와 Excel 분석은 같은 결과 CSV를 덮어쓰므로 필요한 결과를 따로 보관한다.

GitHub의 Code → Codespaces에서 브라우저 환경을 열고 위 설치 및 실행 명령을 사용한다. 원본 엑셀은 별도로 준비한다. 기존 결과 CSV가 있는 경우 그래프 코드부터 실행할 수도 있다.

## 진행 기록

- **11월 23일 · 명세 개정** — FASTA 조성 분석에서 실제 뇌척수액의 집단별 검출 기록을 비교하는 주제로 확장했다. [초기 기획](https://github.com/yewon531/project/blob/a42643c34a2ba36e883521be4ac541529dbed3c9/README.md) · [개정 명세](docs/specification.md)
- **11월 24일 · 1~2단계 완료** — 서열 지표 계산과 엑셀 분석을 구현했다. 원본 보충표로 재실행해 생성한 CSV 3개가 기존 결과 파일과 일치함을 확인했다.
- **11월 27일 · 3단계** — 손검산 테스트 코드 작성 및 두 예제 통과 확인.
- **11월 27일 · 4단계** — 상자그림 및 산점도 생성 코드 추가. FASTA 입력 지원 및 CSV 헤더, 제외 행 번호 오류 수정. 원본 엑셀을 이용한 이번 버전의 재실행은 아직 확인하지 않았다.
