import sys
from openpyxl import load_workbook
from peptide_metrics import calculate_metrics
def preview_excel(path):
    workbook = load_workbook(path, read_only=True)
    for sheet_name in workbook.sheetnames:
        sheet = workbook[sheet_name]
        for row in sheet.iter_rows(min_row=2, max_row=6, min_col=2, max_col=2, values_only = True):
            annotated = row[0]
            if annotated is None:
                continue
            parts = annotated.split(".")
            sequence = parts[1]
            result = calculate_metrics(sequence)
            print(sheet_name, sequence, result["length_aa"], result["hydrophobic_count"], result["max_hydrophobic_run_aa"])
    workbook.close()
if len(sys.argv) == 2:
    preview_excel(sys.argv[1])
else:
    print("사용법: python3 analyze_peptides.py 엑셀 파일경로")

import csv
from pathlib import Path
def analyze_excel(path):
    peptides = {}
    excluded = []
    groups = ["AD", "MCI", "Control"]
    workbook = load_workbook(path, read_only=True)
    try:
        for group in groups:
            if group not in workbook.sheetnames:
                raise ValueError(f"{group} 시트가 없습니다.")

            sheet = workbook[group]
            sheet.reset_dimensions()
            rows = sheet.values
            header = next(rows, ())

            if "Annotated Sequence" not in header:
                raise ValueError(f"{group} 시트에 Annotated Sequence 열이 없습니다.")

            column = header.index("Annotated Sequence")

            # 첫 번째 줄은 열 이름이므로 데이터의 행 번호는 2부터 시작한다.
            for row_number, row in enumerate(rows, start=2):
                if not any(row):
                    continue

                annotated = None
                if column < len(row):
                    annotated = row[column]

                try:
                    if not isinstance(annotated, str):
                        raise ValueError("서열이 비어 있거나 문자열이 아닙니다.")

                    sequence = annotated.strip()
                    if "." in sequence:
                        parts = sequence.split(".")
                        if len(parts) != 3:
                            raise ValueError("서열 표기를 확인하세요.")
                        sequence = parts[1]

                    result = calculate_metrics(sequence)

                except ValueError as error:
                    excluded.append([group, row_number, annotated, str(error)])
                    continue

                # 처음 보는 서열이면 계산값과 빈 집단 목록을 함께 저장한다.
                key = result["sequence"]
                if key not in peptides:
                    peptides[key] = {"result": result, "groups": []}

                # 같은 서열이 다른 집단에서도 나오면 집단 이름을 추가한다.
                detected = peptides[key]["groups"]
                if group not in detected:
                    detected.append(group)
    finally:
        workbook.close()

    folder = Path(__file__).resolve().parent / "results"
    folder.mkdir(exist_ok=True)

    with open(folder / "peptide_summary.csv", "w",
              newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "sequence", "AD", "MCI", "Control", "detection_pattern",
            "length_aa", "hydrophobic_count", "hydrophilic_count",
            "hydrophobic_percent", "hydrophilic_percent",
            "max_hydrophobic_run_aa"
        ])

        for sequence in sorted(peptides):
            result = peptides[sequence]["result"]
            detected = peptides[sequence]["groups"]

            if detected == ["AD"]:
                pattern = "AD_only"
            elif len(detected) == 3:
                pattern = "all_three"
            else:
                pattern = "+".join(detected)

            writer.writerow([
                sequence,
                int("AD" in detected),
                int("MCI" in detected),
                int("Control" in detected),
                pattern,
                result["length_aa"],
                result["hydrophobic_count"],
                result["hydrophilic_count"],
                result["hydrophobic_percent"],
                result["hydrophilic_percent"],
                result["max_hydrophobic_run_aa"]
            ])

    with open(folder / "excluded_rows.csv", "w",
              newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sheet", "excel_row", "original_value", "reason"])
        writer.writerows(excluded)

    print(f"고유 서열 {len(peptides)}개, 제외 행 {len(excluded)}개")
    print("results 폴더에 peptide_summary.csv와 excluded_rows.csv를 저장했습니다.")


if __name__ == "__main__" and len(sys.argv) == 2:
    analyze_excel(sys.argv[1])
