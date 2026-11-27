import csv
import sys
from pathlib import Path
from openpyxl import load_workbook
from peptide_metrics import calculate_metrics

GROUPS = ["AD", "MCI", "Control"]

def add_sequence(peptides, sequence, group=None):
    result = calculate_metrics(sequence)
    key = result["sequence"]
    if key not in peptides:
        peptides[key] = {"result": result, "groups": []}
    detected = peptides[key]["groups"]
    if group is not None and group not in detected:
        detected.append(group)
    return result

def analyze_excel(path):
    peptides = {}
    excluded = []
    workbook = load_workbook(path, read_only=True)
    try:
        for group in GROUPS:
            if group not in workbook.sheetnames:
                raise ValueError(f"{group} 시트가 없습니다.")
            sheet = workbook[group]
            sheet.reset_dimensions()
            rows = sheet.values
            header = next(rows, ())
            if "Annotated Sequence" not in header:
                raise ValueError(f"{group} 시트에 Annotated Sequence 열이 없습니다.")
            column = header.index("Annotated Sequence")
            shown = 0
            # 첫 번째 줄은 열 이름이므로 데이터의 행 번호는 2부터 시작한다.
            for row_number, row in enumerate(rows, start=2):
                if not any(row):
                    continue
                raw = None
                if column < len(row):
                    raw = row[column]
                try:
                    if not isinstance(raw, str):
                        raise ValueError("서열이 비어 있거나 문자열이 아닙니다.")
                    sequence = raw.strip()
                    if "." in sequence:
                        parts = sequence.split(".")
                        if len(parts) != 3:
                            raise ValueError("서열 표기를 확인하세요.")
                        sequence = parts[1]
                    result = add_sequence(peptides, sequence, group)
                except ValueError as error:
                    excluded.append([group, row_number, raw, str(error)])
                    continue
                if shown < 5:
                    print(group, result["sequence"], result["length_aa"], result["hydrophobic_count"], result["max_hydrophobic_run_aa"])
                    shown += 1
    finally:
        workbook.close()
    return peptides, excluded

def save_results(peptides, excluded):
    folder = Path(__file__).resolve().parent / "results"
    folder.mkdir(exist_ok=True)
    fields = ["length_aa", "hydrophobic_count", "hydrophilic_count", "hydrophobic_percent", "hydrophilic_percent", "max_hydrophobic_run_aa"]
    with open(folder / "peptide_summary.csv", "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sequence"] + GROUPS + ["detection_pattern"] + fields)
        for sequence in sorted(peptides):
            result = peptides[sequence]["result"]
            detected = peptides[sequence]["groups"]
            if not detected:
                row = [sequence, "", "", "", "unknown"]
            else:
                pattern = "+".join(detected)
                if detected == ["AD"]:
                    pattern = "AD_only"
                elif len(detected) == 3:
                    pattern = "all_three"
                row = [sequence]
                for group in GROUPS:
                    row.append(int(group in detected))
                row.append(pattern)
            for field in fields:
                row.append(result[field])
            writer.writerow(row)
    with open(folder / "excluded_rows.csv", "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sheet", "excel_row", "original_value", "reason"])
        writer.writerows(excluded)
    print(f"고유 서열 {len(peptides)}개, 제외 행 {len(excluded)}개")
    print("results 폴더에 peptide_summary.csv와 excluded_rows.csv를 저장했습니다.")

def read_fasta(path):
    from Bio import SeqIO
    peptides = {}
    excluded = []
    with open(path, encoding="utf-8-sig") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence = str(record.seq)
            try:
                add_sequence(peptides, sequence)
            except ValueError as error:
                # 기존 제외표의 sheet, excel_row는 FASTA에서는 빈칸으로 둔다.
                excluded.append(["", "", sequence, f"{record.id}: {error}"])
    if not peptides and not excluded:
        raise ValueError("FASTA 기록이 없습니다. >로 시작하는 제목 줄을 확인하세요.")
    return peptides, excluded

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("사용법: python3 analyze_peptides.py 입력파일경로")
    path = Path(sys.argv[1])
    extension = path.suffix.lower()
    if extension == ".xlsx":
        peptides, excluded = analyze_excel(path)
    elif extension in [".fasta", ".fa", ".faa"]:
        peptides, excluded = read_fasta(path)
    else:
        sys.exit("지원하는 파일 형식: .xlsx, .fasta, .fa, .faa")
    save_results(peptides, excluded)
