import csv
from pathlib import Path
from peptide_metrics import AMINO_ACIDS, calculate_metrics


def save_counts(folder):
    input_path = folder / "peptide_summary.csv"
    output_path = folder / "amino_acid_counts.csv"

    if not input_path.is_file():
        raise SystemExit("먼저 analyze_peptides.py를 실행해 주세요.")

    with open(input_path, encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)

        if reader.fieldnames is None or "sequence" not in reader.fieldnames:
            raise ValueError("요약 CSV에 sequence 열이 없습니다.")

        with open(output_path, "w", encoding="utf-8-sig", newline="") as output:
            writer = csv.writer(output)

            # 첫 줄에 sequence, count_A, count_C 등의 열 이름을 쓴다.
            header = ["sequence"]
            for amino_acid in AMINO_ACIDS:
                header.append("count_" + amino_acid)
            writer.writerow(header)

            total = 0
            for row in reader:
                sequence = row["sequence"]
                result = calculate_metrics(sequence)
                counts = result["amino_acid_counts"]

                # 열 이름을 쓴 순서대로 각 아미노산의 개수를 넣는다.
                output_row = [sequence]
                for amino_acid in AMINO_ACIDS:
                    output_row.append(counts[amino_acid])

                writer.writerow(output_row)
                total += 1

    print(f"서열 {total}개의 아미노산 개수를 저장했습니다.")
    print(output_path)


if __name__ == "__main__":
    folder = Path(__file__).resolve().parent / "results"
    save_counts(folder)

