
import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from matplotlib import font_manager

HYDROPHOBIC = set("ACFILMV")
VALID = set("ACDEFGHIKLMNPQRSTVWY")


def longest_positions(sequence):
    runs = []
    start = None

    for i, aa in enumerate(sequence):
        if aa in HYDROPHOBIC:
            if start is None:
                start = i
        elif start is not None:
            runs.append((start, i))
            start = None

    if start is not None:
        runs.append((start, len(sequence)))

    longest = max(
        (end - begin for begin, end in runs),
        default=0
    )

    positions = set()

    for begin, end in runs:
        if end - begin == longest and longest > 0:
            positions.update(range(begin, end))

    return longest, positions


def load_data(csv_path):
    groups = defaultdict(list)

    required = {
        "sequence",
        "detection_pattern",
        "hydrophobic_percent",
        "max_hydrophobic_run_aa"
    }

    with open(csv_path, encoding="utf-8-sig", newline="") as file:
        reader = csv.DictReader(file)

        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"CSV 필수 열이 없습니다: {sorted(missing)}")

        for csv_line, row in enumerate(reader, 2):
            seq = (row["sequence"] or "").strip().upper()
            group = (
                row["detection_pattern"] or ""
            ).strip() or "unknown"

            if not seq or any(aa not in VALID for aa in seq):
                raise ValueError(
                    f"CSV {csv_line}행의 서열을 확인하세요: {seq}"
                )

            run, marked = longest_positions(seq)

            percent = (
                100 * sum(aa in HYDROPHOBIC for aa in seq)
                / len(seq)
            )

            if run != int(row["max_hydrophobic_run_aa"]):
                raise ValueError(
                    f"CSV {csv_line}행: 최대 연속 길이가 다릅니다."
                )

            if abs(percent - float(row["hydrophobic_percent"])) > 0.001:
                raise ValueError(
                    f"CSV {csv_line}행: 소수성 비율이 다릅니다."
                )

            groups[group].append({
                "sequence": seq,
                "length": len(seq),
                "percent": percent,
                "run": run,
                "marked": marked
            })

    return groups


def make_png(group, batch, page_no,
             total_pages, first_no, total_count, output):

    cell_w = 22
    cell_h = 29
    row_h = 37
    seq_x = 590

    max_len = max(len(r["sequence"]) for r in batch)

    width = max(1400, seq_x + max_len * cell_w + 36)
    height = 170 + len(batch) * row_h + 40

    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)

    font_path = font_manager.findfont("DejaVu Sans Mono")
    regular = ImageFont.truetype(font_path, 18)
    small = ImageFont.truetype(font_path, 17)
    title_font = ImageFont.truetype(font_path, 27)

    draw.text(
        (24, 18),
        f"{group} | {total_count} peptides | page {page_no}/{total_pages}",
        font=title_font,
        fill="#111827"
    )

    draw.text(
        (24, 66),
        "Red: hydrophobic (A,C,F,I,L,M,V) | Dark red: longest run",
        font=small,
        fill="#475569"
    )

    draw.text((24, 116), "No.", font=regular, fill="#111827")
    draw.text((116, 116), "Length", font=regular, fill="#111827")
    draw.text((228, 116), "Hydro %", font=regular, fill="#111827")
    draw.text((386, 116), "Max run", font=regular, fill="#111827")
    draw.text(
        (seq_x, 116),
        "Peptide sequence (N to C)",
        font=regular,
        fill="#111827"
    )

    draw.line(
        (24, 150, width - 24, 150),
        fill="#CBD5E1",
        width=2
    )

    for j, row in enumerate(batch):
        y = 161 + j * row_h

        draw.text(
            (24, y + 4), str(first_no + j),
            font=regular, fill="#334155"
        )
        draw.text(
            (116, y + 4), str(row["length"]),
            font=regular, fill="#334155"
        )
        draw.text(
            (228, y + 4), f'{row["percent"]:.1f}%',
            font=regular, fill="#334155"
        )
        draw.text(
            (386, y + 4), f'{row["run"]} aa',
            font=regular, fill="#334155"
        )

        for i, aa in enumerate(row["sequence"]):
            x = seq_x + i * cell_w

            if i in row["marked"]:
                bg = "#FCA5A5"
                fg = "#7F1D1D"

            elif aa in HYDROPHOBIC:
                bg = "#FFE4E6"
                fg = "#B91C1C"

            else:
                bg = "#F1F5F9"
                fg = "#334155"

            draw.rectangle(
                (x, y, x + cell_w - 2, y + cell_h),
                fill=bg
            )

            draw.text(
                (x + 3, y + 3),
                aa,
                font=regular,
                fill=fg
            )

    img.save(output)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--groups",
        nargs="+",
        help="선택 집단만 시각화"
    )

    parser.add_argument(
        "--per-page",
        type=int,
        default=40,
        help="PNG 한 장당 펩타이드 수"
    )

    args = parser.parse_args()

    if not 1 <= args.per_page <= 80:
        parser.error("--per-page는 1부터 80까지 가능합니다.")

    base = Path(__file__).resolve().parent
    input_csv = base / "results" / "peptide_summary.csv"

    if not input_csv.exists():
        raise SystemExit(
            "먼저 analyze_peptides.py를 실행해 주세요."
        )

    groups = load_data(input_csv)
    selected = args.groups or sorted(groups)

    out = base / "results" / "sequence_visualization_png"

    total_png = 0
    total_seqs = 0

    for group in selected:

        if group not in groups:
            print(f"집단을 찾을 수 없습니다: {group}")
            continue

        rows = sorted(
            groups[group],
            key=lambda r: (
                -r["run"],
                -r["percent"],
                r["sequence"]
            )
        )

        safe_name = re.sub(
            r"[^A-Za-z0-9_-]", "_", group
        )

        folder = out / safe_name
        folder.mkdir(parents=True, exist_ok=True)

        pages = (
            len(rows) + args.per_page - 1
        ) // args.per_page

        for page in range(pages):
            start = page * args.per_page
            batch = rows[start:start + args.per_page]

            output = folder / f"page_{page + 1:03d}.png"

            make_png(
                group,
                batch,
                page + 1,
                pages,
                start + 1,
                len(rows),
                output
            )

            total_png += 1

        print(
            f"{group}: {len(rows)}개 서열, PNG {pages}장"
        )

        total_seqs += len(rows)

    print(f"완료: {total_seqs}개 서열, PNG {total_png}장")
    print(f"저장 위치: {out}")


if __name__ == "__main__":
    main()
