AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
hydrophobic = "ACFILMV"
def calculate_metrics(sequence):
    if not isinstance(sequence, str):
        raise ValueError("서열은 문자열이어야 합니다.")
    cleaned = []
    counts = {}
    for amino_acid in AMINO_ACIDS:
        counts[amino_acid] = 0
    hydrophobic_count = 0
    current_run = 0
    longest_run = 0
    for character in sequence:
        if character.isspace():
            continue
        amino_acid = character.upper()
        if not character.isascii() or amino_acid not in counts:
            raise ValueError(f"표준 아미노산 기호를 확인하세요: {character!r}")
        cleaned.append(amino_acid)
        counts[amino_acid] = counts[amino_acid] + 1
        if amino_acid in hydrophobic:
            hydrophobic_count = hydrophobic_count + 1
            current_run = current_run + 1
            if current_run > longest_run:
                longest_run = current_run
        else:
            current_run = 0
    length = len(cleaned)
    if length == 0:
        raise ValueError("서열이 비어 있습니다.")
    hydrophilic_count = length - hydrophobic_count
    return{"sequence": "".join(cleaned), "length_aa": length, "amino_acid_counts": counts, "hydrophobic_count": hydrophobic_count, "hydrophilic_count": hydrophilic_count, "hydrophobic_percent": 100 * hydrophobic_count / length, "hydrophilic_percent": 100 * hydrophilic_count / length, "max_hydrophobic_run_aa": longest_run}
        