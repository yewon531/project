from peptide_metrics import calculate_metrics


def test_sequence(sequence, expected_run):
    result = calculate_metrics(sequence)

    assert result["length_aa"] == 8, "서열 길이가 예상값과 다릅니다."
    assert result["amino_acid_counts"]["L"] == 4, "L 개수가 예상값과 다릅니다."
    assert result["amino_acid_counts"]["S"] == 4, "S 개수가 예상값과 다릅니다."
    assert result["hydrophobic_percent"] == 50.0, "소수성 비율이 예상값과 다릅니다."
    assert result["max_hydrophobic_run_aa"] == expected_run, "최대 연속 길이가 예상값과 다릅니다."

    print(f"통과: {sequence} | 소수성 50% | 최대 연속 길이 {expected_run}개")


if __name__ == "__main__":
    test_sequence("LLLLSSSS", 4)
    test_sequence("LSLSLSLS", 1)
