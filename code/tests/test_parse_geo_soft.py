from pathlib import Path

from parse_geo_soft import parse_soft


def test_parse_two_samples(tmp_path: Path) -> None:
    source = tmp_path / "GSE1.samples.soft.txt"
    source.write_text(
        "^SAMPLE = GSM1\n"
        "!Sample_title = first\n"
        "!Sample_characteristics_ch1 = age: old\n"
        "!Sample_characteristics_ch1 = treatment: control\n"
        "^SAMPLE = GSM2\n"
        "!Sample_title = second\n",
        encoding="utf-8",
    )
    rows = parse_soft(source)
    assert [row["accession"] for row in rows] == ["GSM1", "GSM2"]
    assert rows[0]["characteristics"] == "age: old | treatment: control"
    assert rows[1]["title"] == "second"

