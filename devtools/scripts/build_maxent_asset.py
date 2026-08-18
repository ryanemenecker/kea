"""Rebuild Kea's compressed MaxEntScan data asset from maxentpy text tables."""

import argparse
from pathlib import Path

import numpy as np


ACCEPTOR_TABLE_LENGTHS = (
    4**7,
    4**7,
    4**7,
    4**7,
    4**7,
    4**3,
    4**4,
    4**3,
    4**4,
)


def build_asset(score5_path: Path, score3_path: Path, output_path: Path) -> None:
    donor = np.empty(4**7, dtype=np.float64)
    with score5_path.open() as handle:
        for index, line in enumerate(handle):
            _sequence, value = line.split()
            donor[index] = float(value)

    acceptor = {
        index: np.empty(length, dtype=np.float64)
        for index, length in enumerate(ACCEPTOR_TABLE_LENGTHS)
    }
    with score3_path.open() as handle:
        for line in handle:
            table, index, value = line.split()
            acceptor[int(table)][int(index)] = float(value)

    np.savez_compressed(
        output_path,
        score5=donor,
        **{f"score3_{index}": acceptor[index] for index in range(9)},
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("score5", type=Path, help="maxentpy score5_matrix.txt")
    parser.add_argument("score3", type=Path, help="maxentpy score3_matrix.txt")
    parser.add_argument("output", type=Path, help="output maxent_scan.npz")
    args = parser.parse_args()
    build_asset(args.score5, args.score3, args.output)


if __name__ == "__main__":
    main()

