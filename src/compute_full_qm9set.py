from pathlib import Path
from collections import Counter
import time

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist


def parse_xyz(filepath):

    with open(filepath, "r") as f:
        lines = f.readlines()

    n_atoms = int(lines[0].strip())

    metadata = lines[1].split()
    qm9_id = int(metadata[1])

    elements = []
    coordinates = []

    for line in lines[2:2 + n_atoms]:

        parts = line.split()

        element = parts[0]

        x = float(parts[1].replace("*^", "e"))
        y = float(parts[2].replace("*^", "e"))
        z = float(parts[3].replace("*^", "e"))

        elements.append(element)
        coordinates.append([x, y, z])

    coords = np.array(coordinates, dtype=float)

    return qm9_id, elements, coords


def hill_composition(elements):

    counts = Counter(elements)
    formula = []

    if "C" in counts:

        n = counts.pop("C")
        formula.append("C" if n == 1 else f"C{n}")

        if "H" in counts:

            n = counts.pop("H")
            formula.append("H" if n == 1 else f"H{n}")

    for element in sorted(counts):

        n = counts[element]

        formula.append(
            element if n == 1 else f"{element}{n}"
        )

    return "".join(formula)


def load_disjoint_ids(filepath):

    df = pd.read_csv(
        filepath,
        encoding="utf-8-sig"
    )

    return set(
        df["QM9_id"].astype(int)
    )



def compute_srd(coords):

    center = np.mean(coords, axis=0)

    radial_distances = np.linalg.norm(
        coords - center,
        axis=1
    )

    return np.sort(radial_distances)[::-1]


def compute_spd(coords):
    return np.sort(pdist(coords))


def build_table(descriptor, xyz_dir, disjoint_ids):

    if descriptor not in ("SRD", "SPD"):
        raise ValueError(
            "descriptor must be 'SRD' or 'SPD'"
        )

    rows = []

    files = sorted(xyz_dir.glob("*.xyz"))

    print(f"\nStarting {descriptor}")
    print(f"Found {len(files)} XYZ files")
    print(f"Excluding {len(disjoint_ids)} disjoint molecules")

    processed = 0
    excluded = 0
    failed = 0

    for index, filepath in enumerate(files, start=1):

        try:

            qm9_id, elements, coords = parse_xyz(filepath)

            if qm9_id in disjoint_ids:
                excluded += 1
                continue

            n_atoms = len(coords)

            if descriptor == "SRD":
                values = compute_srd(coords)

            else:
                values = compute_spd(coords)

            values = np.round(values, 4)

            row = {
                "QM9_id": qm9_id,
                "num_atoms": n_atoms,
                "composition": hill_composition(elements)
            }

            for i, value in enumerate(values, start=1):
                row[f"{descriptor}{i}"] = value

            rows.append(row)

            processed += 1

        except Exception as e:

            failed += 1

            print(
                f"Failed on {filepath}: {e}"
            )

            continue

        if index % 10000 == 0:

            print(
                f"{descriptor}: "
                f"scanned {index}/{len(files)} files | "
                f"processed {processed} | "
                f"excluded {excluded} | "
                f"failed {failed}"
            )


    if descriptor == "SRD":

        value_cols = [
            f"SRD{i}"
            for i in range(1, 30)
        ]

    else:

        value_cols = [
            f"SPD{i}"
            for i in range(1, 407)
        ]

    ordered_columns = [
        "QM9_id",
        "num_atoms",
        "composition"
    ] + value_cols

    df = pd.DataFrame(rows).reindex(
        columns=ordered_columns
    )

    print(f"\n{descriptor} complete")
    print(f"Processed: {processed}")
    print(f"Excluded:  {excluded}")
    print(f"Failed:    {failed}")
    print(f"Rows:      {len(df)}")

    return df


def main():

    project_dir = Path(__file__).resolve().parent.parent

    data_dir = project_dir / "data"

    xyz_dir = data_dir / "xyz_files"

    disjoint_file = data_dir / "QM9_disjoint_entries.csv"

    disjoint_ids = load_disjoint_ids(
        disjoint_file
    )

    print(
        f"Loaded {len(disjoint_ids)} "
        f"disjoint QM9 IDs"
    )

    # SRD

    start = time.perf_counter()

    df_srd = build_table(
        "SRD",
        xyz_dir,
        disjoint_ids
    )

    srd_time = time.perf_counter() - start

    srd_output = data_dir / "QM9_SRD.csv"

    df_srd.to_csv(
        srd_output,
        index=False,
        float_format="%.4f"
    )

    print(
        f"\nSRD full QM9 runtime: "
        f"{srd_time:.2f} seconds"
    )


    # spd
    start = time.perf_counter()

    df_spd = build_table(
        "SPD",
        xyz_dir,
        disjoint_ids
    )

    spd_time = time.perf_counter() - start

    spd_output = data_dir / "QM9_SPD.csv"

    df_spd.to_csv(
        spd_output,
        index=False,
        float_format="%.4f"
    )

    print(
        f"\nSPD full QM9 runtime: "
        f"{spd_time:.2f} seconds"
    )



    print(f"SRD runtime: {srd_time:.2f} seconds")
    print(f"SPD runtime: {spd_time:.2f} seconds")

    print()
    print(f"SRD output: {srd_output}")
    print(f"SPD output: {spd_output}")

    print()
    print("SRD shape:", df_srd.shape)
    print("SPD shape:", df_spd.shape)

    print()
    print("SRD columns:")
    print(df_srd.columns.tolist())

    print()
    print("SPD columns:")
    print(df_spd.columns.tolist())


if __name__ == "__main__":
    main()
