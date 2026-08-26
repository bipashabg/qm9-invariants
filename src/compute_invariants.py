from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd


def get_atom_count(filepath):
    
    with open(filepath, "r") as f:
        return int(f.readline().strip())


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

        formula.append(
            "C" if n == 1 else f"C{n}"
        )

        if "H" in counts:

            n = counts.pop("H")

            formula.append(
                "H" if n == 1 else f"H{n}"
            )

    for element in sorted(counts):

        n = counts[element]

        formula.append(
            element if n == 1 else f"{element}{n}"
        )

    return "".join(formula)


def compute_spd(coords):

    n_atoms = len(coords)

    distances = []

    for i in range(n_atoms):

        for j in range(i + 1, n_atoms):

            distance = np.linalg.norm(
                coords[i] - coords[j]
            )

            distances.append(distance)

    return np.sort(distances)


def compute_srd(coords):

    centroid = coords.mean(axis=0)

    radial_distances = np.linalg.norm(
        coords - centroid,
        axis=1
    )

    return np.sort(radial_distances)[::-1]

def main():

    xyz_dir = Path("xyz_files")

    records_3_srd = []
    records_3_spd = []

    records_4_srd = []
    records_4_spd = []

    files = sorted(xyz_dir.glob("*.xyz"))

    print(f"Found {len(files)} XYZ files")


    for filepath in files:

        n_atoms = get_atom_count(filepath)

        if n_atoms not in (3, 4):
            continue

        qm9_id, elements, coords = parse_xyz(filepath)

        composition = hill_composition(elements)

        spd = compute_spd(coords)
        srd = compute_srd(coords)

        spd_record = {
            "QM9_id": qm9_id,
            "composition": composition
        }

        for i, value in enumerate(spd, start=1):
            spd_record[f"SPD{i}"] = value

        srd_record = {
            "QM9_id": qm9_id,
            "composition": composition
        }

        for i, value in enumerate(srd, start=1):
            srd_record[f"SRD{i}"] = value

        if n_atoms == 3:

            records_3_spd.append(spd_record)
            records_3_srd.append(srd_record)

        elif n_atoms == 4:

            records_4_spd.append(spd_record)
            records_4_srd.append(srd_record)


    df_3_srd = pd.DataFrame(records_3_srd)
    df_3_spd = pd.DataFrame(records_3_spd)

    df_4_srd = pd.DataFrame(records_4_srd)
    df_4_spd = pd.DataFrame(records_4_spd)

    df_3_srd.to_csv(
        "QM9_3atoms_SRD.csv",
        index=False
    )

    df_3_spd.to_csv(
        "QM9_3atoms_SPD.csv",
        index=False
    )

    df_4_srd.to_csv(
        "QM9_4atoms_SRD.csv",
        index=False
    )

    df_4_spd.to_csv(
        "QM9_4atoms_SPD.csv",
        index=False
    )

    print()
    print(f"3-atom molecules: {len(df_3_srd)}")
    print(f"4-atom molecules: {len(df_4_srd)}")

    print()
    print("QM9_3atoms_SRD.csv")
    print("QM9_3atoms_SPD.csv")
    print("QM9_4atoms_SRD.csv")
    print("QM9_4atoms_SPD.csv")

    print()
    print("3-atom SPD columns:")
    print(df_3_spd.columns.tolist())

    print()
    print("4-atom SPD columns:")
    print(df_4_spd.columns.tolist())


if __name__ == "__main__":
    main()
