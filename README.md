# QM9 Isometry Invariants (SRD & SPD)

Computes two isometry-invariant descriptors from Vitaliy Kurlin's
*Geometric Data Science*:  Sorted Radial Distances (SRD, Def. 4.4.6) and
Sorted Pairwise Distances (SPD, Def. 4.4.5)- for the 3- and 4-atom
molecules of the QM9 dataset.

Both descriptors are invariant under isometry (translation, rotation,
reflection) because they are built only from distances.

## Steps

- Parses QM9 `.xyz` files (element + x, y, z coordinates per atom).
- Filters to molecules with exactly 3 or 4 atoms.
- Computes SPD (all pairwise distances, increasing order) and SRD
  (distances from the unweighted geometric centroid, decreasing order).
- Writes four CSV files: `QM9_{3,4}atoms_{SRD,SPD}.csv`.

## Testing

`src/test_examples.py` reproduces the trapezium and kite values from
Example 4.4.7 of the book, confirming both descriptors (including the
fact that SPD cannot distinguish the two clouds while SRD can).

## Data

The QM9 dataset is not included in this repository. Download it from
figshare (Ramakrishnan et al., 2014) and extract the `.xyz` files into
`data/xyz_files/`.

## Usage

```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python src/compute_invariants.py
```

## Reference

V. Kurlin, *Geometric Data Science*, https://kurlin.org/Geometric-Data-Science-book.pdf
