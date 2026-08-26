import numpy as np
from compute_invariants import compute_spd, compute_srd

T = np.array([[2, 1, 0], [-2, 1, 0], [4, -1, 0], [-4, -1, 0]], dtype=float)
K = np.array([[5, 0, 0], [-3, 0, 0], [-1, 2, 0], [-1, -2, 0]], dtype=float)

print("SPD(T):", compute_spd(T))
print("SPD(K):", compute_spd(K))
print("SRD(T):", compute_srd(T))
print("SRD(K):", compute_srd(K))
