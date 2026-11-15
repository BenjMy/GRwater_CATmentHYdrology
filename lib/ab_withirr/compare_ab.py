"""Compare A (rain-ETp in atmbc) and B (rain via irrigation_map, atmbc=-ETp): python compare_ab.py DIR_A DIR_B"""
import sys
import numpy as np

def load(d):
    return np.loadtxt(f"{d}/ab/output/dtcoupling", comments="#")

if __name__ == "__main__":
    A, B = load(sys.argv[1]), load(sys.argv[2])
    assert A.shape == B.shape, "different number of time steps"
    for name, col in (("Atmpot-v", 9), ("Atmact-v", 13)):
        a, b = A[:, col].sum(), B[:, col].sum()
        print(f"{name}: A={a:.1f} B={b:.1f} rel.diff={abs(a-b)/abs(a):.2e}")
