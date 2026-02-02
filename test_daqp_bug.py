import casadi as ca
import numpy as np

rng = np.random.default_rng(0)
n = 100
m = 80
ms = 55
x = ca.SX.sym("x", n)
# --- Hessian: dense symmetric matrix ---
M = rng.standard_normal((n, n))
H = 0.5 * (M @ M.T)                         # symmetrize

# --- Linear term ---
f = 100 * rng.standard_normal((n, 1))

# --- Constraint matrix and bounds ---
A = rng.standard_normal((m, n))
bupper = 20 * rng.random((m, 1))
blower = -20 * rng.random((m, 1))

solver = ca.conic('solver', 'daqp',
                  {'h': ca.DM(H).sparsity(), "a": ca.DM(A).sparsity()},
                #   {'discrete': [1] * ms + [0] * (n-ms)}
                  )

daqp_sol = solver(h=H, g=f, a=A, lbx=[0] * ms + [-10] * (n-ms), ubx=[1] * ms + [10] * (n-ms), lba=blower, uba=bupper)

print(f"Optimal solution: {daqp_sol['x'].full().squeeze()}")
print(f"Optimal objective: {float(daqp_sol['f'])}")
