# Issue 3876: Opti('conic') + fatrop crash (reported in MATLAB on Windows). Runs child process to catch hard crashes.
import subprocess, sys, faulthandler
faulthandler.enable()
child = r'''
import faulthandler; faulthandler.enable()
import casadi as ca
for mode in ["conic", "sqp"]:
    opti = ca.Opti("conic") if mode == "conic" else ca.Opti()
    N = 5
    x = opti.variable(2, N+1); u = opti.variable(1, N)
    J = 0
    for k in range(N):
        opti.subject_to(x[:, k+1] == ca.DM([[1, 0.1], [0, 1]]) @ x[:, k] + ca.DM([0, 0.1]) * u[k])
        J += ca.sumsqr(x[:, k]) + u[k]**2
    opti.subject_to(x[:, 0] == ca.DM([1, 0]))
    if mode == "conic":
        opti.subject_to(opti.bounded(-1, u, 1))
        opti.solver("fatrop")
    else:
        opti.solver("sqpmethod", {"qpsol": "fatrop"})
    opti.minimize(J)
    sol = opti.solve()
    print(mode, sol.value(u))
print("CHILD_OK")
'''
r = subprocess.run([sys.executable, "-c", child], capture_output=True, text=True)
print(r.stdout[-2000:]); print(r.stderr[-2000:])
ok = r.returncode == 0 and "CHILD_OK" in r.stdout
print("BUG ABSENT" if ok else "BUG PRESENT")
sys.exit(0)
