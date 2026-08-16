"""MWE: on macOS, repeatedly solving an infeasible LP with Xpress + compute_iis
corrupts the process heap.

Crashes within a handful of iterations, as SIGABRT
(`nanov2_guard_corruption_detected` in libsystem_malloc, tripped by whoever
allocates next) or SIGSEGV in an unrelated casadi destructor.  Set
CASADI_MWE_IIS=0 for the control: identical run without compute_iis is clean.

Requires the real FICO Xpress runtime (casadi/commercial_solvers).
"""
import os

import casadi as ca

IIS = os.environ.get("CASADI_MWE_IIS", "1") == "1"
N = int(os.environ.get("CASADI_MWE_N", "50"))

x = ca.MX.sym("x")
for i in range(N):
    # Infeasible: x >= 5 and -x >= 5.  Bounds stay free so casadi's pre-check passes.
    solver = ca.qpsol("s", "xpress", {"x": x, "f": x, "g": ca.vertcat(x, -x)},
                      {"error_on_fail": False, "compute_iis": IIS,
                       "xpress": {"OUTPUTLOG": 0}})
    solver(x0=0, lbx=-1e10, ubx=1e10, lbg=ca.DM([5, 5]), ubg=ca.DM([1e10, 1e10]))
    print("iteration %d ok" % i, flush=True)

print("NO CRASH", flush=True)
