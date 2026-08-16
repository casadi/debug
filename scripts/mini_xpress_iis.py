"""Minimal repro attempt for the macOS `Abort trap: 6` in casadi's --run_slow suite.

Xpress' XPRSiisfirst (reached via casadi's `compute_iis` option on an infeasible
problem) is the prime suspect for corrupting the macOS nano malloc zone.  The
abort is *deferred*: it fires at whatever allocates next, so the victim test
looks unrelated (in CI: function.Functiontests.test_ad_weight_sp, ~7 s after the
IIS tests ran).  This script plants the corruption repeatedly and churns
allocations of many sizes to detonate it.

  CASADI_MINI_IIS=1 (default) -> expect SIGABRT (nanov2_guard_corruption_detected)
  CASADI_MINI_IIS=0           -> control, must run to completion
  CASADI_MINI_ROUNDS=N        -> how many LP+MILP IIS rounds (default 50)
  CASADI_MINI_CHURN=N         -> allocations churned per round (default 20000)

Run with MallocNanoZone=0 MallocGuardEdges=1 to fault at the corrupting write
itself instead of at the deferred detection point.
"""
import ctypes
import os
import sys

import casadi as ca

IIS = os.environ.get("CASADI_MINI_IIS", "1") == "1"
ROUNDS = int(os.environ.get("CASADI_MINI_ROUNDS", "50"))
CHURN = int(os.environ.get("CASADI_MINI_CHURN", "20000"))
MODE = os.environ.get("CASADI_MINI_MODE", "both")   # lp | milp | both

_libc = ctypes.CDLL(None)
try:
    _libc.malloc_zone_check.restype = ctypes.c_bool
    _libc.malloc_zone_check.argtypes = [ctypes.c_void_p]
except AttributeError:
    _libc = None


def zone_ok():
    if _libc is None:
        return "n/a"
    return _libc.malloc_zone_check(None)


def iis_lp():
    """Infeasible LP: x >= 5 and -x >= 5.  Bounds free so casadi's pre-check passes."""
    x = ca.MX.sym("x")
    s = ca.qpsol("s", "xpress",
                 {"x": x, "f": x, "g": ca.vertcat(x, -x)},
                 {"error_on_fail": False, "compute_iis": IIS,
                  "xpress": {"OUTPUTLOG": 0}})
    s(x0=0, lbx=-1e10, ubx=1e10, lbg=ca.DM([5, 5]), ubg=ca.DM([1e10, 1e10]))
    return s.stats()


def iis_milp():
    """Infeasible MILP: binary x in [0,1] with x >= 2."""
    x = ca.MX.sym("x")
    s = ca.qpsol("s", "xpress",
                 {"x": x, "f": x, "g": x},
                 {"discrete": [1], "error_on_fail": False, "compute_iis": IIS,
                  "xpress": {"OUTPUTLOG": 0}})
    s(x0=0, lbx=ca.DM([0]), ubx=ca.DM([1]), lbg=2, ubg=1e10)
    return s.stats()


def churn(n):
    """Allocate many casadi objects across a spread of size classes."""
    keep = []
    for i in range(n):
        keep.append(ca.MX.sym("y%d" % i, 1 + (i % 17)))
        keep.append(ca.DM.zeros(1 + (i % 13)))
        keep.append(ca.Sparsity.dense(1 + (i % 7), 1 + (i % 5)))
        if len(keep) > 3000:
            keep = []


print("START compute_iis=%s mode=%s rounds=%d churn=%d zone_ok=%s"
      % (IIS, MODE, ROUNDS, CHURN, zone_ok()), flush=True)

for r in range(ROUNDS):
    n_lp = len(iis_lp().get("iis_rows", [])) if MODE in ("lp", "both") else -1
    n_mip = len(iis_milp().get("iis_rows", [])) if MODE in ("milp", "both") else -1
    if CHURN:
        churn(CHURN)
    print("ROUND %d mode=%s iis_rows_lp=%s iis_rows_mip=%s" % (r, MODE, n_lp, n_mip), flush=True)

print("NO CRASH", flush=True)
sys.exit(0)
