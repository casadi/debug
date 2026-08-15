"""Minimal repro of the macOS `Abort trap: 6` in casadi's --run_slow suite.

Xpress' XPRSiisfirst (reached via casadi's `compute_iis` option on an infeasible
LP) corrupts the macOS nano malloc zone.  The abort is *deferred*: it fires at
whatever allocates next, so the victim test looks unrelated (in CI it was
function.Functiontests.test_ad_weight_sp / test_blazing_spline).  This script
plants the corruption and then churns small allocations to detonate it.

  CASADI_MINI_IIS=1 (default) -> expect SIGABRT (nanov2_guard_corruption_detected)
  CASADI_MINI_IIS=0           -> control, must run to completion
"""
import os
import sys

import casadi as ca

IIS = os.environ.get("CASADI_MINI_IIS", "1") == "1"
CHURN = int(os.environ.get("CASADI_MINI_CHURN", "200000"))

x = ca.MX.sym("x")
# Infeasible: x >= 5 and -x >= 5.  Bounds stay free so casadi's pre-check passes.
solver = ca.qpsol("s", "xpress",
                  {"x": x, "f": x, "g": ca.vertcat(x, -x)},
                  {"error_on_fail": False,
                   "compute_iis": IIS,
                   "xpress": {"OUTPUTLOG": 0}})
solver(x0=0, lbx=-1e10, ubx=1e10, lbg=ca.DM([5, 5]), ubg=ca.DM([1e10, 1e10]))
print("SOLVED compute_iis=%s stats=%s" % (IIS, sorted(solver.stats())), flush=True)

print("CHURN %d" % CHURN, flush=True)
keep = []
for i in range(CHURN):
    keep.append(ca.MX.sym("y%d" % i))
    if len(keep) > 1000:
        keep = []
    if i % 20000 == 0:
        print("  churn %d" % i, flush=True)

print("NO CRASH", flush=True)
sys.exit(0)
