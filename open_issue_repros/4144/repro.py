# Issue 4144: fatrop as conic solver (Opti('conic')) / sqpmethod qpsol=fatrop crashes kernel on Windows.
# A crash kills this process; run in a subprocess to classify.
import subprocess, sys
child = r'''
import faulthandler; faulthandler.enable()
from casadi import *
for variant in ['conic','sqp']:
    opti = Opti('conic') if variant=='conic' else Opti()
    x = opti.variable(1)
    opti.minimize((x-5)**2)
    if variant=='conic': opti.solver('fatrop')
    else: opti.solver('sqpmethod',{'qpsol':'fatrop'})
    sol = opti.solve()
    v = float(sol.value(x))
    assert abs(v-5)<1e-6, v
print("CHILD_OK")
'''
r = subprocess.run([sys.executable, "-c", child], capture_output=True, text=True)
print(r.stdout[-2000:]); print(r.stderr[-2000:])
print("returncode", r.returncode)
print("BUG ABSENT" if (r.returncode==0 and "CHILD_OK" in r.stdout) else "BUG PRESENT")
