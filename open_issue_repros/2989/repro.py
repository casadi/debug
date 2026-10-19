# Issue #2989: intermittent access violation (0xC0000005) in cbc on Windows. Downloads the issue's attachment
# (pickled casadi Function + numpy data, loaded with a whitelisting unpickler) and solves 20 times in a child process.
import sys, os, io, zipfile, tempfile, subprocess, urllib.request, faulthandler
faulthandler.enable()
URL = "https://github.com/casadi/casadi/files/10401994/664_casadi_crash_example.zip"
d = tempfile.mkdtemp()
zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(URL).read())).extractall(d)
pkl = os.path.join(d, "664_casadi_crash_example", "solver_data.pickle")
child = r'''
import sys, pickle, faulthandler; faulthandler.enable()
import casadi as ca, numpy as np
ALLOWED={('numpy.core.multiarray','_reconstruct'),('numpy','dtype'),('numpy.core.multiarray','scalar'),('numpy','ndarray'),('casadi.casadi','Function')}
class U(pickle.Unpickler):
    def find_class(self, m, n):
        assert (m,n) in ALLOWED, (m,n)
        if m=='numpy.core.multiarray' and np.__version__>='2': m='numpy._core.multiarray'
        if m=='casadi.casadi': return ca.Function
        return super().find_class(m,n)
d=U(open(sys.argv[1],'rb')).load()
X=ca.SX.sym('X',*d['nlp_x_shape'])
f,g=d['expand_f_g'](X)
opts=dict(d['nlpsol_options'])
solver=getattr(ca,d['casadi_solver_name'])('nlp',d['my_solver'],{'f':f,'g':g,'x':X},opts)
for i in range(20):
    r=solver(x0=d['x0'],lbx=d['lbx'],ubx=d['ubx'],lbg=ca.veccat(*d['lbg']),ubg=ca.veccat(*d['ubg']))
    print('solve', i, solver.stats()['return_status'], float(r['f']), flush=True)
print('ALLDONE')
'''
r = subprocess.run([sys.executable, "-c", child, pkl], capture_output=True, text=True, timeout=1800)
lines = [l for l in r.stdout.splitlines() if l.startswith("solve") or l == "ALLDONE"]
print("\n".join(lines)); print("exit code", r.returncode); print(r.stderr[-3000:])
print("BUG ABSENT" if r.returncode == 0 and "ALLDONE" in r.stdout else "BUG PRESENT")
