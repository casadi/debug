# Issue #1288: "Fatal Python error: PyThreadState_Get: no current thread" on `import casadi` (macOS, Anaconda)
import subprocess, sys
code = "import faulthandler; faulthandler.enable(); import casadi; x=casadi.SX.sym('x'); print(casadi.Function('f',[x],[x**2])(3)); print('IMPORT_OK', casadi.__version__)"
p = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
print(p.stdout); print(p.stderr[-2000:])
ok = p.returncode == 0 and "IMPORT_OK" in p.stdout and "PyThreadState_Get" not in p.stderr
print("BUG ABSENT" if ok else "BUG PRESENT")
