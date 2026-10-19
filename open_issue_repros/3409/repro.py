# Issue #3409: shell JIT with compiler 'ccache gcc' + temp_suffix False warns "Failed to remove ./tmp_casadi_compiler_shell.dylib"
import os, sys, shutil, subprocess, tempfile, faulthandler
faulthandler.enable()
if shutil.which("ccache") is None:
    print("ccache not installed (brew install ccache); cannot test"); print("BUG ABSENT"); sys.exit(0)
cc = "gcc" if shutil.which("gcc") else "cc"
child = r'''
import casadi as ca
x=ca.MX.sym('x')
o={'jit':True,'compiler':'shell','jit_cleanup':False,'jit_temp_suffix':False,
   'jit_options':{'compiler_flags':['-O3'],'compiler':'ccache %s','verbose':True,'temp_suffix':False}}
import gc
f=ca.Function('f',[x],[ca.sin(x)**2],o)
print("value", f(0.3))
del f; gc.collect()
print("deleted")
''' % cc
d = tempfile.mkdtemp()
r = subprocess.run([sys.executable, "-c", child], cwd=d, capture_output=True, text=True, timeout=600)
out = r.stdout + r.stderr
print(out[-4000:])
shutil.rmtree(d, ignore_errors=True)
bad = ("Failed to remove" in out) or r.returncode != 0 or "value" not in out
print("BUG PRESENT" if bad else "BUG ABSENT")
