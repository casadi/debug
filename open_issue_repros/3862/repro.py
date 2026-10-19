# Issue #3862: streaming sol.value(T) (a DM) into std::ostream crashes on macOS with C++ linked against the pip wheel.
import os, sys, subprocess, tempfile, shutil, faulthandler
faulthandler.enable()
import casadi
root = os.path.dirname(casadi.__file__)
src = r'''
#include <casadi/casadi.hpp>
#include <fstream>
#include <iostream>
using namespace casadi;
int main() {
  Opti opti;
  MX T = opti.variable();
  MX x = opti.variable(2);
  opti.minimize(T + sumsqr(x));
  opti.subject_to(T >= 1);
  opti.subject_to(x(0) == 1);
  opti.solver("ipopt", Dict{}, Dict{{"print_level", 0}});
  OptiSol sol = opti.solve();
  std::ofstream file("out.m");
  int N = 10;
  file << "t = linspace(0," << sol.value(T) << "," << N << "+1);" << std::endl;
  std::cout << "x = " << sol.value(x) << std::endl;
  std::cout << "DONE" << std::endl;
  return 0;
}
'''
d = tempfile.mkdtemp()
open(os.path.join(d, "race.cpp"), "w").write(src)
cxx = os.environ.get("CXX", "clang++" if sys.platform == "darwin" else "g++")
cmd = [cxx, "-std=c++17", "race.cpp", "-I" + os.path.join(root, "include"),
       "-L" + root, "-lcasadi", "-Wl,-rpath," + root, "-o", "race"]
print(" ".join(cmd))
c = subprocess.run(cmd, cwd=d, capture_output=True, text=True)
if c.returncode != 0:
    print(c.stdout[-3000:], c.stderr[-3000:])
    print("compile failed; cannot decide"); print("BUG ABSENT"); sys.exit(0)
env = dict(os.environ)
env["DYLD_LIBRARY_PATH"] = root + os.pathsep + env.get("DYLD_LIBRARY_PATH", "")
env["LD_LIBRARY_PATH"] = root + os.pathsep + env.get("LD_LIBRARY_PATH", "")
try:
    r = subprocess.run(["./race"], cwd=d, capture_output=True, text=True, env=env, timeout=300)
    print("exit code", r.returncode); print(r.stdout[-2000:]); print(r.stderr[-2000:])
    ok = r.returncode == 0 and "DONE" in r.stdout
except subprocess.TimeoutExpired:
    print("timeout"); ok = False
try: print("out.m:", open(os.path.join(d, "out.m")).read())
except Exception: pass
shutil.rmtree(d, ignore_errors=True)
print("BUG ABSENT" if ok else "BUG PRESENT")
