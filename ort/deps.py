# Where can a plugin's dependency live and still be found? (libcasadi_conic_highs -> libhighs.dll)
import os, sys, ctypes, shutil
import casadi as ca

mode = sys.argv[1]
casdir = os.path.dirname(ca.__file__)
depdir = os.path.join(os.path.dirname(casdir), "moved_deps")
print("mode:", mode)
print("CASADI_PLUGIN_SEARCH_PATH:", os.environ.get("CASADI_PLUGIN_SEARCH_PATH"))
print("depdir on PATH:", any(os.path.normcase(x.rstrip("\\")) == os.path.normcase(depdir) for x in os.environ.get("PATH","").split(os.pathsep)))
print("libhighs in casadi dir:", os.path.exists(os.path.join(casdir, "libhighs.dll")))
print("libhighs in moved dir :", os.path.exists(os.path.join(depdir, "libhighs.dll")))
try:
  ca.load_conic("highs")
  print("RESULT: load_conic('highs') OK")
except Exception as e:
  print("RESULT: FAILED:", str(e).strip().splitlines()[0][:150])
