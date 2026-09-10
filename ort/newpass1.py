# Faithful ctypes replica of the patched casadi_os.cpp pass 1, run against the SHIPPED
# plugins so the algorithm can be judged before any casadi rebuild.
#   old = today's code : bare name + DLL_LOAD_DIR, one cookie at a time (removed each turn)
#   new = the patch    : fully qualified candidate, cumulative prefix of registered dirs
import ctypes, os, sys, shutil
from ctypes import wintypes

k32 = ctypes.WinDLL("kernel32", use_last_error=True)
k32.AddDllDirectory.restype = ctypes.c_void_p
k32.AddDllDirectory.argtypes = [ctypes.c_wchar_p]
k32.RemoveDllDirectory.argtypes = [ctypes.c_void_p]
k32.LoadLibraryExW.restype = ctypes.c_void_p
k32.LoadLibraryExW.argtypes = [ctypes.c_wchar_p, ctypes.c_void_p, wintypes.DWORD]
k32.GetModuleHandleW.restype = ctypes.c_void_p
k32.GetModuleHandleW.argtypes = [ctypes.c_wchar_p]
k32.GetModuleFileNameW.argtypes = [ctypes.c_void_p, ctypes.c_wchar_p, wintypes.DWORD]

USER_DIRS, DEFAULT_DIRS, DLL_LOAD_DIR = 0x400, 0x1000, 0x100
FLAGS = USER_DIRS | DEFAULT_DIRS | DLL_LOAD_DIR

def bound(name):
  h = k32.GetModuleHandleW(name)
  if not h: return None
  b = ctypes.create_unicode_buffer(1024)
  k32.GetModuleFileNameW(ctypes.c_void_p(h), b, 1024)
  return b.value

def open_old(lib, search_paths):
  for sp in search_paths:
    if not sp: continue
    c = k32.AddDllDirectory(os.path.abspath(sp))
    h = k32.LoadLibraryExW(lib, None, FLAGS)          # BARE name
    err = ctypes.get_last_error()
    if c: k32.RemoveDllDirectory(c)
    if h: return h, sp, 0
  return None, None, err

def open_new(lib, search_paths):
  cookies, handle, resultpath, err = [], None, None, 0
  for k in range(len(search_paths)):
    if not search_paths[k]: continue
    c = k32.AddDllDirectory(os.path.abspath(search_paths[k]))
    if c: cookies.append(c)
    for i in range(k + 1):
      if not search_paths[i]: continue
      cand = os.path.abspath(os.path.join(search_paths[i], lib))   # FULLY QUALIFIED
      handle = k32.LoadLibraryExW(cand, None, FLAGS)
      err = ctypes.get_last_error()
      if handle:
        resultpath = search_paths[i]; break
    if handle: break
  for c in cookies: k32.RemoveDllDirectory(c)
  return handle, resultpath, err

scenario = sys.argv[1]
algo = sys.argv[2]
import casadi as ca
casdir = os.path.dirname(ca.__file__)
ortdir = os.environ["ORT_DIR"]
root = os.path.dirname(casdir)
depdir = os.path.join(root, "moved_deps")
overdir = os.path.join(root, "override")

if scenario == "ort":
  lib, dep, paths = "libcasadi_onnx_ort.dll", "onnxruntime.dll", [ortdir, casdir]
elif scenario == "highs-dep-moved":
  os.makedirs(depdir, exist_ok=True)
  if os.path.exists(os.path.join(casdir, "libhighs.dll")):
    shutil.move(os.path.join(casdir, "libhighs.dll"), os.path.join(depdir, "libhighs.dll"))
  lib, dep, paths = "libcasadi_conic_highs.dll", "libhighs.dll", [depdir, casdir]
elif scenario == "highs-dep-moved-nopath":
  lib, dep, paths = "libcasadi_conic_highs.dll", "libhighs.dll", [casdir]
elif scenario == "clean-override":
  # a self-contained override dir: plugin AND its dependency
  os.makedirs(overdir, exist_ok=True)
  for f in ("libcasadi_conic_highs.dll", "libhighs.dll"):
    for src in (os.path.join(casdir, f), os.path.join(depdir, f)):
      if os.path.exists(src): shutil.copy(src, os.path.join(overdir, f)); break
  lib, dep, paths = "libcasadi_conic_highs.dll", "libhighs.dll", [overdir, casdir]
elif scenario == "control":
  lib, dep, paths = "libcasadi_conic_highs.dll", "libhighs.dll", [casdir]

print("scenario=%s algo=%s" % (scenario, algo))
print("  search_paths:")
for p in paths: print("    -", p)
h, rp, err = (open_new if algo == "new" else open_old)(lib, paths)
print("  load %s -> %s%s" % (lib, "OK" if h else "FAILED", "" if h else " err=%d" % err))
print("  plugin from :", rp)
print("  dependency %s bound from: %s" % (dep, bound(dep)))
