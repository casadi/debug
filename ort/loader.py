# Measure the Windows loader precedence for the ort plugin's onnxruntime.dll dependency.
# One combination per process (a bound module stays bound).
import ctypes, os, sys
from ctypes import wintypes

k32 = ctypes.WinDLL("kernel32", use_last_error=True)
k32.AddDllDirectory.restype = ctypes.c_void_p
k32.AddDllDirectory.argtypes = [ctypes.c_wchar_p]
k32.LoadLibraryExW.restype = ctypes.c_void_p
k32.LoadLibraryExW.argtypes = [ctypes.c_wchar_p, ctypes.c_void_p, wintypes.DWORD]
k32.GetModuleHandleW.restype = ctypes.c_void_p
k32.GetModuleHandleW.argtypes = [ctypes.c_wchar_p]
k32.GetModuleFileNameW.argtypes = [ctypes.c_void_p, ctypes.c_wchar_p, wintypes.DWORD]

DLL_LOAD_DIR   = 0x100
APPLICATION_DIR= 0x200
USER_DIRS      = 0x400
SYSTEM32       = 0x800
DEFAULT_DIRS   = 0x1000

def bound():
  h = k32.GetModuleHandleW("onnxruntime.dll")
  if not h: return None
  b = ctypes.create_unicode_buffer(1024)
  k32.GetModuleFileNameW(ctypes.c_void_p(h), b, 1024)
  return b.value

case = sys.argv[1]
ortdir = os.environ["ORT_DIR"]
import casadi                      # brings in libcasadi.dll, not the ort plugin
plugin = os.path.join(os.path.dirname(casadi.__file__), "libcasadi_onnx_ort.dll")

flags = {
  "flags+cookie":   (DLL_LOAD_DIR | USER_DIRS | DEFAULT_DIRS, True),   # casadi pass 1, dir kept
  "flags-nocookie": (DLL_LOAD_DIR | USER_DIRS | DEFAULT_DIRS, False),
  "nosystem32":     (DLL_LOAD_DIR | USER_DIRS | APPLICATION_DIR, True),
  "legacy":         (0, True),                                          # casadi pass 2 shape
}[case]
dwFlags, want_cookie = flags

if want_cookie:
  c = k32.AddDllDirectory(ortdir)
  print("AddDllDirectory(%s) -> %s" % (ortdir, "ok" if c else "FAILED"))
print("case=%s flags=0x%x" % (case, dwFlags))
h = k32.LoadLibraryExW(plugin, None, dwFlags)
print("LoadLibraryExW ->", "ok" if h else "FAILED err=%d" % ctypes.get_last_error())
print("BOUND:", bound())
