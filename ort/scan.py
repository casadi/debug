# List every onnxruntime*.dll reachable from PATH and ask each what API versions it serves.
import ctypes, os, glob, sys

class Base(ctypes.Structure):
  _fields_ = [("GetApi", ctypes.CFUNCTYPE(ctypes.c_void_p, ctypes.c_uint32)),
              ("GetVersionString", ctypes.CFUNCTYPE(ctypes.c_char_p))]

def probe(path):
  try:
    lib = ctypes.WinDLL(path)
    lib.OrtGetApiBase.restype = ctypes.c_void_p
    b = ctypes.cast(lib.OrtGetApiBase(), ctypes.POINTER(Base)).contents
    print("      version=%s  GetApi(22)=%s" % (b.GetVersionString(),
          "NULL <-- CANNOT SERVE casadi" if not b.GetApi(22) else "ok"))
  except Exception as e:
    print("      probe failed:", str(e)[:150])

seen = set()
for i, d in enumerate(os.environ.get("PATH", "").split(os.pathsep)):
  if not d or not os.path.isdir(d):
    continue
  for f in glob.glob(os.path.join(d, "onnxruntime*.dll")):
    if f.lower() in seen:
      continue
    seen.add(f.lower())
    print("[PATH #%d] %s (%d bytes)" % (i, f, os.path.getsize(f)))
    probe(f)
print("total onnxruntime*.dll on PATH:", len(seen))
