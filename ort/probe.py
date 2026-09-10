# Which onnxruntime.dll does casadi's ort plugin bind, and what does GetApi() answer?
import ctypes, os, sys

k32 = ctypes.WinDLL("kernel32", use_last_error=True)
k32.GetModuleHandleW.restype = ctypes.c_void_p
k32.GetModuleHandleW.argtypes = [ctypes.c_wchar_p]
k32.GetModuleFileNameW.argtypes = [ctypes.c_void_p, ctypes.c_wchar_p, ctypes.c_uint32]

class Base(ctypes.Structure):
  _fields_ = [("GetApi", ctypes.CFUNCTYPE(ctypes.c_void_p, ctypes.c_uint32)),
              ("GetVersionString", ctypes.CFUNCTYPE(ctypes.c_char_p))]

def loaded_path(name):
  h = k32.GetModuleHandleW(name)
  if not h:
    return None
  buf = ctypes.create_unicode_buffer(1024)
  k32.GetModuleFileNameW(ctypes.c_void_p(h), buf, 1024)
  return buf.value

def probe(path):
  try:
    lib = ctypes.WinDLL(path)
  except OSError as e:
    print("   LOAD FAILED:", e); return
  try:
    lib.OrtGetApiBase.restype = ctypes.c_void_p
    b = ctypes.cast(lib.OrtGetApiBase(), ctypes.POINTER(Base)).contents
  except Exception as e:
    print("   NO OrtGetApiBase:", e); return
  try:
    ver = b.GetVersionString()
  except Exception as e:
    ver = "<GetVersionString failed: %s>" % e
  print("   version string:", ver)
  print("   GetApi:", {v: ("NULL" if not b.GetApi(v) else "ok") for v in (1, 11, 20, 22, 23, 25)})

mode = sys.argv[1]
print("=" * 70); print("MODE:", mode); print("=" * 70)

if mode == "preimport":
  import onnxruntime as ort
  print("pip onnxruntime", ort.__version__, "imported first")

print("-- candidate DLLs on disk:")
import glob
cands = []
for pat in (os.path.join(os.environ.get("ORT_DIR", ""), "onnxruntime.dll"),
            os.path.join(os.path.dirname(__file__), "..", "ortpip", "onnxruntime", "capi", "onnxruntime.dll")):
  cands += glob.glob(pat)
try:
  import onnxruntime as _o
  cands += glob.glob(os.path.join(os.path.dirname(_o.__file__), "capi", "onnxruntime.dll"))
except Exception:
  pass
for c in dict.fromkeys(cands):
  print(" *", c)
  probe(c)

print("-- before importing casadi, loaded onnxruntime.dll:", loaded_path("onnxruntime.dll"))
import casadi as ca
print("-- casadi", ca.__version__)
print("-- has_onnxbackend('ort'):", ca.has_onnxbackend("ort"))
try:
  ca.load_onnxbackend("ort"); print("-- load_onnxbackend: ok")
except Exception as e:
  print("-- load_onnxbackend FAILED:", str(e)[:300])
p = loaded_path("onnxruntime.dll")
print("-- after loading the plugin, loaded onnxruntime.dll:", p)
if p:
  print("   probing THAT module:"); probe(p)

import onnx
from onnx import helper, TensorProto, numpy_helper
import numpy
g = helper.make_graph([helper.make_node("Add", ["x", "a"], ["y"])], "probe",
      [helper.make_tensor_value_info("x", TensorProto.FLOAT, [1])],
      [helper.make_tensor_value_info("y", TensorProto.FLOAT, [1])],
      [numpy_helper.from_array(numpy.full((1,), 1.0, numpy.float32), "a")])
m = helper.make_model(g, opset_imports=[helper.make_opsetid("", 13)]); m.ir_version = 8
onnx.save(m, "probe.onnx")
try:
  f = ca.GraphBuilder("probe.onnx").create("f")
  print("-- CREATE OK ->", f(0.5))
except Exception as e:
  print("-- CREATE FAILED:", str(e)[:300])
print("-- final loaded onnxruntime.dll:", loaded_path("onnxruntime.dll"))
