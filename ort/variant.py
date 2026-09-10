# Can we steer casadi's ort plugin away from C:\Windows\System32\onnxruntime.dll (1.5.2)?
import ctypes, os, sys

k32 = ctypes.WinDLL("kernel32", use_last_error=True)
k32.GetModuleHandleW.restype = ctypes.c_void_p
k32.GetModuleHandleW.argtypes = [ctypes.c_wchar_p]
k32.GetModuleFileNameW.argtypes = [ctypes.c_void_p, ctypes.c_wchar_p, ctypes.c_uint32]

def loaded_path(name="onnxruntime.dll"):
  h = k32.GetModuleHandleW(name)
  if not h:
    return None
  buf = ctypes.create_unicode_buffer(1024)
  k32.GetModuleFileNameW(ctypes.c_void_p(h), buf, 1024)
  return buf.value

mode = sys.argv[1]
ortdir = os.environ["ORT_DIR"]
print("=" * 60); print("VARIANT:", mode); print("=" * 60)
print("CASADI_PLUGIN_SEARCH_PATH =", os.environ.get("CASADI_PLUGIN_SEARCH_PATH"))

if mode == "adddll":
  os.add_dll_directory(ortdir); print("os.add_dll_directory(%s)" % ortdir)
elif mode == "preload":
  ctypes.WinDLL(os.path.join(ortdir, "onnxruntime.dll")); print("ctypes-preloaded by full path")

import casadi as ca
try:
  ca.load_onnxbackend("ort"); print("load_onnxbackend: ok")
except Exception as e:
  print("load_onnxbackend FAILED:", str(e)[:200])

import numpy, onnx
from onnx import helper, TensorProto, numpy_helper
g = helper.make_graph([helper.make_node("Add", ["x", "a"], ["y"])], "probe",
      [helper.make_tensor_value_info("x", TensorProto.FLOAT, [1])],
      [helper.make_tensor_value_info("y", TensorProto.FLOAT, [1])],
      [numpy_helper.from_array(numpy.full((1,), 1.0, numpy.float32), "a")])
m = helper.make_model(g, opset_imports=[helper.make_opsetid("", 13)]); m.ir_version = 8
onnx.save(m, "probe.onnx")
try:
  f = ca.GraphBuilder("probe.onnx").create("f")
  print("RESULT: CREATE OK ->", f(0.5))
except Exception as e:
  print("RESULT: CREATE FAILED:", str(e).strip().splitlines()[-1][:120])
print("bound onnxruntime.dll:", loaded_path())
