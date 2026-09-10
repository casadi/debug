"""Probe the onnxruntime adaptor's CASADI_ONNXRUNTIME_LIB contract.

    python ort/adaptor_probe.py <adaptor> [preload-full-path]

Optionally loads an ONNX Runtime by full path first -- standing in for a host that
already has one in the process -- then asks the adaptor to bind whatever
CASADI_ONNXRUNTIME_LIB names, and reports which runtime it ended up with.
"""
import ctypes
import os
import sys


class OrtApiBase(ctypes.Structure):
    _fields_ = [("GetApi", ctypes.CFUNCTYPE(ctypes.c_void_p, ctypes.c_uint32)),
                ("GetVersionString", ctypes.CFUNCTYPE(ctypes.c_char_p))]


ORT_API_VERSION = 22

adaptor, preload = sys.argv[1], (sys.argv[2] if len(sys.argv) > 2 else None)
print("CASADI_ONNXRUNTIME_LIB=%s" % os.environ.get("CASADI_ONNXRUNTIME_LIB", "<unset>"))

if preload:
    try:
        ctypes.CDLL(preload, mode=getattr(ctypes, "RTLD_GLOBAL", 0))
        print("preload %s -> ok" % preload)
    except OSError as e:
        print("preload %s -> FAILED %s" % (preload, e))
        sys.exit(2)

a = ctypes.CDLL(adaptor)
a.onnxruntime_adaptor_load.argtypes = [ctypes.c_char_p, ctypes.c_uint]
a.onnxruntime_adaptor_load.restype = ctypes.c_int
a.OrtGetApiBase.restype = ctypes.POINTER(OrtApiBase)

buf = ctypes.create_string_buffer(512)
ret = a.onnxruntime_adaptor_load(buf, len(buf))
print("load ret=%d err=%s" % (ret, buf.value.decode(errors="replace")))

base = a.OrtGetApiBase().contents
print("RESULT version=%s api=%s" % (base.GetVersionString().decode(),
                                    "ok" if base.GetApi(ORT_API_VERSION) else "NULL"))
sys.exit(0 if ret == 0 else 1)
