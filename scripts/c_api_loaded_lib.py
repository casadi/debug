# The C API of libcasadi through ctypes on macOS, from a wheel: which libcasadi file to open.
# The wheel carries libcasadi.3.7.dylib and libcasadi.dylib as copies (no symlinks); opening the
# copy python did not load gives a second libcasadi with its own types.
import casadi as ca, ctypes, os, sys, tempfile
import numpy as np

def loaded_libcasadi():
  dyld = ctypes.CDLL("/usr/lib/libSystem.B.dylib")
  dyld._dyld_get_image_name.restype = ctypes.c_char_p
  for i in range(dyld._dyld_image_count()):
    path = dyld._dyld_get_image_name(i).decode()
    if os.path.basename(path).startswith("libcasadi."): return path

class CStats(ctypes.Structure):
  _fields_ = [("p", ctypes.POINTER(ctypes.c_ubyte)), ("cap", ctypes.c_longlong),
              ("needed", ctypes.c_longlong), ("reserve", ctypes.c_void_p), ("data", ctypes.c_void_p)]

def run(path):
  lib = ctypes.CDLL(path)
  lib.casadi_c_stats_sink.restype = CStats
  lib.casadi_c_stats_sink.argtypes = [ctypes.POINTER(ctypes.c_ubyte), ctypes.c_longlong]
  x = ca.SX.sym("x"); f = ca.Function("f", [x], [ca.sin(x)])
  y = ca.MX.sym("y", 1, 4); g = ca.Function("g", [y], [f.map(4, "thread")(y) * 2])
  d = tempfile.mkdtemp(); fn = os.path.join(d, "g.casadi"); g.save(fn)
  ret = lib.casadi_c_push_file(fn.encode())
  if ret != 0: return "push_file -> %d" % ret
  try:
    gid = lib.casadi_c_id(b"g")
    ci = ctypes.c_longlong; sz = [ci() for _ in range(4)]
    lib.casadi_c_work_id(gid, *[ctypes.byref(v) for v in sz])
    P = ctypes.POINTER(ctypes.c_double)
    yin = np.arange(4.); out = np.zeros(4)
    arg = (P * max(sz[0].value, 1))(yin.ctypes.data_as(P)); res = (P * max(sz[1].value, 1))(out.ctypes.data_as(P))
    iw = (ci * max(sz[2].value, 1))(); w = (ctypes.c_double * max(sz[3].value, 1))()
    buf = (ctypes.c_ubyte * 4096)()
    sink = lib.casadi_c_stats_sink(ctypes.cast(buf, ctypes.POINTER(ctypes.c_ubyte)), 4096)
    mem = lib.casadi_c_checkout_id(gid)
    flag = lib.casadi_c_eval_with_stats_id(gid, arg, res, iw, w, mem, ctypes.byref(sink), ctypes.c_longlong(-1))
    lib.casadi_c_release_id(gid, mem)
    S = ca.StatsRecorder.from_bytes(bytes(buf)[:sink.needed])
    ok = flag == 0 and np.allclose(out, 2 * np.sin(yin)) and len(S.to_native()) == 1
    return "flag %d, %d bytes, ok=%s" % (flag, sink.needed, ok)
  finally:
    lib.casadi_c_pop()

print("casadi", ca.__version__, sys.platform)
libdir = ca.GlobalOptions.getCasadiPath()
print("files:", sorted(n for n in os.listdir(libdir) if n.startswith("libcasadi.")))
loaded = loaded_libcasadi()
print("loaded by python:", loaded)
print("NEW (loaded image):", run(loaded))
print("OLD (libcasadi.dylib):", run(os.path.join(libdir, "libcasadi.dylib")))
