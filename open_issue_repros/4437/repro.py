# Issue #4437: map(n,'thread') workers must inherit caller's x87 control word.
import sys, ctypes, platform, faulthandler
faulthandler.enable()
import casadi as ca

if platform.machine().lower() not in ("x86_64", "amd64"):
    print("x87 only on x86_64; not applicable"); print("BUG ABSENT"); sys.exit(0)

win = sys.platform.startswith("win")
# fnstcw [arg]; ret   /   fldcw [arg]; ret
code_get = bytes([0xD9, 0x39 if win else 0x3F, 0xC3])
code_set = bytes([0xD9, 0x29 if win else 0x2F, 0xC3])
if win:
    k32 = ctypes.windll.kernel32
    k32.VirtualAlloc.restype = ctypes.c_void_p
    k32.VirtualAlloc.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_ulong, ctypes.c_ulong]
    buf = k32.VirtualAlloc(None, 4096, 0x3000, 0x40)
else:
    libc = ctypes.CDLL(None)
    libc.mmap.restype = ctypes.c_void_p
    libc.mmap.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_long]
    buf = libc.mmap(None, 4096, 7, 0x22, -1, 0)
ctypes.memmove(buf, code_get, 3)
ctypes.memmove(buf + 16, code_set, 3)
proto = ctypes.CFUNCTYPE(None, ctypes.POINTER(ctypes.c_ushort))
get_cw_ = proto(buf); set_cw_ = proto(buf + 16)
def get_cw():
    v = ctypes.c_ushort(0); get_cw_(ctypes.byref(v)); return v.value
def set_cw(x):
    v = ctypes.c_ushort(x); set_cw_(ctypes.byref(v))

class X87CW(ca.Callback):
    def __init__(self):
        ca.Callback.__init__(self); self.construct("x87cw", {})
    def get_n_in(self): return 1
    def get_n_out(self): return 1
    def eval(self, arg):
        return [ca.DM(get_cw())]

print("casadi", ca.__version__, "initial cw: 0x%X" % get_cw())
f = X87CW()
bad = False
for cw in (0x27F, 0x37F):
    set_cw(cw)
    direct = int(f(0))
    r_ser = [int(v) for v in ca.DM(f.map(2, "serial")(ca.DM.zeros(1, 2))).full().ravel()]
    r_thr = [int(v) for v in ca.DM(f.map(2, "thread", 2)(ca.DM.zeros(1, 2))).full().ravel()]
    print("caller 0x%X direct: 0x%X" % (cw, direct), "serial:", ["0x%X" % v for v in r_ser], "thread:", ["0x%X" % v for v in r_thr])
    bad = bad or any(v != direct for v in r_thr)
set_cw(0x27F)
print("BUG PRESENT" if bad else "BUG ABSENT")
