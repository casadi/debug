"""Which ONNX Runtime modules are actually loaded in this process?

    python ort/loaded_modules.py [--import-onnxruntime]

Decides whether the pip `onnxruntime` package loads the standalone runtime it ships
(so a bare CASADI_ONNXRUNTIME_LIB could reuse it) or only its statically-linked
pybind extension (so it could not).
"""
import ctypes
import sys


def loaded():
    if sys.platform == "win32":
        from ctypes import wintypes
        psapi, k32 = ctypes.WinDLL("psapi"), ctypes.WinDLL("kernel32")
        # HANDLE is pointer-sized: without these the -1 pseudo-handle is passed as a
        # 32-bit int and EnumProcessModules quietly returns nothing.
        psapi.EnumProcessModules.argtypes = [wintypes.HANDLE,
                                             ctypes.POINTER(wintypes.HMODULE),
                                             wintypes.DWORD,
                                             ctypes.POINTER(wintypes.DWORD)]
        k32.GetModuleFileNameW.argtypes = [wintypes.HMODULE, wintypes.LPWSTR, wintypes.DWORD]
        n = wintypes.DWORD()
        arr = (wintypes.HMODULE * 2048)()
        proc = wintypes.HANDLE(-1)  # GetCurrentProcess(), without a lossy round-trip
        if not psapi.EnumProcessModules(proc, arr, ctypes.sizeof(arr), ctypes.byref(n)):
            raise OSError("EnumProcessModules failed: %d" % ctypes.get_last_error())
        out = []
        for i in range(n.value // ctypes.sizeof(wintypes.HMODULE)):
            buf = ctypes.create_unicode_buffer(260)
            k32.GetModuleFileNameW(arr[i], buf, 260)
            out.append(buf.value)
        return out
    if sys.platform == "darwin":
        libc = ctypes.CDLL(None)
        libc._dyld_get_image_name.restype = ctypes.c_char_p
        return [libc._dyld_get_image_name(i).decode()
                for i in range(libc._dyld_image_count())]
    with open("/proc/self/maps") as f:
        return sorted({ln.split()[-1] for ln in f if ln.split()[-1].startswith("/")})


before = [m for m in loaded() if "onnxruntime" in m.lower()]
print("before import: %s" % (before or "none"))
if "--import-onnxruntime" in sys.argv:
    import onnxruntime
    print("imported onnxruntime %s" % onnxruntime.__version__)
    after = [m for m in loaded() if "onnxruntime" in m.lower()]
    print("after import:")
    for m in after:
        print("   %s" % m)
    # the runtime proper, not its pybind wrapper nor the providers shim
    def is_runtime(m):
        b = m.replace("\\", "/").rsplit("/", 1)[-1].lower()
        return b.startswith(("libonnxruntime.", "onnxruntime.dll"))
    print("VERDICT standalone-runtime-loaded=%s" % bool([m for m in after if is_runtime(m)]))
