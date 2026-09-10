"""What does the Windows loader do with a BARE dll name?

Separates the two questions the adaptor design turns on:
  1. LoadLibraryExW(bare, LOAD_LIBRARY_SEARCH_*) -- documented to need a full path.
  2. GetModuleHandleExW(bare)  -- does it find a module loaded from another directory?
Run with the full path of a real onnxruntime.dll.
"""
import ctypes
import sys
from ctypes import wintypes

k32 = ctypes.WinDLL("kernel32", use_last_error=True)
k32.LoadLibraryExW.argtypes = [wintypes.LPCWSTR, wintypes.HANDLE, wintypes.DWORD]
k32.LoadLibraryExW.restype = wintypes.HMODULE
k32.GetModuleHandleExW.argtypes = [wintypes.DWORD, wintypes.LPCWSTR,
                                   ctypes.POINTER(wintypes.HMODULE)]
k32.GetModuleFileNameW.argtypes = [wintypes.HMODULE, wintypes.LPWSTR, wintypes.DWORD]

FLAGS = 0x400 | 0x1000 | 0x100  # USER_DIRS | DEFAULT_DIRS | DLL_LOAD_DIR
NAME = "onnxruntime.dll"


def where(h):
    buf = ctypes.create_unicode_buffer(260)
    k32.GetModuleFileNameW(h, buf, 260)
    return buf.value


def report(what, h):
    print("  %-42s -> %s" % (what, where(h) if h else "NULL (err %d)" % ctypes.get_last_error()))


def probe():
    ctypes.set_last_error(0)
    report("LoadLibraryExW(bare, SEARCH_ flags)", k32.LoadLibraryExW(NAME, None, FLAGS))
    h = wintypes.HMODULE()
    ctypes.set_last_error(0)
    ok = k32.GetModuleHandleExW(0, NAME, ctypes.byref(h))
    report("GetModuleHandleExW(bare)", h if ok else None)


if "--after-import-onnxruntime" in sys.argv:
    # Does the pip package leave a DLL a bare CASADI_ONNXRUNTIME_LIB could reuse?
    import onnxruntime
    print("imported onnxruntime %s" % onnxruntime.__version__)
    probe()
    raise SystemExit(0)

print("A. before loading anything")
probe()
if len(sys.argv) > 1:
    print('B. after LoadLibraryW("%s")' % sys.argv[1])
    if not k32.LoadLibraryExW(sys.argv[1], None, FLAGS):
        print("  preload FAILED err %d" % ctypes.get_last_error())
        sys.exit(1)
    probe()
