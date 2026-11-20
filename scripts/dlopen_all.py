import ctypes, os, sys
from macho import machos

ok = fail = 0
for p in machos(sys.argv[1]):
    if not p.endswith(".dylib"):
        continue
    try:
        ctypes.CDLL(p)
        print("OK  ", os.path.basename(p)); ok += 1
    except OSError as e:
        print("FAIL", os.path.basename(p), str(e).split("\n")[0][-120:]); fail += 1
print("dlopen ok=%d fail=%d" % (ok, fail), file=sys.stderr)
