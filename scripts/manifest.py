# sha256 + mode of every file under a directory, sorted by relative path
import hashlib, os, stat, sys

root = sys.argv[1]
for d, _, fs in sorted(os.walk(root)):
    for f in sorted(fs):
        p = os.path.join(d, f)
        rel = os.path.relpath(p, root)
        if os.path.islink(p):
            print(rel, "symlink", os.readlink(p))
            continue
        h = hashlib.sha256(open(p, "rb").read()).hexdigest()
        print(rel, oct(stat.S_IMODE(os.stat(p).st_mode)), h)
