# Load commands of every Mach-O + hashes of everything else in an unpacked wheel
import hashlib, os, subprocess, sys
from macho import is_macho

root = sys.argv[1]
KEEP = {"LC_ID_DYLIB", "LC_LOAD_DYLIB", "LC_LOAD_WEAK_DYLIB", "LC_REEXPORT_DYLIB", "LC_RPATH"}
for d, _, fs in sorted(os.walk(root)):
    for f in sorted(fs):
        p = os.path.join(d, f)
        rel = os.path.relpath(p, root)
        if rel.endswith("RECORD"):
            continue
        if not is_macho(p):
            print(rel, "sha", hashlib.sha256(open(p, "rb").read()).hexdigest()[:16])
            continue
        out = subprocess.check_output(["otool", "-l", p], text=True).splitlines()
        cmd = None
        for l in out:
            l = l.strip()
            if l.startswith("cmd "):
                cmd = l.split()[1]
            elif cmd in KEEP and (l.startswith("name ") or l.startswith("path ")):
                print(rel, cmd, l.split()[1])
