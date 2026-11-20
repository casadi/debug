# Make an already-delocated tree look like a fresh CI install tree again
import os, subprocess, sys
from macho import machos

root = sys.argv[1]
fake = "/Users/runner/work/casadi/casadi/install/casadi/"
names = set(os.listdir(root))
stats = {"full": 0, "partial": 0, "untouched": 0}
for p in machos(root):
    out = subprocess.check_output(["otool", "-l", p], text=True).splitlines()
    rpaths, deps = [], []
    for i, l in enumerate(out):
        l = l.strip()
        if l.startswith("path "):
            rpaths.append(l.split()[1])
        if l == "cmd LC_LOAD_DYLIB":
            deps.append(out[i + 2].split()[1])
    args = []
    for dep in deps:
        b = dep.split("/")[-1]
        if dep.startswith("@rpath/") and b in names:
            args += ["-change", dep, fake + b]
    extra = ["-add_rpath", "/Users/runner/work/casadi/casadi/build/lib"]
    if "@loader_path" in rpaths:
        extra += ["-add_rpath", "@loader_path/"]
    # Fall back to fewer edits when the header padding is too small
    for attempt in (args + extra, args, extra):
        if attempt and subprocess.run(["install_name_tool"] + attempt + [p],
                                      capture_output=True).returncode == 0:
            stats[len(attempt) == len(args + extra) and "full" or "partial"] += 1
            break
    else:
        stats["untouched"] += 1
print("perturbed", stats)
