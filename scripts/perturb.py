# Make an already-delocated tree look like a fresh CI install tree again
import os, subprocess, sys
from macho import machos

root = sys.argv[1]
fake = "/Users/runner/work/casadi/casadi/install/casadi/"
names = set(os.listdir(root))
n_changed = 0
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
    args += ["-add_rpath", "/Users/runner/work/casadi/casadi/build/lib"]
    if "@loader_path" in rpaths:
        args += ["-add_rpath", "@loader_path/"]
    subprocess.check_call(["install_name_tool"] + args + [p], stderr=subprocess.DEVNULL)
    n_changed += 1
print("perturbed", n_changed, "files")
