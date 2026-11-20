# Remove Mach-O files that also exist in the conda compiler env lib dir,
# so delocate must find them via DYLD_LIBRARY_PATH and copy them back in
import os, sys
from macho import machos

root, conda_lib = sys.argv[1], sys.argv[2]
conda = set(os.listdir(conda_lib))
removed = []
for p in list(machos(root)):
    if os.path.basename(p) in conda:
        os.remove(p)
        removed.append(os.path.relpath(p, root))
print("removed", len(removed), sorted(removed))
