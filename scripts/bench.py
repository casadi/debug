# Time subprocess calls of otool / install_name_tool under env variants
import os, shutil, subprocess, sys, time

lib = sys.argv[1]
N = int(sys.argv[2])
work = lib + ".bench.dylib"
shutil.copy(lib, work)
dyld = os.environ.get("CI_DYLD", "")
sdk = os.environ.get("CI_SDKROOT", "")
base = {k: v for k, v in os.environ.items() if k not in ("DYLD_LIBRARY_PATH", "SDKROOT")}
variants = {
    "none": {},
    "dyld": {"DYLD_LIBRARY_PATH": dyld},
    "sdkroot": {"SDKROOT": sdk},
    "both(CI)": {"DYLD_LIBRARY_PATH": dyld, "SDKROOT": sdk},
}
cmds = {
    "otool -L": ["otool", "-L", work],
    "otool -l": ["otool", "-l", work],
    "install_name_tool -id": ["install_name_tool", "-id", "@rpath/x.dylib", work],
}
for vname, extra in variants.items():
    env = dict(base, **extra)
    for cname, cmd in cmds.items():
        subprocess.run(cmd, env=env, capture_output=True)  # warm-up
        t = time.perf_counter()
        for _ in range(N):
            subprocess.run(cmd, env=env, capture_output=True)
        dt = (time.perf_counter() - t) / N * 1000
        print("%-10s %-22s %7.1f ms/call" % (vname, cname, dt), flush=True)
