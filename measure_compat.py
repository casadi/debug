# Max mismatch per function of serialize_3.5.5 references vs this platform
import os, sys, math
from collections import defaultdict
import casadi as ca
d = sys.argv[1]
names = os.listdir(d)
inputs, outputs = defaultdict(list), defaultdict(list)
for n in names:
    (inputs if "_in" in n else outputs)[n.split("_")[0]].append(n)
res = []
for fun in [n for n in names if n.endswith(".casadi")]:
    try:
        f = ca.Function.load(os.path.join(d, fun))
    except Exception:
        continue
    key = fun.split(".")[0]
    worst = 0.0
    for i_n, o_n in zip(sorted(inputs[key]), sorted(outputs[key])):
        try:
            outp = f.call(f.generate_in(os.path.join(d, i_n)))
        except Exception:
            continue
        ref = f.generate_out(os.path.join(d, o_n))
        for name, o, r in zip(f.name_out(), outp, ref):
            if "SuperscsInterface" in str(f) and name.startswith("lam_"): continue
            o, r = ca.DM(o), ca.DM(r)
            if o.nnz() == 0 or o.shape != r.shape: continue
            diff = float(ca.mmax(ca.fabs(ca.densify(o) - ca.densify(r))))
            scale = max(1.0, float(ca.mmax(ca.fabs(ca.densify(r)))))
            worst = max(worst, diff / scale)
    res.append((worst, f.class_name(), key))
res.sort(reverse=True)
for w, c, k in res[:15]:
    print("%.3e  digits_ok=%s  %s %s" % (w, "inf" if w == 0 else int(-math.log10(2 * w)), c, k))
print("functions:", len(res))
