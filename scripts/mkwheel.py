# Same layout as universal_grafter's create_wheel_local.py
import os, sys, zipfile

src, whl = sys.argv[1], sys.argv[2]
info = "dummy-1.0.dist-info/"
meta = {info + "WHEEL": "Wheel-Version: 1.0\nRoot-Is-Purelib: false\n"
                        "Tag: cp39-none-macosx_10_13_x86_64\nTag: cp39-none-macosx_10_13_intel\n",
        info + "METADATA": "Metadata-Version: 2.0\nName: dummy\nVersion: 1.0\n"}
arcs = []
with zipfile.ZipFile(whl, "w", zipfile.ZIP_DEFLATED) as z:
    for d, _, fs in os.walk(src):
        for f in fs:
            p = os.path.join(d, f)
            arc = "dummy/" + os.path.relpath(p, src)
            zi = zipfile.ZipInfo(arc)
            zi.external_attr = (os.stat(p).st_mode & 0xFFFF) << 16
            with open(p, "rb") as fh:
                z.writestr(zi, fh.read(), zipfile.ZIP_DEFLATED)
            arcs.append(arc)
    for arc, data in meta.items():
        z.writestr(arc, data)
        arcs.append(arc)
    arcs.append(info + "RECORD")
    z.writestr(info + "RECORD", "".join(a + ",,\n" for a in arcs))
