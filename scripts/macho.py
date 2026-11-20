import os

MAGIC = {b"\xfe\xed\xfa\xce", b"\xce\xfa\xed\xfe", b"\xfe\xed\xfa\xcf",
         b"\xcf\xfa\xed\xfe", b"\xca\xfe\xba\xbe", b"\xbe\xba\xfe\xca"}


def is_macho(path):
    try:
        with open(path, "rb") as f:
            return f.read(4) in MAGIC
    except OSError:
        return False


def machos(root):
    for d, _, fs in os.walk(root):
        for f in sorted(fs):
            p = os.path.join(d, f)
            if not os.path.islink(p) and is_macho(p):
                yield p
