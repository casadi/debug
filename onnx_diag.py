import sys, os, subprocess

MODEL = "/tmp/m.onnx"

def build_model():
    import onnx, numpy
    from onnx import helper, TensorProto, numpy_helper
    a = numpy_helper.from_array(numpy.full((3,), 2.0, numpy.float32), "a")
    b = numpy_helper.from_array(numpy.full((3,), 1.0, numpy.float32), "b")
    x = helper.make_tensor_value_info("x", TensorProto.FLOAT, [3])
    y = helper.make_tensor_value_info("y", TensorProto.FLOAT, [3])
    g = helper.make_graph([helper.make_node("Mul", ["x", "a"], ["t"]),
                           helper.make_node("Add", ["t", "b"], ["y"])],
                          "affine", [x], [y], [a, b])
    m = helper.make_model(g, opset_imports=[helper.make_opsetid("", 13)])
    m.ir_version = 8
    onnx.save(m, MODEL)

def maps():
    return sorted({l.split()[-1] for l in open('/proc/self/maps') if 'libonnxruntime' in l})

def core():
    import casadi as ca
    print("  has_onnx(ort)=", ca.has_onnx("ort"))
    f = ca.GraphBuilder(MODEL).create("f")
    out = f(ca.DM([0.5, 1.0, -2.0]))
    print("  EVAL OK", out.T)
    print("  maps:", maps())

def scenario(name):
    if name == "s1":                       # CI order: nothing pre-loaded
        core()
    elif name == "s2":                     # import onnxruntime (no session) first
        import onnxruntime
        print("  imported onnxruntime; maps now:", maps())
        core()
    elif name == "s3":                     # real InferenceSession first (forces native load)
        import onnxruntime as ort
        ort.InferenceSession(MODEL)
        print("  made InferenceSession; maps now:", maps())
        core()
    elif name == "s4":                     # ctypes preload the ACTUAL pip libonnxruntime RTLD_GLOBAL
        import ctypes, onnxruntime, glob
        capi = os.path.join(os.path.dirname(onnxruntime.__file__), "capi")
        cands = sorted(glob.glob(os.path.join(capi, "libonnxruntime.so*")))
        print("  pip libonnxruntime candidates:", cands)
        ctypes.CDLL(cands[0], mode=ctypes.RTLD_GLOBAL)
        print("  ctypes-preloaded", cands[0], "; maps now:", maps())
        core()

if __name__ == "__main__":
    if len(sys.argv) == 2:
        scenario(sys.argv[1])
    else:
        build_model()
        print("model written")
        for s in ["s1", "s2", "s3", "s4"]:
            print("==================== %s ====================" % s)
            subprocess.run([sys.executable, __file__, s])
