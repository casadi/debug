import os
def maps():
    return sorted({l.split()[-1] for l in open('/proc/self/maps') if 'libonnxruntime' in l})
import onnx, numpy
from onnx import helper, TensorProto, numpy_helper
a = numpy_helper.from_array(numpy.full((3,), 2.0, numpy.float32), "a")
b = numpy_helper.from_array(numpy.full((3,), 1.0, numpy.float32), "b")
x = helper.make_tensor_value_info("x", TensorProto.FLOAT, [3])
y = helper.make_tensor_value_info("y", TensorProto.FLOAT, [3])
g = helper.make_graph([helper.make_node("Mul", ["x", "a"], ["t"]),
                       helper.make_node("Add", ["t", "b"], ["y"])], "affine", [x], [y], [a, b])
m = helper.make_model(g, opset_imports=[helper.make_opsetid("", 13)]); m.ir_version = 8
onnx.save(m, "/tmp/m.onnx")
print("LD_LIBRARY_PATH =", os.environ.get("LD_LIBRARY_PATH", "<unset>"))
import casadi as ca
print("has_onnx(ort) =", ca.has_onnx("ort"))
try:
    f = ca.GraphBuilder("/tmp/m.onnx").create("f")
    print("EVAL OK", f(ca.DM([0.5, 1.0, -2.0])).T)
except Exception as e:
    print("CREATE FAILED:", str(e).splitlines()[-1])
print("maps:", maps())
