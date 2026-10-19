# Issue #3696: "Cannot load shared library 'libcasadi_interpolant_linear.dll'" (Windows, error 126)
import faulthandler; faulthandler.enable()
import casadi as ca
ok = True
try:
    f = ca.interpolant('f', 'linear', [[0., 1., 2.]], [0., 1., 4.])
    v = float(f(1.5)); print("linear interp(1.5) =", v); ok = ok and abs(v - 2.5) < 1e-12
    g = ca.interpolant('g', 'bspline', [[0., 1., 2., 3.]], [0., 1., 4., 9.])
    print("bspline interp(1.5) =", float(g(1.5)))
except Exception as e:
    print("ERROR:", e); ok = False
print("BUG ABSENT" if ok else "BUG PRESENT")
