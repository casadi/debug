import faulthandler; faulthandler.enable()
import casadi as ca
print("casadi", ca.__version__)
bug = False
try:
    x = ca.MX.sym("x"); y = ca.MX.sym("y"); z = ca.MX.sym("z")
    F = ca.nlpsol("F", "ipopt", {"x": ca.vertcat(x, y, z), "f": x**2 + 100*z**2, "g": z + (1-x)**2 - y}, {"print_time": 0, "ipopt.print_level": 0})
    r = F(x0=[2.5, 3.0, 0.75], ubg=0, lbg=0)
    print("x =", r["x"])
    if abs(float(r["f"])) > 1e-3: bug = True
except Exception as e:
    print("ERROR:", e); bug = True
print("BUG PRESENT" if bug else "BUG ABSENT")
