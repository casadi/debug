import faulthandler; faulthandler.enable()
import casadi as ca, numpy as np
x = ca.MX.sym('x', 1); y = ca.MX.sym('y', 1); x_in = ca.vertcat(x,y)
S = ca.nlpsol('S', 'ipopt', {'x':x_in, 'f':x**2+y**2, 'g':x+y-10},{'ipopt.print_level':0,'print_time':False})
z = S(x0 = x_in, lbg = 0.)["x"]
func = ca.Function("func", [x_in], [z])
n_eval = 100
func_map = func.map(n_eval, "thread", 4)
x_in_val = np.random.multivariate_normal(np.zeros(2), np.eye(2), size = n_eval).T
z_val = np.array(func_map(x_in_val))
ok = np.allclose(z_val, 5.0, atol=1e-6)
print("max err", np.abs(z_val-5).max())
print("BUG ABSENT" if ok else "BUG PRESENT")
