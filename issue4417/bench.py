# issue #4417: Ipopt wall time doubled on Windows 3.8.0 -> 3.8.1 with unchanged nlp_* eval times.
# Dense NLP so MUMPS spends its time in medium-size dense BLAS-3 (where OpenBLAS threads kick in);
# reports wall minus evaluation time = time inside Ipopt+MUMPS+BLAS.
import casadi as ca, numpy as np, os, sys, time, glob, ctypes
n = int(sys.argv[1]) if len(sys.argv) > 1 else 600
reps = int(sys.argv[2]) if len(sys.argv) > 2 else 5
m = n // 2
rng = np.random.default_rng(0)
Q = rng.standard_normal((n, n)); Q = Q @ Q.T / n + np.eye(n)
A = rng.standard_normal((m, n)) / np.sqrt(n)
x = ca.MX.sym('x', n)
f = 0.5 * ca.bilin(ca.DM(Q), x, x) + ca.sum1(x**4) - ca.sum1(x)
g = ca.mtimes(ca.DM(A), x) + x[:m]**2
solver = ca.nlpsol('solver', 'ipopt', {'x': x, 'f': f, 'g': g},
                   {'ipopt.print_level': 0, 'print_time': False, 'ipopt.max_iter': 200})
wall = ev = iters = 0
for r in range(reps):
    t0 = time.perf_counter()
    sol = solver(x0=np.zeros(n), lbg=-1, ubg=1, lbx=-0.5, ubx=0.5)
    wall += time.perf_counter() - t0
    st = solver.stats()
    ev += sum(st['t_wall_' + k] for k in ['nlp_f', 'nlp_g', 'nlp_grad_f', 'nlp_jac_g', 'nlp_hess_l'])
    iters += st['iter_count']
cfg = nthr = '?'
try:
    lib = glob.glob(os.path.join(os.path.dirname(ca.__file__), '*casadi-tp-openblas*'))[0]
    ob = ctypes.CDLL(lib)
    ob.openblas_get_config.restype = ctypes.c_char_p
    cfg = ob.openblas_get_config().decode(); nthr = ob.openblas_get_num_threads()
except Exception as e:
    cfg = repr(e)
print(f"RESULT casadi {ca.__version__} n={n} reps={reps} OPENBLAS_NUM_THREADS={os.environ.get('OPENBLAS_NUM_THREADS','<unset>')} "
      f"iters={iters} status={st['return_status']} wall={wall:.2f}s eval={ev:.2f}s ipopt+linsol={wall-ev:.2f}s "
      f"| {cfg} threads={nthr} cpus={os.cpu_count()}", flush=True)
