function solver = sysid_gauss_newton_patched(e,nlp,V)
  % Helper function for sysid.m
  % Patched version: uses AMD ordering instead of METIS to avoid METIS 4.0.3 bug
  % See https://github.com/coin-or/Ipopt/issues/208
  J = jacobian(e,V);
  H = triu(J'*J);
  sigma = casadi.MX.sym('sigma');

  io = struct;
  io.x = V;
  io.lam_f = sigma;
  io.hess_gamma_x_x = sigma*H;

  opts = struct;
  disp('WARNING: on Windows, JIT may require a special environment cfr https://github.com/casadi/casadi/wiki/FAQ:-how-to-perform-jit-for-function-evaluations-of-my-optimization-problem%3F')
  opts.jit = true;
  opts.compiler='shell';
  opts.jit_options.verbose = true;
  hessLag = casadi.Function('nlp_hess_l',io,{'x','p','lam_f','lam_g'}, {'hess_gamma_x_x'},opts);
  opts.hess_lag = hessLag;
  % Use AMD ordering (0) instead of automatic (7) to avoid METIS 4.0.3 bug
  % ICNTL(7) values: 0=AMD, 5=METIS, 7=automatic
  % opts.ipopt.mumps_pivot_order = 0;
  solver = casadi.nlpsol('solver','ipopt', nlp, opts);
end
