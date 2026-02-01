import casadi.*

% Stripped down sysid example to isolate crash
% Tests hessLag evaluation before Ipopt

N = 10000;
fs = 610.1;

param_truth = [5.625e-6;2.3e-4;1;4.69];
param_guess = [5;2;1;5];
scale = [1e-6;1e-4;1;1];

y  = MX.sym('y');
dy = MX.sym('dy');
u  = MX.sym('u');

states = [y;dy];
controls = u;

M = MX.sym('x');
c = MX.sym('c');
k = MX.sym('k');
k_NL = MX.sym('k_NL');

params = [M;c;k;k_NL];

rhs = [dy; (u-k_NL*y.^3-k*y-c*dy)/M];

ode = Function('ode',{states,controls,params},{rhs});

N_steps_per_sample = 10;
dt = 1/fs/N_steps_per_sample;

k1 = ode(states,controls,params);
k2 = ode(states+dt/2.0*k1,controls,params);
k3 = ode(states+dt/2.0*k2,controls,params);
k4 = ode(states+dt*k3,controls,params);

states_final = states+dt/6.0*(k1+2*k2+2*k3+k4);

one_step = Function('one_step',{states, controls, params},{states_final});

X = states;
for i=1:N_steps_per_sample
    X = one_step(X, controls, params);
end

one_sample = Function('one_sample',{states, controls, params}, {X});
one_sample = one_sample.expand();

all_samples = one_sample.mapaccum('all_samples', N);

u_data = 0.1*rand(N,1);

x0 = DM([0,0]);
X_measured = all_samples(x0, u_data, repmat(param_truth,1,N));
y_data = X_measured(1,:)';

% Multiple shooting setup (the large problem that crashes)
X = MX.sym('X', 2, N);

res = one_sample.map(N, 'thread', 4);
Xn = res(X, u_data', repmat(params.*scale,1,N));

gaps = Xn(:,1:end-1)-X(:,2:end);
e = y_data-Xn(1,:)';

V = veccat(params, X);

nlp = struct('x',V, 'f',0.5*dot(e,e),'g',vec(gaps));

% Build hessLag (from sysid_gauss_newton)
J = jacobian(e,V);
H = triu(J'*J);
sigma = MX.sym('sigma');

io = struct;
io.x = V;
io.lam_f = sigma;
io.hess_gamma_x_x = sigma*H;

opts = struct;
opts.jit = true;
opts.compiler = 'shell';
opts.jit_options.verbose = true;

disp('Creating JIT-compiled hessLag...');
hessLag = Function('nlp_hess_l', io, {'x','p','lam_f','lam_g'}, {'hess_gamma_x_x'}, opts);

disp('hessLag created successfully');
disp(hessLag);

% Test evaluation BEFORE Ipopt
disp('Testing hessLag evaluation...');

n_x = size(V, 1);

x_test = rand(n_x, 1);

disp(['n_x = ' num2str(n_x)]);

disp('Calling hessLag with keyword args...');
H_result = hessLag('x', x_test, 'lam_f', 1.0);
disp('hessLag evaluation succeeded!');
H_mat = H_result.hess_gamma_x_x;
disp(['Result size: ' num2str(size(H_mat,1)) ' x ' num2str(size(H_mat,2))]);
disp(['Result nnz: ' num2str(nnz(H_mat))]);

disp('Test complete - no crash during hessLag evaluation');
