clc, clear, close all

%% Declare system dynamics
max_N = 300;

F_max = 4 * 0.35e-3;
m = 1;

kb = F_max / m;
ka = 1;

Ac = [0, 0;
    ka, 0];
Bc = [kb; 0];
Cc = eye(2);
Dc = 0;

n = size(Ac, 1);
m = size(Bc, 2);

Sys_c = ss(Ac, Bc, Cc, Dc);

%% Discretize the system and randomize initial state

ws = 2 * 0.5;
Ts = 2 * pi / ws;

Sys_d = c2d(Sys_c, Ts, 'tustin');

Ad = Sys_d.A;
Bd = Sys_d.B;

% x0 = [0.5, 3]';
x0 = randn(2, 1);
xf = [0; 0];

%% Compute analytical Minimum Time

tic
q0 = x0(2) - xf(2);
w0 = x0(1);

if q0 > -0.5 * w0 * abs(w0) * ka / kb
    t_min = w0 / kb + sqrt(2 * (w0 / kb) ^ 2 + 4 * q0 / (ka * kb));
elseif q0 < -0.5 * w0 * abs(w0) * ka / kb
    t_min = -w0 / kb + sqrt(2 * (w0 / kb) ^ 2 - 4 * q0 / (ka * kb));
else
    t_min = abs(w0 / kb);
end

analytical_n_min = t_min / Ts;

display(['Analytical Minimum Time: ', num2str(analytical_n_min)]);

n_min = ceil(analytical_n_min);

r = 0.5;
% n_L1 = ceil(2 * n_min);
n_L1 = ceil(t_min / (Ts * r));

display(['L1 Time Horizon: ', num2str(n_L1)])

if n_L1 > max_N
    display(['Time horizon exceeds Maximum Time Horizon']);
    return
end

%% Create L1 problem parameters & Solve it

[E, Q, q, ~] = create_params(Ad, Bd, n_L1, x0, xf);

E_max = zeros(n, max_N * m);
E_max(:, 1 : (n_L1 * m)) = E;

F = Ad ^ n_L1 * x0 - xf;
lb = [-2 * ones(max_N * m, 1);
    zeros(max_N * m, 1);
    zeros(max_N * m, 1);
    -F];
ub = [zeros(max_N * m, 1);
    2 * ones(max_N * m, 1);
    ones(max_N * m, 1);
    -F];
A_L1 = [eye(max_N * m), -eye(max_N * m);
            eye(max_N * m), eye(max_N * m);
            zeros(max_N * m), eye(max_N * m);
            E_max, zeros(n, max_N * m)];
Q_L1 = zeros(2 * max_N * m);
q_L1 = [zeros(max_N * m, 1);
        1 * ones(max_N * m, 1);];
x_opt = lp_ADMM_opt(Q_L1, q_L1, lb, ub, A_L1, E_max);

u_opt = zeros(m, n_L1);
for it = 1 : n_L1
    for jt = 1 : m
        u_opt(jt, it) = x_opt((it - 1) * m + jt);
        % u_opt(jt, it) = round(x_opt((it - 1) * m + jt));
    end
end
toc
%% Plot the results

tt = 0:Ts:((n_L1 - 1) * Ts);
[~, ~, x_evol] = lsim(Sys_d, u_opt', tt, x0);

plot_states_input(x_evol, tt, u_opt, tt);

