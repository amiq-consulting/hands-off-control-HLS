clc, clear, close all

max_N = 500;


m = 1;
Fmax = 4 * 0.35e-3;

ka = 1;
kb = Fmax / m;

Ac = [0, 0;
    1, 0];
Bc = [Fmax / m;
    0];

sys_c = ss(Ac, Bc, eye(2), 0);

ws = 2 * 0.5;
% ws = 10;
Ts = 2 * pi / ws;

sys_d = c2d(sys_c, Ts, 'tustin');

% x0i = [0.5, 3]';
x0i = [-0.5, 3]';
xf = [0, 0]';

n = 2;
m = 1;

Ad = sys_d.A;
Bd = sys_d.B;

% method = 'linprog';
% method = 'admm';
method = 'analytical';
%%

Tmin = 5 * Ts;

r = input('Set sparsity rate: r = ');
s = input('Start? ', "s");

total_n_of_searches = 0;
u_L0 = [];
x_rez = [];
k_total = 0;
T_k = [0];
noise_rez = [];

x0 = x0i;
weights = ones(m, 1);

while s == 'y'    
    for itera = 1 : 10
        
        % Min Time computation
        [n_min, opt_point, min_point, n_of_searches, optiones] = min_time(Ad, Bd, x0, xf, 1e-3,2); % Solving Min Time Problem
        total_n_of_searches = total_n_of_searches + n_of_searches;
    
        % Computing horizon length
        T = max(Tmin, (1 / r) * n_min * Ts); % Compute current time horizon length
        T_k = [T_k, T_k(end) + T];
        k_L1 = ceil(T / Ts);
    
        % Solving L1 control problem

        switch method
            case 'linprog'
                [E, ~, ~, ~] = create_params(Ad, Bd, k_L1, x0, xf);
                E = [E, zeros(n, k_L1*m)];
                F = xf - Ad^k_L1 * x0;
        
                f = [zeros(k_L1*m, 1); repmat(weights, k_L1, 1)];
                A = [eye(k_L1*m), -eye(k_L1*m);
                    -eye(k_L1*m), -eye(k_L1*m)];
                b = [zeros(2*k_L1*m, 1)];
                LB = [-ones(k_L1*m, 1); zeros(k_L1 * m, 1)];
                UB = [ones(k_L1*m, 1); ones(k_L1 * m, 1)];
                optiones2 = optimoptions("linprog", 'OptimalityTolerance', 1e-9, 'MaxIterations', 1e6, 'ConstraintTolerance', 1e-9);
                x_opt_L1 = linprog(f, A, b, E, F, LB, UB, optiones2);
                x_opt_L1 = x_opt_L1(1:k_L1 * m);

                % Save computed L1 commands
                u_opt_L1_temp = zeros(m, k_L1);
                for it = 1 : k_L1
                    for jt = 1 : m
                        % u_opt_L1_temp(jt, it) = round(x_opt_L1((it - 1) * m + jt));
                        % u_opt_L1_temp(jt, it) = x_opt_L1((it - 1) * m + jt);
                        u_opt_L1_temp(jt, it) = round(x_opt_L1((it - 1) * m + jt));
                    end
                end

            case 'admm'
                [E, ~, ~, ~] = create_params(Ad, Bd, k_L1, x0, xf);
                F = Ad ^ k_L1 * x0 - xf;
                lb = [-2 * ones(k_L1 * m, 1);
                    zeros(k_L1 * m, 1);
                    zeros(k_L1 * m, 1);
                    -F];
                ub = [zeros(k_L1 * m, 1);
                    2 * ones(k_L1 * m, 1);
                    ones(k_L1 * m, 1);
                    -F];
                A_L1 = [eye(k_L1 * m), -eye(k_L1 * m);
                            eye(k_L1 * m), eye(k_L1 * m);
                            zeros(k_L1 * m), eye(k_L1 * m);
                            E, zeros(n, k_L1 * m)];
                Q_L1 = zeros(2 * k_L1 * m);
                q_L1 = [zeros(k_L1 * m, 1);
                        repmat(weights, k_L1, 1)];
                x_opt_L1 = lp_ADMM_opt(Q_L1, q_L1, lb, ub, A_L1, E);
                x_opt_L1 = x_opt_L1(1:k_L1 * m);

                % Save computed L1 commands
                u_opt_L1_temp = zeros(m, k_L1);
                for it = 1 : k_L1
                    for jt = 1 : m
                        % u_opt_L1_temp(jt, it) = round(x_opt_L1((it - 1) * m + jt));
                        % u_opt_L1_temp(jt, it) = x_opt_L1((it - 1) * m + jt);
                        u_opt_L1_temp(jt, it) = round(x_opt_L1((it - 1) * m + jt));
                    end
                end

            case 'analytical'
                q0 = x0(2);
                w0 = x0(1);

                Tf = k_L1 * Ts;

                if q0 > -0.5 * w0 * abs(w0) * ka / kb
                    first = -1;
                    t_sw_1 = 0.5 * (Tf + w0 / kb - sqrt((Tf - w0 / kb) ^ 2 - 4 * q0 / kb - 2 * (w0 / kb) ^ 2));
                    t_sw_2 = 0.5 * (Tf + w0 / kb + sqrt((Tf - w0 / kb) ^ 2 - 4 * q0 / kb - 2 * (w0 / kb) ^ 2));
                elseif q0 < -0.5 * w0 * abs(w0) * ka / kb
                    first = 1;
                    t_sw_1 = 0.5 * (Tf - w0 / kb - sqrt((Tf + w0 / kb) ^ 2 + 4 * q0 / kb - 2 * (w0 / kb) ^ 2));
                    t_sw_2 = 0.5 * (Tf - w0 / kb + sqrt((Tf + w0 / kb) ^ 2 + 4 * q0 / kb - 2 * (w0 / kb) ^ 2));
                else
                    t_min = abs(w0 / kb);
                end
                
                sample_sw_1 = ceil(t_sw_1 / Ts);
                sample_sw_2 = floor(t_sw_2 / Ts);
                
                u_opt_L1_temp = zeros(k_L1, 1)';
                u_opt_L1_temp(1 : sample_sw_1) = first;
                u_opt_L1_temp(sample_sw_2 : end) = -first;
                u_opt_L1_temp(end) = -first;

            otherwise
        end

        % Apply commands and noise on the system
        % noise = randn(n, k_L1 + 1) / 200;

        noise = [0.005 * (rand(1, k_L1 + 1) * 2 - 1);
            zeros(1, k_L1 + 1)];

        % [x_evol, ~, ~] = lsim(Sys_d, [u_opt_L1_temp, zeros(m, 1)] + noise, (1:(k_L1 + 1)) * Ts, x0);
        x_evol = state_integrator(Ad, Bd, [u_opt_L1_temp, zeros(m, 1)], x0, k_L1 + 1, noise);

        % Save computed values
        u_L0 = [u_L0, u_opt_L1_temp, zeros(m, 1)]; % L0 Optimal Control
        x_rez = [x_rez; x_evol]; % Computed state trajectory
        k_total = k_total + k_L1 + 1; % Final discrete step
        noise_rez = [noise_rez, noise]; % Save generated noise in order to simulate the system under no control law, but with this noise
    
        x0 = x_evol(end, :)'; % Reinitialize intitial state
    end
    s = input('Continue 10 times? ', "s");
end
    
%%
t = (0 : k_total - 1) * Ts;
plot_states_input(x_rez, t, u_L0, t);

%% Output displaying

obt_sp_rate = sum(abs(u_L0)');
str = 'Optained Sparsity rates: ';
str_input = strings(m, 1);
for it = 1 : m
    str_input(it) = ['Input ', num2str(it), ': ', num2str(obt_sp_rate(it) / k_total)];
end

disp(['Desired sparsity rate: ', num2str(r)]);
disp(str);
disp(str_input);
