function [x_k] = lp_ADMM_opt(Q, q, lb, ub, A, E)
    
    Nmax = 500;
    
    % sigma = 0.1;
    % rho = 0.2;
    % alfa = 1.8;

    sigma = 0.01;
    rho = 3;
    alfa = 1.9;

    n = size(Q, 1);
    m = size(A, 1);

    x_0 = zeros(n, 1);
    z_0 = zeros(m, 1);
    y_0 = zeros(m, 1);
 
    x_k = x_0;
    x_k_1 = x_k;

    z_k = z_0;
    z_k_1 = z_k;

    y_k = y_0;
    y_k_1 = y_k;

    kkt_mat_small = E' * E + sigma / rho * eye(n / 2) + 2 * eye(n / 2);
    % inv_kkt_mat_small = inv(kkt_mat_small);
    chol_fac = chol(kkt_mat_small);

    for i = 1 : Nmax
        %% Solve KKT system
        b = [sigma * x_k - q;
            z_k - y_k / rho];

        b1 = b(1         : n/2);
        b2 = b(n/2+1     : n);
        b3 = b(n+1       : n+n/2);
        b4 = b(n+n/2+1   : 2*n);
        b5 = b(2*n+1     : 2*n+n/2);
        b6 = b(2*n+n/2+1 : end);

        % u = inv_kkt_mat_small * (E' * b6 + b3 + b4 + b1 / rho);
        b_vec = (E' * b6 + b3 + b4 + b1 / rho);
        % u = chol_fac \ b_vec;
        % u = chol_fac \ u;
        u = linsolve(chol_fac', b_vec);
        u = linsolve(chol_fac, u);
        t = (rho * (-b3 + b4 + b5) + b2) / (sigma + 3*rho);

        v1 = (u - t - b3) * rho;
        v2 = (u + t - b4) * rho;
        v3 = (t - b5) * rho;
        v4 = (E*u - b6) * rho;

        v_k_1 = [v1;v2;v3;v4];
        temp_x_k_1 = [u;t];
        %% Update variables
        temp_z_k_1 = z_k + (v_k_1 - y_k) / rho;
        x_k_1 = alfa * temp_x_k_1 + (1 - alfa) * x_k;
        z_k_1 = alfa * temp_z_k_1 + (1 - alfa) * z_k + y_k / rho;
        z_k_1 = min(ub, max(lb, z_k_1));
        y_k_1 = y_k + rho * (alfa * temp_z_k_1 + (1 - alfa) * z_k - z_k_1);

        %% Shift variables
        x_k = x_k_1;
        z_k = z_k_1;
        y_k = y_k_1;

        % %% Update Rho
        % r_k_p = A * x_k - z_k;
        % r_k_d = Q * x_k + q + A' * y_k;
        % rho = rho * sqrt(max(abs(r_k_p)) / max(abs(r_k_d)));
        % kkt_mat = [Q + sigma * eye(n), A';
        %         A, -eye(n) / rho];
        % inv_kkt_mat = inv(kkt_mat);
    end
end

