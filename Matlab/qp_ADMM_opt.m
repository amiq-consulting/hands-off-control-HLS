function [x_k] = qp_ADMM_opt(Q, q, lb, ub, A)
    
    Nmax = 500;

    % sigma = 0.05;
    % rho = 0.05;
    % alfa = 1.9;

    sigma = 0.01;
    rho = 0.01;
    alfa = 1.9;
    
    % sigma = 1;
    % rho = 1;
    % alfa = 1;

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

    kkt_mat_small = Q + sigma * eye(n) + rho * eye(n);
    inv_kkt_mat_small = inv(kkt_mat_small);

    for i = 1 : Nmax
        %% Solve KKT system
        b1 = sigma * x_k - q;
        b2 = z_k - y_k / rho;

        temp_x_k_1 = inv_kkt_mat_small * (b1 + rho * b2);
        v_k_1 = rho * (temp_x_k_1 - b2);

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

