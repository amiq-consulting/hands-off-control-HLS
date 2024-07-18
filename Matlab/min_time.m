function [k, optimal_u_s, min_point, n_of_searches, options, n_of_qp_solve] = min_time(Ad, Bd, x0, xf, in_tol, sol_type)
    % sol_type -    1 - Quadprog
    %               2 -ADMM

    n_of_qp_solve = 0;

    k = 5; % Initialize k

    k_min = 0;
    n_of_searches = 0;

    n = size(Ad, 1);
    m = size(Bd, 2);
    q = [];
    [~, Q, q, c] = create_params(Ad, Bd, k, x0, xf);

    if sol_type == 1
        [x_opt, ~, ~, outt] = quadprog(Q, q, [], [], [], [], -ones(k * m, 1), ones(k * m, 1));
        n_of_qp_solve = n_of_qp_solve + 1;
        n_of_searches = n_of_searches + outt.iterations;
    else
        lb = -ones(k * m, 1);
        ub = ones(k * m, 1);
        
        [~, Q, q, c] = create_params(Ad, Bd, k, x0, xf);
        x_opt = qp_ADMM_opt(Q, q, lb, ub, eye(k * m));

        n_of_searches = n_of_searches + 500;

        n_of_qp_solve = n_of_qp_solve + 1;
        % x_opt = sign(x_opt);
    end

    violated_final_state = (x_opt' * Q * x_opt)/2 + q'*x_opt + c;

    while (violated_final_state > in_tol)
        k_min = k;
        k = 2 * k;
        display(['Trying k = ', num2str(k)])
        if (k >= 3e3)
            error('Minimum time is too large. K exceeds 5000. Computations become unstable. Aborting procedure')
        end
        [~, Q, q, c] = create_params(Ad, Bd, k, x0, xf);

        if sol_type == 1
            [x_opt, ~, ~, outt] = quadprog(Q, q, [], [], [], [], -ones(k * m, 1), ones(k * m, 1));
            n_of_searches = n_of_searches + outt.iterations;

            n_of_qp_solve = n_of_qp_solve + 1;
        else
            lb = -ones(k * m, 1);
            ub = ones(k * m, 1);
            
            [~, Q, q, c] = create_params(Ad, Bd, k, x0, xf);
            x_opt = qp_ADMM_opt(Q, q, lb, ub, eye(k * m));
            
            n_of_searches = n_of_searches + 500;
            % x_opt = sign(x_opt);

            n_of_qp_solve = n_of_qp_solve + 1;
        end

        violated_final_state = (x_opt' * Q * x_opt)/2 + q'*x_opt + c;
    end

    k_max = k;

    while(1)
        % Search for the first time we find a solution
        k = round((k_max + k_min) / 2);
        [~, Q, q, c] = create_params(Ad, Bd, k, x0, xf);
        
        if sol_type == 1
            [x_opt, ~, ~, outt] = quadprog(Q, q, [], [], [], [], -ones(k * m, 1), ones(k * m, 1));
            n_of_searches = n_of_searches + outt.iterations;

            n_of_qp_solve = n_of_qp_solve + 1;
        else
            lb = -ones(k * m, 1);
            ub = ones(k * m, 1);
            
            [~, Q, q, c] = create_params(Ad, Bd, k, x0, xf);
            x_opt = qp_ADMM_opt(Q, q, lb, ub, eye(k * m));
            
            n_of_searches = n_of_searches + 500;
            % x_opt = sign(x_opt);

            n_of_qp_solve = n_of_qp_solve + 1;
        end      

        violated_final_state = (x_opt' * Q * x_opt)/2 + q'*x_opt + c;
        if violated_final_state <= in_tol
            if (k_max <= k_min + 1)
                break;
            end
            k_max = k;
        else
            k_min = k;
        end
    end

    optimal_u_s = x_opt;
    min_point = violated_final_state;
    if sol_type == 1
        options = outt;
    else
        options = n_of_searches;
    end
end
