function [E, Q, q, c] = create_params(Ad, Bd, N, x0, xf)

    F = Ad^N * x0 - xf;
    temp_mat = Bd;
    E = temp_mat;
    for i = 1 : N - 1
        temp_mat = Ad * temp_mat;
        E = [temp_mat, E];
    end
    
    Q = 2 * E' * E;
    q = 2 * E' * F;
    c = F' * F;
    
end
