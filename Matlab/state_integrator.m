function [x_evol] = state_integrator(Ad, Bd, u, x0, N, dist)
    n = size(Ad, 1);
    m = size(Bd, 2);
    x_evol = zeros(N, n);
    x_evol(1,:) = x0'; 
    for i = 2 : N
        x_evol(i,:) = (Ad * (x_evol(i-1,:)') + Bd * u(:, i) + dist(:, i))';
    end
end
