clc, clear, close all;

load('state_file.dat');
load('hardware_output_file.dat');

k_total = length(state_file);

Ts = 2 * pi;
r = 0.5
T = (0 : k_total - 1) * Ts;

x_rez = [state_file(:, 2), state_file(:, 1)];

plot_states_input(x_rez, T, hardware_output_file', T);

obt_sp_rate = sum(abs(hardware_output_file)');
str = 'Optained Sparsity rates: ';
str_input = strings(1, 1);
for it = 1 : 1
    str_input(it) = ['Input ', num2str(it), ': ', num2str(obt_sp_rate(it) / k_total)];
end

disp(['Desired sparsity rate: ', num2str(r)]);
disp(str);
disp(str_input);
