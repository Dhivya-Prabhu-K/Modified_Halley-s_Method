clear all
clc;
n_values = linspace(1000000, 1300000,25);
num_runs = 20;

time_method1 = zeros(size(n_values));
time_method2 = zeros(size(n_values));
time_method4 = zeros(size(n_values));

for j = 1:length(n_values)
    n = round(n_values(j));
    time_method1(j) = fourth_her_cpu(n, num_runs);
    time_method2(j) = mhalley_cpu(n, num_runs);
    time_method4(j) = third_her_cpu(n, num_runs);
   
end

%% Plotting the results
figure;
plot(n_values, time_method1, '-r', 'LineWidth', 2); hold on;
plot(n_values, time_method4, '-k', 'LineWidth', 2);
plot(n_values, time_method2, '-b', 'LineWidth', 2);

xlabel('n', 'FontSize', 14);
ylabel('Average CPU Time (seconds)', 'FontSize', 14);

% Create the legend
h_legend = legend('FOM-H', 'TOM-H', 'MHM-H');
set(h_legend, 'Location', 'NorthWest', 'FontSize', 14);

% Move legend a bit downward
pos = get(h_legend, 'Position');   % [left bottom width height]
pos(2) = pos(2) - 0.13;            % decrease Y (move downward)
set(h_legend, 'Position', pos);
title('Hermite Polynomial', 'FontSize', 16);
grid on;
set(gca, 'FontSize', 12);



