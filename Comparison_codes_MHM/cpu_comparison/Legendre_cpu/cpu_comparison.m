clear all
clc;
n_values = linspace(1000000, 1300000, 25);

num_runs = 20;

time_method1 = zeros(size(n_values));
time_method2 = zeros(size(n_values));
time_method4 = zeros(size(n_values));

for j = 1:length(n_values)
    n = round(n_values(j));
    time_method1(j) = cpu_time_fourth_test_leg(n, num_runs);
    time_method2(j) = cpu_time_halley_leg(n, num_runs);
    time_method4(j) = cpu_time_third_leg(n, num_runs);
   
   
   
end

%% Plotting the results
figure;
plot(n_values, time_method1, '-r', 'LineWidth', 2); hold on;
plot(n_values, time_method4, '-k', 'LineWidth', 2);
plot(n_values, time_method2, '-b', 'LineWidth', 2);



xlabel('n', 'FontSize', 14);
ylabel('Average CPU Time (seconds)', 'FontSize', 14);

% Create the legend
h_legend = legend('FOM-L', 'TOM-L', 'MHM-L'); 
set(h_legend, 'Location', 'NorthWest', 'FontSize', 12);

title('Legendre Polynomial', 'FontSize', 16);



