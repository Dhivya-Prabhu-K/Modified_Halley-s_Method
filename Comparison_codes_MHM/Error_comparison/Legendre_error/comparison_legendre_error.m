clear all
clc;
n=10000;
    % Get data from each method
    [x1, y1,x11, y11] = compare_leg_fourth(n);          % First: Red
    [x2, y2, x22, y22] = compare_leg_third(n);          % Second: Black
    [x4, y4, x44, y44] = compare_leg_mhalley(n);          % Fourth: Blue

    % Define index for y3 if x not available
   

    % Plot each curve with specific color and marker
    figure;
    semilogy(x1, y1, 'o-', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'FOM');
    hold on;
    semilogy(x2, y2, 's-', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'SOM');
    semilogy(x4, y4, 'd-', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'TOM');

    % Labels and grid
    grid on;
    xlabel('zeros of P_{10000}(x)');
    ylabel('Relative error');
    % title(sprintf('Relative Error Comparison of Hermite Zeros (n = %d)', n));
    legend('Location', 'best');
    set(gca, 'FontSize', 12);
%%

figure;
plot(x1, y1, '-', 'Color', 'r', 'LineWidth', 1, 'DisplayName', 'FOM-L');  % Red
hold on;
plot(x2, y2, '-', 'Color', 'k', 'LineWidth', 1, 'DisplayName', 'TOM-L');  % Black
plot(x4, y4, '-', 'Color', 'b', 'LineWidth', 1, 'DisplayName', 'MHM-L');  % Blue

% Labels and grid
grid on;
xlabel('zeros of H_{10000}(x)');
ylabel('Relative error');
% title(sprintf('Relative Error Comparison of Hermite Zeros (n = %d)', n));
legend('Location', 'best');
set(gca, 'FontSize', 12);

     % Plot each curve with specific color and marker
    figure;
    semilogy(x11, y11, 'o-', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'FOM-L');
    hold on;
    semilogy(x22, y22, 's-', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'TOM-L');
    semilogy(x44, y44, 'd-', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'MHM-L');

    % Labels and grid
    grid on;
    xlabel('Weights of guass hermite quadrature');
    ylabel('Relative error');
    % title(sprintf('Relative Error Comparison of Hermite Zeros (n = %d)', n));
    legend('Location', 'best');
    set(gca, 'FontSize', 12);
%%

figure;
plot(x11, y11, '-', 'Color', 'r', 'LineWidth', 1, 'DisplayName', 'FOM-L');  % Red
hold on;
plot(x22, y22, '-', 'Color', 'k', 'LineWidth', 1, 'DisplayName', 'TOM-L');  % Black
plot(x44, y44, '-', 'Color', 'b', 'LineWidth', 1, 'DisplayName', 'MHM-L');  % Blue

% Labels and grid
grid on;
xlabel('Weights of guass hermite quadrature');
ylabel('Relative error');
% title(sprintf('Relative Error Comparison of Hermite Zeros (n = %d)', n));
legend('Location', 'best');
set(gca, 'FontSize', 12);