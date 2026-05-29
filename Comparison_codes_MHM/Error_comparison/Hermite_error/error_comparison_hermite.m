clear all
clc;
n=10000;
    % Get data from each method
    [x1, y1,x11, y11] = error_fourth_her(n);          % First: Red
    [x2, y2, x22, y22] = error_third_her(n);          % Second: Black
    [x4, y4, x44, y44] = error_mhalley_her(n);          % Fourth: Blue

    % Define index for y3 if x not available
   

    % Plot each curve with specific color and marker
    figure;
    semilogy(x1, y1, 'o-', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'FOM');
    hold on;
    semilogy(x2, y2, 's-', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'SOM');
    semilogy(x4, y4, 'd-', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'TOM');

    % Labels and grid
    grid on;
    xlabel('zeros of H_{10000}(x)');
    ylabel('Relative error');
    legend('Location', 'best');
    set(gca, 'FontSize', 12);
%%

figure;
plot(x1, y1, '-', 'Color', 'r', 'LineWidth', 1, 'DisplayName', 'FOM-H');  % Red
hold on;
plot(x2, y2, '-', 'Color', 'k', 'LineWidth', 1, 'DisplayName', 'TOM-H');  % Black
% plot(x3, y3, '-', 'Color', 'g', 'LineWidth', 1, 'DisplayName', 'ASY');  % Green
plot(x4, y4, '-', 'Color', 'b', 'LineWidth', 1, 'DisplayName', 'MHM-H');  % Blue

% Labels and grid
grid on;
xlabel('zeros of H_{10000}(x)');
ylabel('Relative error');
% title(sprintf('Relative Error Comparison of Hermite Zeros (n = %d)', n));
legend('Location', 'best');
set(gca, 'FontSize', 12);

     % Plot each curve with specific color and marker
    figure;
    semilogy(x11, y11, 'o-', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'FOM-H');
    hold on;
    semilogy(x22, y22, 's-', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'TOM-H');
    semilogy(x44, y44, 'd-', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'MHM-H');

    % Labels and grid
    grid on;
    xlabel('Weights of guass hermite quadrature');
    ylabel('Relative error');
    % title(sprintf('Relative Error Comparison of Hermite Zeros (n = %d)', n));
    legend('Location', 'best');
    set(gca, 'FontSize', 12);
%%

figure;
plot(x11, y11, '-', 'Color', 'r', 'LineWidth', 1, 'DisplayName', 'FOM-H');  % Red
hold on;
plot(x22, y22, '-', 'Color', 'k', 'LineWidth', 1, 'DisplayName', 'TOM-H');  % Black
plot(x44, y44, '-', 'Color', 'b', 'LineWidth', 1, 'DisplayName', 'MHM-H');  % Blue

% Labels and grid
grid on;
xlabel('Weights of guass hermite quadrature');
ylabel('Relative error');
% title(sprintf('Relative Error Comparison of Hermite Zeros (n = %d)', n));
legend('Location', 'best');
set(gca, 'FontSize', 12);